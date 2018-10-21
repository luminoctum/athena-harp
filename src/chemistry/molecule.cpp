// C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../parameter_input.hpp"
#include "../misc.hpp"
#include "../math_funcs.hpp"
#include "../globals.hpp"
#include "molecule.hpp"

std::ostream& operator<<(std::ostream &os, Molecule const& mol)
{
  os << "name: " << mol.myname << std::endl
     << "molecular weight: " << mol.mu << " kg/mol" << std::endl
     << "heat capacity (cp): " << mol.m_cp << " J/(mol K)" << std::endl
     << "vaporization heat at triple point: " << mol.m_latent << " kJ/mol" << std::endl
     << "vapor pressure coefficients: " << mol.m_beta << " " << mol.m_gamma << std::endl
     << "standard entropy: " << mol.m_entropy << " J/(mol K)" << std::endl
     << "standard enthalpy: " << mol.m_enthalpy << " kJ/mol" << std::endl;
  return os;
}

Molecule::Molecule(std::string name):
  myname(name), mu(0), tr(0), pr(0), tc(0), pc(0), phase(GAS),
  m_cp(0), m_latent(0), m_entropy(0), m_enthalpy(0), m_gibbs(0),
  m_cliq(0), m_enliq(0), m_csld(0), m_ensld(0),
  m_beta(0), m_gamma(0), m_nshomate(0)
{
  prev = NULL;
  next = NULL;

  for (int i = 0; i < MAXSHOMATE; ++i)
    for (int j = 0; j < NSHOMATE; ++j)
      m_shomate[i][j] = 0.;
  for (int i = 0; i < MAXSHOMATE; ++i)
    m_shomate_sp[i] = 0.;
}

Molecule::Molecule(ParameterInput *pin):
  mu(0), tr(0), pr(0), tc(0), pc(0), phase(GAS),
  m_cp(0), m_latent(0), m_entropy(0), m_enthalpy(0), m_gibbs(0),
  m_cliq(0), m_enliq(0), m_csld(0), m_ensld(0),
  m_beta(0), m_gamma(0), m_nshomate(0)
{
  prev = NULL;
  next = NULL;

  std::stringstream msg;

  std::string gas = pin->GetString("chemistry", "gas");
  std::string cloud = pin->GetOrAddString("chemistry", "cloud", "");
  std::string folder = pin->GetOrAddString("chemistry", "folder", "");

  std::vector<std::string> agas, acloud;
  SplitString(gas, agas);
  SplitString(cloud, acloud);
  int ngas = agas.size();
  int ncloud = acloud.size();

  if (ngas == 0) {
    msg << "### FATAL ERROR in Molecule. Number of gas must be at least 1" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  myname = agas[0];
  Molecule *pmol;

  LoadChemistryFile(folder + myname + ".chem");
  phase = GAS;

  for (int i = 1; i < ngas; ++i) {
    std::string name = agas[i];
    pmol = AddMolecule(name);
    pmol->LoadChemistryFile(folder + name + ".chem");
    pmol->phase = GAS;
  }

  for (int i = 0; i < ncloud; ++i) {
    std::string name = acloud[i];
    pmol = AddMolecule(name);
    pmol->LoadChemistryFile(folder + name.substr(0, name.size() - 3) + ".chem");
    std::string p = name.substr(name.size() - 3, 3);
    if (p == "(l)")
      pmol->SetPhase(LIQUID);
    else if (p == "(s)")
      pmol->SetPhase(SOLID);
    else {
      msg << "### FATAL ERROR in Molecule " << name << ". Phase not found"
      << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
}

// destructor
Molecule::~Molecule()
{
  if (prev != NULL) prev->next = next;
  if (next != NULL) next->prev = prev;
}

Molecule* Molecule::AddMolecule(std::string name)
{
  std::stringstream msg;
  Molecule *p = this;

  if (p == NULL) {
    msg << "### FATAL ERROR in Molecule::AddMolecule. Molecule is empty. "
        << "Use new Molecle instred" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  while (p->next != NULL) p = p->next;
  p->next = new Molecule(name);
  p->next->prev = p;
  p->next->next = NULL;

  return p->next;
}

void Molecule::LoadChemistryFile(std::string chemfile)
{
  std::stringstream inp(DecommentFile(chemfile));
  Real junk;

  inp >> myname >> mu
      >> m_entropy >> m_enthalpy >> m_gibbs >> m_cp
      >> tr >> pr >> tc >> pc; 
  inp >> m_nshomate;
  for (int i = 0; i < m_nshomate; i++) {
    inp >> m_shomate_sp[i] >> junk;
    for (int j = 0; j < NSHOMATE; j++)
      inp >> m_shomate[i][j];
  }
  m_shomate_sp[m_nshomate] = junk;

  inp >> m_cliq >> m_enliq
      >> m_csld >> m_ensld;

  mu *= 1.E-3;  // g/mol -> kg/mol
}

int Molecule::TotalNumber()
{
  int ntotal = 1;
  Molecule *p = this;
  while (p->next != NULL) {
    p = p->next;
    ntotal++;
  }
  return ntotal;
}

Real Molecule::Cp(Real T) const 
{
  int i = _locate(m_shomate_sp, T, m_nshomate + 1);
  if (i == m_nshomate) i = m_nshomate - 1; if (i < 0) i = 0;

  Real result;
  T /= 1.E3;
  result = m_shomate[i][0] + m_shomate[i][1]*T + m_shomate[i][2]*T*T + 
    m_shomate[i][3]*T*T*T + m_shomate[i][4]/(T*T);
  return result;
}

Real Molecule::Cv(Real T) const
{
  return Cp(T) - Globals::Rgas;
}

Real Molecule::Enthalpy(Real T) const
{
  int i = _locate(m_shomate_sp, T, m_nshomate + 1);
  if (i == m_nshomate) i = m_nshomate - 1; if (i < 0) i = 0;

  T /= 1.E3;
  return m_shomate[i][0]*T + 0.5 * m_shomate[i][1]*T*T +
    1./3. * m_shomate[i][2]*T*T*T + 1./4. * m_shomate[i][3]*T*T*T*T + 
    - m_shomate[i][4]/T + m_shomate[i][5];
}

Real Molecule::Entropy(Real T) const
{
  int i = _locate(m_shomate_sp, T, m_nshomate + 1);
  if (i == m_nshomate) i = m_nshomate - 1; if (i < 0) i = 0;

  T /= 1.E3;
  return m_shomate[i][0] * log(T) + m_shomate[i][1]*T +
    0.5 * m_shomate[i][2]*T*T + 1./3. * m_shomate[i][3]*T*T*T + 
    - m_shomate[i][4]/(2.*T*T) + m_shomate[i][6];
}

Real Molecule::Latent(Real T) const
{
  return m_latent + Enthalpy(T) - Enthalpy(tr) - m_cp * (T - tr) / 1.E3;
}

Real Molecule::SatVaporPres(Real T) const
{
  return pr * exp((1. - tr / T) * m_beta - m_gamma * log(T / tr));
};

Real Molecule::SatVaporTemp(Real P, Real Tmin, Real Tmax, Real precision) const
{
  Real Tsat;
  SatVaporTempSolver solver;
  solver.pmol = this;
  solver.pres = P;
  int error = _root(Tmin, Tmax, precision, &Tsat, solver);
  std::stringstream msg;
  if (error) {
    msg << "SatVaporTemp failed" << std::endl
        << "Pressure = " << P << " Pa" << std::endl
        << "Tmin = " << Tmin << std::endl
        << "Tmax = " << Tmax << std::endl
        << "Saturation vapor pressure at Tmin = " << SatVaporPres(Tmin) << std::endl
        << "Saturation vapor pressure at Tmax = " << SatVaporPres(Tmax) << std::endl
        << *this << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return Tsat;
}

void Molecule::SetPhase(PhaseID pid)
{
  if (pid == LIQUID) {
    myname += "(l)";
    m_cp = m_cliq;
    m_latent = m_enliq;
    m_beta = (m_latent * 1.E3 - (Cp(tr) - m_cp) * tr) / (Globals::Rgas * tr);
    m_gamma = (m_cp - Cp(tr)) / Globals::Rgas;
    phase = LIQUID;
  } else if (pid == SOLID) {
    myname += "(s)";
    m_cp = m_csld;
    m_latent = m_ensld;
    m_beta = (m_latent * 1.E3 - (Cp(tr) - m_cp) * tr) / (Globals::Rgas * tr);
    m_gamma = (m_cp - Cp(tr)) / Globals::Rgas;
    phase = SOLID;
  }
}

Molecule* Molecule::GetMolecule(std::string name)
{
  std::stringstream msg;
  Molecule *p = this;

  while ((p != NULL) && (p->myname != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in Molecule:GetMolecule " << name << " not found" <<
    std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return p;
}

Molecule const* Molecule::GetMolecule(std::string name) const
{
  std::stringstream msg;
  Molecule const *p = this;

  while ((p != NULL) && (p->myname != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in Molecule:GetMolecule " << name << " not found" <<
    std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return p;
}


int Molecule::GetMoleculeId(std::string name) const
{
  std::stringstream msg;
  Molecule const *p = this;
  int id = 0;

  while ((p != NULL) && (p->myname != name)) {
    p = p->next;
    id++;
  }
  if (p == NULL) {
    msg << "### FATAL ERROR in Molecule:GetMolecule " << name << " not found" <<
    std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return id;
}
