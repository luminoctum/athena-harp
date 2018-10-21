// ! \file radconv.cpp
//   \brief radiative-convective model

// C++ headers
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../chemistry/molecule.hpp"
#include "../chemistry/condensation.hpp"
#include "../chemistry/chem_solver.hpp"
#include "../radiation/absorber.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/disort_wrapper.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../misc.hpp"
#include "../jupiter/jupiter_eos.hpp"

double c_planck_func2(double wnumlo, double wnumhi, double t) { return t;}

Real LogPressureCoordinate(Real x, RegionSize rs)
{
  return pow(rs.x1min,1.-x)*pow(rs.x1max,x);
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserMeshGenerator(X1DIR, LogPressureCoordinate);
}

void Radiation::InitAbsorbers(ParameterInput *pin, int)
{
  Molecule *pmol = pmy_block->pmol;
  pabs = new MwrAbsorberCIA(pmol);
  pabs->AddAbsorber(MwrAbsorberH2OKarpowicz(pmol));
  pabs->AddAbsorber(MwrAbsorberNH3Hanley(pmol));
  pabs->AddAbsorber(MwrAbsorberH2S(pmol));
  pabs->AddAbsorber(MwrAbsorberPH3Hoffman(pmol));
}

// \brief Problem generator for radiative model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::stringstream msg;

	// Step 0. define equation of state
	delete peos;
	peos = new JupiterEOS(this, pin);

  // Step 1. define molecules
  pmol = new Molecule(pin);
  if (pmol->TotalNumber() != NCOMP) {
    msg << "### FATAL ERROR in ProblemGenerator: pmol->TotalNumber() != NMOL" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  Real mixr[NCOMP], sum = 0;
  Molecule *pm = pmol;
  for (int n = 0; n < NGAS; ++n) {
    Real solar = pin->GetReal("chemistry", pm->myname + ".solar");
    Real enrich = pin->GetReal("chemistry", pm->myname + ".enrich");
    mixr[n] = solar*enrich;
    pm = pm->next;
    sum += mixr[n];
  }
  for (int n = NGAS; n < NCOMP; ++n) mixr[n] = 0.;
  for (int n = 0; n < NCOMP; ++n) mixr[n] /= sum;

  // Step 2. define reactions
  prg = new ReactionGroup(this, "condensation");
  prg->AddReaction(pin, "chemistry", "r1", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r2", pmol, &LiquidSolidIdeal);
  prg->AddReaction(pin, "chemistry", "r3", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r4", pmol, &GasGasSolidNH4SH);
  prg->AddReaction(pin, "chemistry", "r5", pmol, &GasCloudIdeal);

  // Step 3. set primitive variables (T/P, mixing ratios, ...)
  Real prim[NHYDRO];
  for (int n = 1; n < NCOMP; ++n)
    prim[n] = mixr[n];
  for (int n = NCOMP; n < NHYDRO; ++n)
    prim[n] = 0.;
  prim[IT] = pin->GetReal("problem", "tref");
  prim[IPR] = pin->GetReal("problem", "pref");
  prg->EquilibrateTP(prim);

  IsentropicSolver solver;
  solver.peos = peos;
  solver.prg = prg;
  solver.prim = prim;
  solver.max_iter = 20;
  solver.entropy = peos->Entropy(prim);

  int i0 = _locate(pcoord->x1v.data(), prim[IPR], pcoord->x1v.GetSize());

  for (int i = i0 + 1; i <= ie; ++i) {
    Real temp;
    prim[IPR] = pcoord->x1v(i);
    int err = _root(prim[IT], prim[IT]*2, 1.E-4, &temp, solver);
    if (err) {
      msg << "### FATAL ERROR: IsentropicSolver doesn't converge at level " << i << std::endl
          << "Try to change the temperature bracket." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    prim[IT] = temp;
    for (int n = 0; n < NHYDRO; ++n)
      phydro->w(n,ks,js,i) = prim[n];
  }

  prim[IT] = phydro->w(IT,ks,js,i0 + 1);
  for (int i = i0; i >= is; --i) {
    Real temp;
    prim[IPR] = pcoord->x1v(i);
    int err = _root(prim[IT]/2, prim[IT], 1.E-4, &temp, solver);
    if (err) {
      msg << "### FATAL ERROR: IsentropicSolver doesn't converge at level " << i << std::endl
          << "Try to change the temperature bracket." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    prim[IT] = temp;
    for (int n = 0; n < NHYDRO; ++n)
      phydro->w(n,ks,js,i) = prim[n];
  }

  // Step 4. convert primitive variables to conservative variables
  // bcc is cell-centered magnetic fields, it is only a place holder here
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  // Step 5. initialize radiation and abosrbers
  prad = new Radiation(this, pin, "mwr");
  prad->InitAbsorbers(pin);

  // Step 6. calculate vertical distance and set optical property (tau, ssalb, pmom, ...)
  Real grav = pin->GetReal("problem", "grav");
  prad->level[0] = 0.;
  for (int i = 1; i < prad->nlevel; ++i) {
    Real rho1 = 0., rho2 = 0.;
    for (int n = 0; n < NCOMP; ++n) {
      rho1 += phydro->u(n,ks,js,is + i - 1);
      rho2 += phydro->u(n,ks,js,is + i);
    }
    prad->level[i] = prad->level[i - 1] -
      2./grav*(phydro->w(IPR,ks,js,is + i) - phydro->w(IPR,ks,js,is + i - 1))/(rho1 + rho2);
  }
  prad->SetOpticalProperties(phydro->w, is, ie, js, js, ks, ks);
  prad->WriteOpticalDepth("tauc.out");

  // Step 7. set boundary condition and run radiative transfer model
  Radiation *pr = prad;
  DisortWrapper disort(pin, pr->nwave);
  while (pr != NULL) {
    disort.Run(pr->rad, pr->uu, pr->oppr, pr->wave);
    pr = pr->next;
  }
  prad->WriteTopRadiance("top_radiance.out");
}
