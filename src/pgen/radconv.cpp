// ! \file radconv.cpp
//   \brief radiative-convective model

// C++ headers
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../chemistry/molecule.hpp"
#include "../chemistry/reaction.hpp"
#include "../chemistry/condensation.hpp"
#include "../radiation/absorber.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/disort_wrapper.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../misc.hpp"

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
  std::string kfile = pin->GetString("radiation", myname);
  std::string fcia1 = pin->GetString("radiation", "fcia1");
  std::string fcia2 = pin->GetString("radiation", "fcia2");
  
  pabs = new HitranAbsorber("CH4", pmol);
  pabs->LoadCoefficient(kfile);
  pabs->AddAbsorber(HitranAbsorber("C2H2", pmol))
      ->LoadCoefficient(kfile);

  pabs->AddAbsorber(CIAXiz("H2 H2", pmol))
      ->LoadCoefficient(fcia1);
  pabs->AddAbsorber(CIAXiz("H2 He", pmol))
      ->LoadCoefficient(fcia2);
}

// \brief Problem generator for radiative model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::stringstream msg;

  // Step 1. define molecules
  pmol = new Molecule(pin);
  if (pmol->TotalNumber() != NCOMP) {
    msg << "### FATAL ERROR in ProblemGenerator: pmol->TotalNumber() != NMOL" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Step 2. define reactions
  prg = new ReactionGroup(this, "condensation");
  prg->AddReaction(pin, "chemistry", "r1", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r2", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r3", pmol, &LiquidSolidIdeal);
  prg->AddReaction(pin, "chemistry", "r4", pmol, &GasGasSolidNH4SH);
  prg->AddReaction(pin, "chemistry", "r5", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r6", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r7", pmol, &LiquidSolidIdeal);

  // Step 3. set primitive variables (T/P, mixing ratios, ...)
  std::string afile = pin->GetString("radiation", "atm_file");
  NamedArray atm = ReadNamedArray(afile);
  std::vector<Real> &paxis = atm["PRE"];
  for (int i = is; i <= ie; ++i) {
    Real pres = pcoord->x1v(i)/100.;
    phydro->w(IT,i) = _interp1(pres, atm["TEM"].data(), paxis.data(), paxis.size(),1);
    Molecule *p = pmol->next; // exclude the mixing ratio of the first one
    for (int n = 1; n < NCOMP; ++n) {
      if (atm.find(p->myname) != atm.end()) // if find this molecule in the table
        phydro->w(n,i) = _interp1(pres, atm[p->myname].data(), paxis.data(), paxis.size(),1)/1.E6;
      else
        phydro->w(n,i) = 0.;
      p = p->next;
    }
    phydro->w(IVX,i) = 0.;
    phydro->w(IVY,i) = 0.;
    phydro->w(IVZ,i) = 0.;
    phydro->w(IPR,i) = pcoord->x1v(i);
  }

  // Step 4. convert primitive variables to conservative variables
  // bcc is cell-centered magnetic fields, it is only a place holder here
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  // Step 5. initialize radiation and abosrbers
  prad = new Radiation(this, pin, "b1");
  prad->InitAbsorbers(pin); 
  
  // Step 6. set solar flux
  std::string sfile = pin->GetString("radiation", "solar_file");
  Real dist = pin->GetReal("radiation", "distance");
  AthenaArray<Real> sdata(sfile);
  std::vector<Real> saxis(sdata.GetDim2()), solar(sdata.GetDim2());
  for (int i = 0; i < sdata.GetDim2(); ++i) {
    saxis[i] = sdata(i, 0);
    solar[i] = sdata(i, 1)/(dist*dist);
  }

  // Step 7. calculate vertical distance and set optical property (tau, ssalb, pmom, ...) 
  //         for each radiation band
  Real grav = pin->GetReal("problem", "grav");
#pragma omp parallel for
  for (int b = 0; b < 1; ++b) {
    Radiation *pb = prad->Next(b);

    pb->level[0] = 0.;
    for (int i = 1; i < pb->nlevel; ++i) {
      Real rho1 = 0., rho2 = 0.;
      for (int n = 0; n < NCOMP; ++n) {
        rho1 += phydro->u(n,is + i - 1);
        rho2 += phydro->u(n,is + i);
      }
      pb->level[i] = pb->level[i - 1] -
        2./grav*(phydro->w(IPR,is + i) - phydro->w(IPR,is + i - 1))/(rho1 + rho2);
    }
    pb->SetOpticalProperties(phydro->w, ks, js, is, ie);

    std::cout << "running disort for band " << b+1 << " ..." << std::endl;
    DisortWrapper disort(pin, pb->nwave, TOP2BOT);
    for (int i = 0; i < pb->nwave; ++i)
      disort.bc[i].fbeam = _interp1(pb->wave[i], solar.data(), saxis.data(), saxis.size());
    disort.Run(pb->rad, pb->uu, pb->oppr, pb->wave);
  }
  
  // Step 8. calculate flux and heating rate
  AthenaArray<Real> flux, hrate;
  prad->TotalFlux(flux);
  hrate.NewAthenaArray(flux.GetDim2(), flux.GetDim1());
  for (int i = 0; i < flux.GetDim2() - 1; ++i)
    for (int b = 0; b < flux.GetDim1(); ++b)
      hrate(i,b) = grav*(flux(i+1,b) - flux(i,b))/(phydro->w(IPR,is+i+1) - phydro->w(IPR,is+i));
  for (int b = 0; b < flux.GetDim1(); ++b)
    hrate(flux.GetDim2() - 1, b) = 0.;
  WriteHeatingRate("heating_rate.out", flux, hrate, prad->level);
}
