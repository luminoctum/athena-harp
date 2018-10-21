//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_hydro.cpp
//  \brief implements functions in class HeterogeneousHydro for adiabatic hydrodynamics`

// C/C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN
#include <cstdlib>  // atof

// Athena++ headers
#include "eos.hpp"
#include "../hydro/hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../field/field.hpp"
#include "../math_funcs.hpp"
#include "../misc.hpp"
#include "../chemistry/reaction.hpp"

// HeterogeneousHydro constructor

HeterogeneousHydro::HeterogeneousHydro(MeshBlock *pmb, ParameterInput *pin):
  EquationOfState(pmb)
{
  std::stringstream msg;

  std::vector<std::string> str;
  // read mu, size should be NCOMP, unit is [g/mol], convert to [kg/mol]
  SplitString(pin->GetString("hydro", "mu"), str);
  if (str.size() != NCOMP) {
    msg << "### FATAL ERROR in heterogeneous_hydro.cpp::HeterogeneousHydro: number of gases in 'mu' does not equal NGAS" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  for (int i = 0; i < NCOMP; ++i)
    mu_[i] = atof(str[i].c_str())*1.E-3;

  // read cv, size should be NCOMP, unit is [J/(mol K)]
  SplitString(pin->GetString("hydro", "cv"), str);
  if (str.size() != NCOMP) {
    msg << "### FATAL ERROR in heterogeneous_hydro.cpp::HeterogeneousHydro: number of species in 'cv' does not equal NCOMP" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  for (int i = 0; i < NCOMP; ++i)
    cv_[i] = atof(str[i].c_str());

  // read latent, size should be NCOMP - NGAS, unit is [kJ/mol]
  if (NCOMP - NGAS > 0) {
    SplitString(pin->GetString("hydro", "latent"), str);
    if (str.size() != NCOMP - NGAS) {
      msg << "### FATAL ERROR in heterogeneous_hydro.cpp::HeterogeneousHydro: number of clouds in 'latent' does not equal " << NCOMP - NGAS << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
  for (int i = 0; i < NGAS; ++i)
    latent_[i] = 0.;
  for (int i = NGAS; i < NCOMP; ++i)
    latent_[i] = atof(str[i-NGAS].c_str())*1.E3;

  //density_floor_  = pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  density_floor_  = 0.;
  pressure_floor_ = pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));

  // default gamma
  gamma_ = (cv_[0] + Globals::Rgas)/cv_[0];
}

// destructor

HeterogeneousHydro::~HeterogeneousHydro()
{
}

//----------------------------------------------------------------------------------------
// \!fn void HeterogeneousHydro::ConservedToPrimitive(AthenaArray<Real> &cons,
//           const AthenaArray<Real> &prim_old, const FaceField &b,
//           AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke)
// \brief Converts conserved into primitive variables in adiabatic hydro.

void HeterogeneousHydro::ConservedToPrimitive(AthenaArray<Real> &cons,
  const AthenaArray<Real> &prim_old, const FaceField &b, AthenaArray<Real> &prim,
  AthenaArray<Real> &bcc, Coordinates *pco, int is, int ie, int js, int je, int ks, int ke)
{
  int nthreads = pmy_block_->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  ReactionGroup *prg = NULL;
  if (pmy_block_->prg != NULL)
    prg = pmy_block_->prg->GetReactionGroup("condensation");

  Real prim1[NHYDRO];
  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
  for (int j=js; j<=je; ++j){
    for (int i=is; i<=ie; ++i){
      Real& m1 = cons(IM1,k,j,i);
      Real& m2 = cons(IM2,k,j,i);
      Real& m3 = cons(IM3,k,j,i);
      Real rho = 0., rck = 0., rc = 0., ohr;

      // apply density floor, without changing momentum or energy
      for (int n = 0; n < NCOMP; ++n) {
        cons(n,k,j,i) = _max(cons(n,k,j,i), density_floor_);
        rck += cons(n,k,j,i)*Globals::Rgas/mu_[n];
        rc  += cons(n,k,j,i)*cv_[n]/mu_[n];
        rho += cons(n,k,j,i);
      }
      ohr = 1./rho;

      // molar mixing ratio
      for (int n = 1; n < NCOMP; ++n)
        prim1[n] = cons(n,k,j,i)*Globals::Rgas/(mu_[n]*rck);

      // subtract cloud components when calculate pressure
      Real latent = 0.;
      for (int n = NGAS; n < NCOMP; ++n) {
        rck -= cons(n,k,j,i)*Globals::Rgas/mu_[n];
        latent += cons(n,k,j,i)*latent_[n]/mu_[n];
      }

      // apply pressure floor, correct total energy
      prim1[IPR] = rck/rc*(cons(IEN,k,j,i)-latent-0.5*ohr*(_sqr(m1)+_sqr(m2)+_sqr(m3)));
      cons(IEN,k,j,i) = (prim1[IPR] > pressure_floor_) ?  cons(IEN,k,j,i) :
        (pressure_floor_/rck*rc + 0.5*ohr*(_sqr(m1)+_sqr(m2)+_sqr(m3) + latent));
      prim1[IPR] = _max(prim1[IPR], pressure_floor_);

      // temperature
      prim1[IT] = prim1[IPR]/rck;

      // velocity
      prim1[IVX] = ohr*m1;
      prim1[IVY] = ohr*m2;
      prim1[IVZ] = ohr*m3;

      // chemical equilibrium
      if (prg != NULL) {
        // for debug
        //Real prim2[NHYDRO];
        //for (int n = 0; n < NHYDRO; ++n) prim2[n] = prim1[n];
        int status = prg->EquilibrateUV(prim1, this);
        if (status == 0) { // not converged
          std::stringstream msg;
          msg << "### FATAL ERROR in ConservedToPrimitive: EquilibrateUV does not converge." 
              << std::endl;
          // for debug
          //msg << "prim before: ";
          //for (int n = 0; n < NHYDRO; ++n)
          //  msg << prim2[n] << " ";
          //msg << std::endl;
          msg << "prim after: ";
          for (int n = 0; n < NHYDRO; ++n)
            msg << prim1[n] << " ";
          throw std::runtime_error(msg.str().c_str());
        } else if (status == 2) { // converged and reacted
          // total gas mixing ratios
          Real xt = 1.;
          for (int n = NGAS; n < NCOMP; ++n) xt -= prim1[n];

          // gas density
          Real x1 = xt;
          for (int n = 1; n < NGAS; ++n) {
            cons(n,k,j,i) = mu_[n]*prim1[n]/xt*prim1[IPR]/(Globals::Rgas*prim1[IT]);
            x1 -= prim1[n];
          }
          cons(0,k,j,i) = mu_[0]*x1/xt*prim1[IPR]/(Globals::Rgas*prim1[IT]);

          // cloud density
          for (int n = NGAS; n < NCOMP; ++n)
            cons(n,k,j,i) = prim1[n]*mu_[n]/(x1*mu_[0])*cons(0,k,j,i);
        }
      }

      // update primitive variables
      for (int n = 0; n < NHYDRO; ++n) prim(n,k,j,i) = prim1[n];
    }
  }}
}

  return;
}


//----------------------------------------------------------------------------------------
// \!fn void HeterogeneousHydro::PrimitiveToConserved(const AthenaArray<Real> &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke);
// \brief Converts primitive variables into conservative variables

void HeterogeneousHydro::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
     int is, int ie, int js, int je, int ks, int ke)
{
  int nthreads = pmy_block_->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
  for (int j=js; j<=je; ++j){
    for (int i=is; i<=ie; ++i){
      Real const& v1 = prim(IVX,k,j,i);
      Real const& v2 = prim(IVY,k,j,i);
      Real const& v3 = prim(IVZ,k,j,i);

      // total gas mixing ratios
      Real xt = 1.;
      for (int n = NGAS; n < NCOMP; ++n)
        xt -= prim(n,k,j,i);

      // gas density
      Real x1 = xt, rho = 0., rc = 0.;
      for (int n = 1; n < NGAS; ++n) {
        cons(n,k,j,i) = mu_[n]*prim(n,k,j,i)/xt*prim(IPR,k,j,i)/(Globals::Rgas*prim(IT,k,j,i));
        x1 -= prim(n,k,j,i);
        rho += cons(n,k,j,i);
        rc += cons(n,k,j,i)*cv_[n]/mu_[n];
      }
      cons(0,k,j,i) = mu_[0]*x1/xt*prim(IPR,k,j,i)/(Globals::Rgas*prim(IT,k,j,i));
      rho += cons(0,k,j,i);
      rc += cons(0,k,j,i)*cv_[0]/mu_[0];

      // cloud density
      for (int n = NGAS; n < NCOMP; ++n) {
        cons(n,k,j,i) = prim(n,k,j,i)*mu_[n]/(x1*mu_[0])*cons(0,k,j,i);
        rho += cons(n,k,j,i);
        rc += cons(n,k,j,i)*cv_[n]/mu_[n];
      }

      // momentum
      cons(IM1,k,j,i) = rho*v1;
      cons(IM2,k,j,i) = rho*v2;
      cons(IM3,k,j,i) = rho*v3;

      // total energy
      cons(IEN,k,j,i) = prim(IT,k,j,i)*rc + 0.5*rho*(_sqr(v1)+_sqr(v2)+_sqr(v3));
      for (int n = NGAS; n < NCOMP; ++n)
        cons(IEN,k,j,i) += cons(n,k,j,i)*latent_[n]/mu_[n];
    }
  }}
}
  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real HeterogeneousHydro::SoundSpeed(Real const prim[])
// \brief returns adiabatic sound speed given vector of primitive variables

Real HeterogeneousHydro::SoundSpeed(Real const prim[]) const
{
  Real xt = 1., cv = 0., rho = 0.;
  for (int n = NGAS; n < NCOMP; ++n)
    xt -= prim[n];
  Real x1 = xt;
  for (int n = 1; n < NGAS; ++n) {
    rho += prim[n]/xt*prim[IPR]*mu_[n]/(Globals::Rgas*prim[IT]);
    cv += cv_[n]*prim[n];
    x1 -= prim[n];
  }
  rho += x1/xt*prim[IPR]*mu_[0]/(Globals::Rgas*prim[IT]);
  cv += cv_[0]*x1;
  for (int n = NGAS; n < NCOMP; ++n) {
    rho += prim[n]/x1*mu_[n]/mu_[0];
    cv += cv_[n]*prim[n];
  }
  Real kappa = xt*Globals::Rgas/cv;
  return sqrt((kappa + 1.)*prim[IPR]/rho);
}

Real HeterogeneousHydro::Cp(Real const prim[]) const
{
  Real x1 = 1., cp = 0;
  for (int n = 1; n < NGAS; ++n) {
    cp += (cv_[n] + Globals::Rgas)*prim[n];
    x1 -= prim[n];
  }
  for (int n = NGAS; n < NCOMP; ++n) {
    cp += cv_[n]*prim[n];
    x1 -= prim[n];
  }
  cp += (cv_[0] + Globals::Rgas)*x1;
  return cp;
}

Real HeterogeneousHydro::Cv(Real const prim[]) const
{
  Real x1 = 1., cv = 0;
  for (int n = 1; n < NCOMP; ++n) {
    cv += cv_[n]*prim[n];
    x1 -= prim[n];
  }
  cv += cv_[0]*x1;
  return cv;
}

Real HeterogeneousHydro::Mass(Real const prim[]) const
{
  Real x1 = 1., mu = 0;
  for (int n = 1; n < NCOMP; ++n) {
    mu += mu_[n]*prim[n];
    x1 -= prim[n];
  }
  mu += mu_[0]*x1;
  return mu;
}

Real HeterogeneousHydro::Entropy(Real const prim[]) const
{
  Real entropy = 0., x1 = 1., xt = 1.;
  // total gas mixing ratios
  for (int n = NGAS; n < NCOMP; ++n)
    xt -= prim[n];
  // gas entropy
  for (int n = 1; n < NGAS; ++n) {
    entropy += (cv_[n] + Globals::Rgas)*log(prim[IT])*prim[n]
      - Globals::Rgas*log(prim[n]/xt*prim[IPR])*prim[n];
    x1 -= prim[n];
  }
  // cloud entropy
  for (int n = NGAS; n < NCOMP; ++n) {
    entropy += ((cv_[n] + Globals::Rgas)*log(prim[IT])+ latent_[n]/prim[IT])*prim[n];
    x1 -= prim[n];
  }
  entropy += (cv_[0] + Globals::Rgas)*log(prim[IT])*x1
    - Globals::Rgas*log(x1/xt*prim[IPR])*x1;
  return entropy;
}

Real HeterogeneousHydro::Enthalpy(Real const prim[]) const
{
  Real enthalpy = 0., x1 = 1.;
  for (int n = 1; n < NGAS; ++n) {
    enthalpy += (cv_[n] + Globals::Rgas)*prim[IT]*prim[n];
    x1 -= prim[n];
  }
  for (int n = NGAS; n < NCOMP; ++n) {
    enthalpy += ((cv_[n] + Globals::Rgas)*prim[IT] + latent_[n])*prim[n];
    x1 -= prim[n];
  }
  enthalpy += (cv_[0] + Globals::Rgas)*x1*prim[IT];
  return enthalpy;
}

/*Real HeterogeneousHydro::Energy(Real const prim[]) const
{
  Real energy = 0., x1 = 1.;
  for (int n = 1; n < NGAS; ++n) {
    energy += cv_[n]*prim[IT]*prim[n];
    x1 -= prim[n];
  }
  for (int n = NGAS; n < NCOMP; ++n) {
    energy += (cv_[n]*prim[IT] + latent_[n])*prim[n];
    x1 -= prim[n];
  }
  energy += cv_[0]*x1*prim[IT];
  return energy; 
}

void HeterogeneousHydro::EquilibrateUV(
  Real prim[], Real u0, Real xt0, Real *cons) const
{
  Real xt = 1., cv = 0.;
  for (int n = NGAS; n < NCOMP; ++n) {
    u0 -= latent_[n]*prim[n];
    cv += cv_[n]*prim[n];
    xt -= prim[n];
  }
  Real x1 = xt;
  for (int n = 1; n < NGAS; ++n) {
    cv += cv_[n]*prim[n];
    x1 -= prim[n];
  }
  cv += cv_[0]*x1;
  prim[IPR] *= xt*(u0/cv)/(xt0*prim[IT]);
  prim[IT] = u0/cv;

  // conservative variable
  if (cons != NULL) {
    cons[0] = mu_[0]*x1/xt*prim[IPR]/(Globals::Rgas*prim[IT]);
    for (int n = 1; n < NGAS; ++n) // gas
      cons[n] = mu_[n]*prim[n]/xt*prim[IPR]/(Globals::Rgas*prim[IT]);
    for (int n = NGAS; n < NCOMP; ++n) // cloud
      cons[n] = prim[n]*mu_[n]/(x1*mu_[0])*cons[0];
  }
}*/
