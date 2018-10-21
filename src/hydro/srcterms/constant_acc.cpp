//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief source terms due to constant acceleration (e.g. for RT instability)

// Athena++ headers
#include "hydro_srcterms.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../eos/eos.hpp"
#include "../../globals.hpp"

inline Real TotalDensity(AthenaArray<Real> const& prim, int i, int j, int k, Real const mu[])
{
  Real xt = 1., rho = 0.;
  for (int n = NGAS; n < NCOMP; ++n)
    xt -= prim(n,k,j,i);
  Real x1 = xt;
  for (int n = 1; n < NGAS; ++n) {
    rho += prim(n,k,j,i)*mu[n];
    x1 -= prim(n,k,j,i);
  }
  Real rho0 = x1/xt*prim(IPR,k,j,i)*mu[0]/(Globals::Rgas*prim(IT,k,j,i));
  rho *= prim(IPR,k,j,i)/(xt*Globals::Rgas*prim(IT,k,j,i));
  rho += rho0;
  for (int n = NGAS; n < NCOMP; ++n)
    rho += prim(n,k,j,i)*mu[n]/(x1*mu[0])*rho0;
  return rho;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::ConstantAcceleration
//  \brief Adds source terms for constant acceleration to conserved variables

void HydroSourceTerms::ConstantAcceleration(const Real dt,const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  for (int k=pmb->ks; k<=pmb->ke; ++k)
    #pragma omp parallel for schedule(static)
    for (int j=pmb->js; j<=pmb->je; ++j)
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real s1, s2, s3;
        if (EQUATION_OF_STATE == "adiabatic") {
          s1 = dt*prim(IDN,k,j,i)*g1_;
          s2 = dt*prim(IDN,k,j,i)*g2_;
          s3 = dt*prim(IDN,k,j,i)*g3_;
        } else if (EQUATION_OF_STATE == "heterogeneous") {
          Real rho = TotalDensity(prim, i, j, k, pmb->peos->mu_);
          s1 = dt*rho*g1_;
          s2 = dt*rho*g2_;
          s3 = dt*rho*g3_;
        }
        cons(IM1,k,j,i) += s1;
        cons(IM2,k,j,i) += s2;
        cons(IM3,k,j,i) += s3;
        if (NON_BAROTROPIC_EOS)
          cons(IEN,k,j,i) += s1*prim(IVX,k,j,i) + s2*prim(IVY,k,j,i) + s3*prim(IVZ,k,j,i);
      }

  return;
}
