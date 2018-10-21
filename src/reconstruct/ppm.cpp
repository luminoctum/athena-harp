//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ppm.cpp
//  \brief piecewise parabolic reconstruction, adapted from remap in CMHOG code

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
//! \fn HydroIntegrator::ReconstructionFuncX1()
//  \brief 

/*
void HydroIntegrator::PiecewiseParabolicX1(const int n, const int m, const int k,
  const int j, const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{

// Compute average linear slopes (eqn 1.7)
#pragma simd
  for (int i=is; i<=(ie+1); ++i){
    Real dplus = q(n,i+1,j,k) - q(n,i  ,j,k);
    Real dmnus = q(n,i  ,j,k) - q(n,i-1,j,k);
    dd(i) = c1_(i)*dplus + c2_(i)*dmnus;

// Monotonize (eqn 1.8)
    Real qa = std::min( fabs(dd(i)), 2.0*std::min(fabs(dmnus),fabs(dplus)) ); 
    dd(i) = 0.0;
    if (dplus*dmnus > 0.0) dd(i) = qa*SIGN(dd(i));
  }
 
// construct interface values (eqn 1.6)
 
  for (int i=is-1; i<=(ie+2); ++i){
    dph(i)= c3_(i)*q(i-1,j,k) + c4_(i)*q(i,j,k) + c5_(i)*dd(i-1) + c6_(i)*dd(i);
  }
 
// left and right values

  for (int i=is-1; i<=ie+1; ++i){
    Real qa = (dph(i+1) - q(i,j,k))*(q(i,j,k) - dph(i))
    Real qd = dph(i+1)-dph(i)
    Real qe = 6.0*(q(i,j,k) - 0.5*(dph(i+1) + dph(i)))
    if (qa <= 0.0) {
      ql(i) = q(i,j,k);
      qr(i) = q(i,j,k);
    } else {
      ql(i) = dph(i  );
      qr(i) = dph(i+1);
    }
    if ( (qd*(qd - qe)) < 0.0) ql(i) = 3.0*q(i,j,k) - 2.0*dr(i);
    if ( (qd*(qd + qe)) < 0.0) qr(i) = 3.0*q(i,j,k) - 2.0*dl(i);
  }
      
  return;
}

//----------------------------------------------------------------------------------------
//! \fn HydroIntegrator::ReconstructionFuncX2()
//  \brief 

void HydroIntegrator::PiecewiseParabolicX2(const int n, const int m, const int k,
  const int j, const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  return;
}

//----------------------------------------------------------------------------------------
//! \fn HydroIntegrator::ReconstructionFuncX3()
//  \brief 

void HydroIntegrator::PiecewiseParabolicX3(const int n, const int m, const int k,
  const int j, const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  return;
}
*/
