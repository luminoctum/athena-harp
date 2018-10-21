//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file plm.cpp
//  \brief  piecewise linear reconstruction

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

// WENO 3 interpolation
inline Real interp_weno3(Real phim1, Real phi, Real phip1) {
    Real p0 = (-1.0/2.0) * phim1 + (3.0/2.0) * phi;
    Real p1 = (1.0/2.0) * phi + (1.0/2.0) * phip1;

    Real beta1 = (phip1 - phi) * (phip1 - phi);
    Real beta0 = (phi - phim1) * (phi - phim1);

    Real alpha0 = (1.0/3.0) /((beta0 + 1e-10) * (beta0 + 1.0e-10));
    Real alpha1 = (2.0/3.0)/((beta1 + 1e-10) * (beta1 + 1.0e-10));

    Real alpha_sum_inv = 1.0/(alpha0 + alpha1);

    Real w0 = alpha0 * alpha_sum_inv;
    Real w1 = alpha1 * alpha_sum_inv;

    return w0 * p0 + w1 * p1;
};

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX1()
//  \brief 

void Reconstruction::HighResFuncX1(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int n=0; n<NHYDRO; ++n) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      ql(n,i) = interp_weno3(q(n,k,j,i-2), q(n,k,j,i-1), q(n,k,j,i));
      qr(n,i) = interp_weno3(q(n,k,j,i+1), q(n,k,j,i), q(n,k,j,i-1));
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX2()
//  \brief 

void Reconstruction::HighResFuncX2(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int n=0; n<NHYDRO; ++n) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      ql(n,i) = interp_weno3(q(n,k,j-2,i), q(n,k,j-1,i), q(n,k,j,i));
      qr(n,i) = interp_weno3(q(n,k,j+1,i), q(n,k,j,i), q(n,k,j-1,i));
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::ReconstructionFuncX3()
//  \brief 

void Reconstruction::HighResFuncX3(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &ql, AthenaArray<Real> &qr)
{
  for (int n=0; n<NHYDRO; ++n) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      ql(n,i) = interp_weno3(q(n,k-2,j,i), q(n,k-1,j,i), q(n,k,j,i));
      qr(n,i) = interp_weno3(q(n,k+1,j,i), q(n,k,j,i), q(n,k-1,j,i));
    }
  }

  return;
}
