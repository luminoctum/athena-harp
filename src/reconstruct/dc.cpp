//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file dc.cpp
//  \brief piecewise constant (donor cell) reconstruction

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX1()
//  \brief 

void Reconstruction::DonorCellX1(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  for (int n=0; n<NHYDRO; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(n,i) = w(n,k,j,i-1);
      wr(n,i) = w(n,k,j,i  );
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(IBY,i) = bcc(IB2,k,j,i-1);
      wl(IBZ,i) = bcc(IB3,k,j,i-1);
      wr(IBY,i) = bcc(IB2,k,j,i  );
      wr(IBZ,i) = bcc(IB3,k,j,i  );
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX2()
//  \brief 

void Reconstruction::DonorCellX2(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  for (int n=0; n<NHYDRO; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(n,i) = w(n,k,j-1,i);
      wr(n,i) = w(n,k,j  ,i);
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(IBY,i) = bcc(IB3,k,j-1,i);
      wl(IBZ,i) = bcc(IB1,k,j-1,i);
      wr(IBY,i) = bcc(IB3,k,j  ,i);
      wr(IBZ,i) = bcc(IB1,k,j  ,i);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX3()
//  \brief 

void Reconstruction::DonorCellX3(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  for (int n=0; n<NHYDRO; ++n){
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(n,i) = w(n,k-1,j,i);
      wr(n,i) = w(n,k  ,j,i);
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
#pragma simd
    for (int i=il; i<=iu; ++i){
      wl(IBY,i) = bcc(IB1,k-1,j,i);
      wl(IBZ,i) = bcc(IB2,k-1,j,i);
      wr(IBY,i) = bcc(IB1,k  ,j,i);
      wr(IBZ,i) = bcc(IB2,k  ,j,i);
    }
  }

  return;
}
