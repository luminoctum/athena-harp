//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file minkowski.cpp
//  \brief implements functions for Minkowski (flat) spacetime and Cartesian (t,x,y,z)
//  coordinates in a derived class of the Coordinates abstract base class.
//
// Notes:
//   coordinates: t, x, y, z
//   metric: ds^2 = -dt^2 + dx^2 + dy^2 + dz^2

// C++ headers
#include <cmath>  // sqrt

// Athena++ headers
#include "coordinates.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../math_funcs.hpp"

//----------------------------------------------------------------------------------------
// Minkowski constructor
// Inputs:
//   pmb: pointer to block containing this grid
//   pin: pointer to runtime inputs (not used)
//   flag: true if object is for coarse grid only in an AMR calculation

Minkowski::Minkowski(MeshBlock *pmb, ParameterInput *pin, bool flag)
  : Coordinates(pmb, pin, flag)
{
  // Set indices
  pmy_block = pmb;
  coarse_flag = flag;
  int il, iu, jl, ju, kl, ku, ng;
  if (coarse_flag == true) {
    il = pmb->cis;
    iu = pmb->cie;
    jl = pmb->cjs;
    ju = pmb->cje;
    kl = pmb->cks;
    ku = pmb->cke;
    ng = pmb->cnghost;
  } else {
    il = pmb->is;
    iu = pmb->ie;
    jl = pmb->js;
    ju = pmb->je;
    kl = pmb->ks;
    ku = pmb->ke;
    ng = NGHOST;
  }
  Mesh *pm = pmy_block->pmy_mesh;
  RegionSize& mesh_size = pmy_block->pmy_mesh->mesh_size;
  RegionSize& block_size = pmy_block->block_size;

  // Allocate arrays for volume-centered coordinates and positions of cells
  int ncells1 = (iu-il+1) + 2*ng;
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = (ju-jl+1) + 2*ng;
  if (block_size.nx3 > 1) ncells3 = (ku-kl+1) + 2*ng;
  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);
  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);

  // Allocate arrays for area weighted positions for AMR/SMR MHD
  if (pm->multilevel && MAGNETIC_FIELDS_ENABLED) {
    x1s2.NewAthenaArray(ncells1);
    x1s3.NewAthenaArray(ncells1);
    x2s1.NewAthenaArray(ncells2);
    x2s3.NewAthenaArray(ncells2);
    x3s1.NewAthenaArray(ncells3);
    x3s2.NewAthenaArray(ncells3);
  }

  // Initialize volume-averaged coordinates and spacings: x-direction
  for (int i = il-ng; i <= iu+ng; ++i) {
    x1v(i) = 0.5 * (x1f(i) + x1f(i+1));
  }
  for (int i = il-ng; i <= iu+ng-1; ++i) {
    dx1v(i) = x1v(i+1) - x1v(i);
  }

  // Initialize volume-averaged coordinates and spacings: y-direction
  if (pmb->block_size.nx2 == 1) {
    x2v(jl) = 0.5 * (x2f(jl) + x2f(jl+1));
    dx2v(jl) = dx2f(jl);
  } else {
    for (int j = jl-ng; j <= ju+ng; ++j) {
      x2v(j) = 0.5 * (x2f(j) + x2f(j+1));
    }
    for (int j = jl-ng; j <= ju+ng-1; ++j) {
      dx2v(j) = x2v(j+1) - x2v(j);
    }
  }

  // Initialize volume-averaged coordinates and spacings: z-direction
  if (pmb->block_size.nx3 == 1) {
    x3v(kl) = 0.5 * (x3f(kl) + x3f(kl+1));
    dx3v(kl) = dx3f(kl);
  } else {
    for (int k = kl-ng; k <= ku+ng; ++k) {
      x3v(k) = 0.5 * (x3f(k) + x3f(k+1));
    }
    for (int k = kl-ng; k <= ku+ng-1; ++k) {
      dx3v(k) = x3v(k+1) - x3v(k);
    }
  }

  // Initialize area-averaged coordinates used with MHD AMR
  if (pmb->pmy_mesh->multilevel && MAGNETIC_FIELDS_ENABLED) {
    for (int i = il-ng; i <= iu+ng; ++i) {
      x1s2(i) = x1s3(i) = x1v(i);
    }
    if (pmb->block_size.nx2 == 1) {
      x2s1(jl) = x2s3(jl) = x2v(jl);
    } else {
      for (int j = jl-ng; j <= ju+ng; ++j) {
        x2s1(j) = x2s3(j) = x2v(j);
      }
    }
    if (pmb->block_size.nx3 == 1) {
      x3s1(kl) = x3s2(kl) = x3v(kl);
    } else {
      for (int k = kl-ng; k <= ku+ng; ++k) {
        x3s1(k) = x3s2(k) = x3v(k);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// Destructor

Minkowski::~Minkowski()
{
  dx1v.DeleteAthenaArray();
  dx2v.DeleteAthenaArray();
  dx3v.DeleteAthenaArray();
  x1v.DeleteAthenaArray();
  x2v.DeleteAthenaArray();
  x3v.DeleteAthenaArray();
  if (pmy_block->pmy_mesh->multilevel && MAGNETIC_FIELDS_ENABLED) {
    x1s2.DeleteAthenaArray();
    x1s3.DeleteAthenaArray();
    x2s1.DeleteAthenaArray();
    x2s3.DeleteAthenaArray();
    x3s1.DeleteAthenaArray();
    x3s2.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
// Function for computing cell-centered metric coefficients
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void Minkowski::CellMetric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i) {
    g(I00,i) = -1.0;
    g(I11,i) = 1.0;
    g(I22,i) = 1.0;
    g(I33,i) = 1.0;
    g_inv(I00,i) = -1.0;
    g_inv(I11,i) = 1.0;
    g_inv(I22,i) = 1.0;
    g_inv(I33,i) = 1.0;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Functions for computing face-centered metric coefficients
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void Minkowski::Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i) {
    g(I00,i) = -1.0;
    g(I11,i) = 1.0;
    g(I22,i) = 1.0;
    g(I33,i) = 1.0;
    g_inv(I00,i) = -1.0;
    g_inv(I11,i) = 1.0;
    g_inv(I22,i) = 1.0;
    g_inv(I33,i) = 1.0;
  }
  return;
}

void Minkowski::Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i) {
    g(I00,i) = -1.0;
    g(I11,i) = 1.0;
    g(I22,i) = 1.0;
    g(I33,i) = 1.0;
    g_inv(I00,i) = -1.0;
    g_inv(I11,i) = 1.0;
    g_inv(I22,i) = 1.0;
    g_inv(I33,i) = 1.0;
  }
  return;
}

void Minkowski::Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i) {
    g(I00,i) = -1.0;
    g(I11,i) = 1.0;
    g(I22,i) = 1.0;
    g(I33,i) = 1.0;
    g_inv(I00,i) = -1.0;
    g_inv(I11,i) = 1.0;
    g_inv(I22,i) = 1.0;
    g_inv(I33,i) = 1.0;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Functions for transforming face-centered primitives to locally flat frame
// Inputs
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   bb1: 3D array of normal components B^1 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   transformation is trivial

void Minkowski::PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb1, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx)
{
  if (MAGNETIC_FIELDS_ENABLED) {
    #pragma simd
    for (int i = il; i <= iu; ++i) {
      bbx(i) = bb1(k,j,i);
    }
  }
  return;
}

void Minkowski::PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb2, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx)
{
  if (MAGNETIC_FIELDS_ENABLED) {
    #pragma simd
    for (int i = il; i <= iu; ++i) {
      bbx(i) = bb2(k,j,i);
    }
  }
  return;
}

void Minkowski::PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb3, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx)
{
  if (MAGNETIC_FIELDS_ENABLED) {
    #pragma simd
    for (int i = il; i <= iu; ++i) {
      bbx(i) = bb3(k,j,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming fluxes to global frame: X-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates (unused)
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index

void Minkowski::FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux)
{
  #pragma simd
  for (int i = il; i <= iu; ++i) {
    const Real &txt = flux(IEN,i);
    Real &t10 = flux(IEN,i);
    t10 = -txt;
  }
  return;
}

void Minkowski::FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux)
{
  #pragma simd
  for (int i = il; i <= iu; ++i) {
    const Real &tyt = flux(IEN,i);
    Real &t20 = flux(IEN,i);
    t20 = -tyt;
  }
  return;
}

void Minkowski::FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux)
{
  #pragma simd
  for (int i = il; i <= iu; ++i) {
    const Real &tzt = flux(IEN,i);
    Real &t30 = flux(IEN,i);
    t30 = -tzt;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming cell-centered 4-vector from Minkowski to global
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j,i: z-, y-, and x-indices (unused)
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial

void Minkowski::TransformVectorCell(Real at, Real ax, Real ay, Real az, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = at;
  *pa1 = ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

//----------------------------------------------------------------------------------------
// Functions for transforming face-centered 4-vectors from Minkowski to global
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j,i: z-, y-, and x-indices (unused)
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial

void Minkowski::TransformVectorFace1(Real at, Real ax, Real ay, Real az, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = at;
  *pa1 = ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

void Minkowski::TransformVectorFace2(Real at, Real ax, Real ay, Real az, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = at;
  *pa1 = ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

void Minkowski::TransformVectorFace3(Real at, Real ax, Real ay, Real az, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = at;
  *pa1 = ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

//----------------------------------------------------------------------------------------
// Function for raising covariant components of a vector
// Inputs:
//   a_0,a_1,a_2,a_3: covariant components of vector
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to contravariant 4-vector components

void Minkowski::RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = -a_0;
  *pa1 = a_1;
  *pa2 = a_2;
  *pa3 = a_3;
  return;
}

//----------------------------------------------------------------------------------------
// Function for lowering contravariant components of a vector
// Inputs:
//   a0,a1,a2,a3: contravariant components of vector
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa_0,pa_1,pa_2,pa_3: pointers to covariant 4-vector components

void Minkowski::LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j,
    int i, Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3)
{
  *pa_0 = -a0;
  *pa_1 = a1;
  *pa_2 = a2;
  *pa_3 = a3;
  return;
}

//----------------------------------------------------------------------------------------
// Function for returning Minkowski coordinates of given cell in GR
// Inputs:
//   x0,x1,x2,x3: global coordinates to be converted
// Outputs:
//   pt,px,py,pz: variables pointed to set to Minkowski coordinates
// Notes:
//   conversion is trivial
//   useful to have if other coordinate systems for Minkowski space are developed

void Minkowski::GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
    Real *px, Real *py, Real *pz)
{
  *pt = x0;
  *px = x1;
  *py = x2;
  *pz = x3;
  return;
}

//----------------------------------------------------------------------------------------
// Function for returning spatial separation between points at same time
// Inputs:
//   x1,x2,x3: spatial coordinates of one point
//   y1,y2,y3: spatial coordinates of other point
// Outputs:
//   returned value: spatial separation between x and y
// Notes:
//   distance function is Euclidean

Real Minkowski::DistanceBetweenPoints(Real x1, Real x2, Real x3, Real y1, Real y2,
    Real y3)
{
  return std::sqrt(_sqr(x1-y1) + _sqr(x2-y2) + _sqr(x3-y3));
}
