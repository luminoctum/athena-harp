//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file spherical_latlon.cpp
//  \brief implements functions for spherical latlon (lambda-theta-r) coordinates in a 
//  derived class of the Coordinates abstract base class.

// C/C++ headers
#include <math.h>  // pow, trig functions

// Athena++ headers
#include "coordinates.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../math_funcs.hpp"   // _sqr, _cub
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"

//----------------------------------------------------------------------------------------
// Spherical latlon coordinates constructor

SphericalLatlon::SphericalLatlon(MeshBlock *pmb, ParameterInput *pin, bool flag)
  : Coordinates(pmb, pin, flag)
{
  pmy_block = pmb;
  coarse_flag=flag;
  int il, iu, jl, ju, kl, ku, ng;
  if(coarse_flag==true) {
    il = pmb->cis; jl = pmb->cjs; kl = pmb->cks;
    iu = pmb->cie; ju = pmb->cje; ku = pmb->cke;
    ng=pmb->cnghost;
  } else {
    il = pmb->is; jl = pmb->js; kl = pmb->ks;
    iu = pmb->ie; ju = pmb->je; ku = pmb->ke;
    ng=NGHOST;
  }
  Mesh *pm=pmy_block->pmy_mesh;
  RegionSize& mesh_size  = pmy_block->pmy_mesh->mesh_size;
  RegionSize& block_size = pmy_block->block_size;

  // allocate arrays for volume-centered coordinates and positions of cells
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

  // allocate arrays for area weighted positions for AMR/SMR MHD
  if((pm->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    x1s2.NewAthenaArray(ncells1);
    x1s3.NewAthenaArray(ncells1);
    x2s1.NewAthenaArray(ncells2);
    x2s3.NewAthenaArray(ncells2);
    x3s1.NewAthenaArray(ncells3);
    x3s2.NewAthenaArray(ncells3);
  }

  // initialize volume-averaged coordinates and spacing
  // x1-direction: x1v = (\int phi dV / \int dV) = dphi/2
  for (int i=il-ng; i<=iu+ng; ++i) {
    x1v(i) = 0.5*(x1f(i+1) + x1f(i));
  }
  for (int i=il-ng; i<=iu+ng-1; ++i) {
    dx1v(i) = x1v(i+1) - x1v(i);
  }

  // x2-direction: x2v = (\int theta dV / \int dV) =
  //   d(cos[theta] + theta sin[theta])/d(sin[theta])
  if (pmb->block_size.nx2 == 1) {
    x2v(jl) = 0.5*(x2f(jl+1) + x2f(jl));
    dx2v(jl) = dx2f(jl);
  } else {
    for (int j=jl-ng; j<=ju+ng; ++j) {
      x2v(j) = ((cos(x2f(j+1)) + x2f(j+1)*sin(x2f(j+1))) -
                (cos(x2f(j  )) + x2f(j  )*sin(x2f(j  ))))/
                (sin(x2f(j+1)) - sin(x2f(j)));
    }
    for (int j=jl-ng; j<=ju+ng-1; ++j) {
      dx2v(j) = x2v(j+1) - x2v(j);
    }
  }

  // x3-direction: x3v = (\int r dV / \int dV) = d(r^4/4)/d(r^3/3)
  if (pmb->block_size.nx3 == 1) {
    x3v(kl) = 0.75*(pow(x3f(kl+1),4) - pow(x3f(kl),4))/(pow(x3f(kl+1),3) - pow(x3f(kl),3));
    dx3v(kl) = dx3f(kl);
  } else {
    for (int k=kl-ng; k<=ku+ng; ++k) {
      x3v(k) = 0.75*(pow(x3f(k+1),4) - pow(x3f(k),4))/(pow(x3f(k+1),3) - pow(x3f(k),3));
    }
    for (int k=kl-ng; k<=ku+ng-1; ++k) {
      dx3v(k) = x3v(k+1) - x3v(k);
    }
  }

  // initialize area-averaged coordinates used with MHD AMR
  if((pmb->pmy_mesh->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    for (int i=il-ng; i<=iu+ng; ++i) {
      x1s2(i) = x1s3(i) = x1v(i);
    }
    if (pmb->block_size.nx2 == 1) {
      x2s1(jl) = x2s3(jl) = x2v(jl);
    } else {
      for (int j=jl-ng; j<=ju+ng; ++j) {
        x2s1(j) = x2s3(j) = x2v(j);
      }
    }
    if (pmb->block_size.nx3 == 1) {
      x3s1(kl) = x3s2(kl) = (2.0/3.0)*(pow(x3f(kl+1),3) - pow(x3f(kl),3))
                          / (_sqr(x3f(kl+1)) - _sqr(x3f(kl)));
    } else {
      for (int k=kl-ng; k<=ku+ng; ++k) {
        x3s1(k) = x3s2(k) = (2.0/3.0)*(pow(x3f(k+1),3) - pow(x3f(k),3))
                            / (_sqr(x3f(k+1)) - _sqr(x3f(k)));
      }
    }
  }
}

// destructor

SphericalLatlon::~SphericalLatlon()
{
  dx1v.DeleteAthenaArray();
  dx2v.DeleteAthenaArray();
  dx3v.DeleteAthenaArray();
  x1v.DeleteAthenaArray();
  x2v.DeleteAthenaArray();
  x3v.DeleteAthenaArray();
  if((pmy_block->pmy_mesh->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    x1s2.DeleteAthenaArray();
    x1s3.DeleteAthenaArray();
    x2s1.DeleteAthenaArray();
    x2s3.DeleteAthenaArray();
    x3s1.DeleteAthenaArray();
    x3s2.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
// Edge?Length functions: compute physical length at cell edge-? as vector
// Edge1(i,j,k) located at (i,j-1/2,k-1/2), i.e. (x1v(i), x2f(j), x3f(k))
//
void SphericalLatlon::Edge1Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // length1 = r cos(theta) d(phi)
    len(i) = x3f(k)*cos(x2f(j))*dx1f(i);
  }
  return;
}

// Edge2(i,j,k) located at (i-1/2,j,k-1/2), i.e. (x1f(i), x2v(j), x3f(k))

void SphericalLatlon::Edge2Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &len)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // length2 = r d(theta)
    len(i) = x3f(k)*dx2f(j);
  }
  return;
}

// Edge3(i,j,k) located at (i-1/2,j-1/2,k), i.e. (x1f(i), x2f(j), x3v(k))

//----------------------------------------------------------------------------------------
// GetEdge?Length functions: return length of edge-? at (i,j,k)
//
Real SphericalLatlon::GetEdge1Length(const int k, const int j, const int i)
{
  return x3f(k)*cos(x2f(j))*dx1f(i);
}

Real SphericalLatlon::GetEdge2Length(const int k, const int j, const int i)
{
  return x3f(k)*dx2f(j);
}

//----------------------------------------------------------------------------------------
// CenterWidthX functions: return physical width in X-dir at (i,j,k) cell-center

void SphericalLatlon::CenterWidth1(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx1)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    dx1(i) = x3v(k)*cos(x2v(j))*dx1f(i);
  }
  return;
}

void SphericalLatlon::CenterWidth2(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx2)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    dx2(i) = x3v(k)*dx2f(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector

void SphericalLatlon::Face1Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area1 = dr r dtheta = d(r^2/2) dtheta
    area(i) = 0.5*(_sqr(x3f(k+1)) - _sqr(x3f(k)))*dx2f(j);
  }
  return;
}

void SphericalLatlon::Face2Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area2 = dr r cos[theta] dphi = d(r^2/2) cos[theta] dphi
    area(i) = 0.5*(_sqr(x3f(k+1)) - _sqr(x3f(k)))*cos(x2f(j))*dx1f(i);
  }
  return;
}

void SphericalLatlon::Face3Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // area3 = r^2 cos[theta] dtheta dphi = r^2 d(sin[theta]) dphi
    area(i) = _sqr(x3f(k))*(sin(x2f(j+1)) - sin(x2f(j)))*dx1f(i); 
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)

Real SphericalLatlon::GetFace1Area(const int k, const int j, const int i)
{
  return 0.5*(_sqr(x3f(k+1)) - _sqr(x3f(k)))*dx2f(j);
}

Real SphericalLatlon::GetFace2Area(const int k, const int j, const int i)
{
  return 0.5*(_sqr(x3f(k+1)) - _sqr(x3f(k)))*cos(x2f(j))*dx1f(i);
}

Real SphericalLatlon::GetFace3Area(const int k, const int j, const int i)
{
  return _sqr(x3f(k))*(sin(x2f(j+1)) - sin(x2f(j)))*dx1f(i); 
}

//----------------------------------------------------------------------------------------
// Cell Volume function: compute volume of cell as vector

void SphericalLatlon::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &vol)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    // volume = r^2 cos(theta) dr dtheta dphi = d(r^3/3) d(sin theta) dphi
    vol(i) = 1./3.*(_cub(x3f(k+1)) - _cub(x3f(k)))*(sin(x2f(j+1)) - sin(x2f(j)))*dx1f(i);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)

Real SphericalLatlon::GetCellVolume(const int k, const int j, const int i)
{
  return 1./3.*(_cub(x3f(k+1)) - _cub(x3f(k)))*(sin(x2f(j+1)) - sin(x2f(j)))*dx1f(i);
}

//----------------------------------------------------------------------------------------
// Coordinate (Geometric) source term function

void SphericalLatlon::CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)
{
  Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();
  bool use_x2_fluxes = pmy_block->block_size.nx2 > 1;

  // Go through cells
  for (int k=pmy_block->ks; k<=pmy_block->ke; ++k) {
#pragma omp parallel for schedule(static)
    for (int j=pmy_block->js; j<=pmy_block->je; ++j) {
#pragma simd
      for (int i=pmy_block->is; i<=pmy_block->ie; ++i) {
        // src_3 = < M_{theta theta} + M_{phi phi} ><1/r>
        Real m_kk = prim(IDN,k,j,i)*(_sqr(prim(IM1,k,j,i)) + _sqr(prim(IM2,k,j,i)));

        if (EQUATION_OF_STATE == "adiabatic") {
          m_kk += 2.0*prim(IPR,k,j,i);
        } else if (EQUATION_OF_STATE == "shallow_water") {
          m_kk = 0.;
        } else if (EQUATION_OF_STATE == "isothermal") {
          m_kk += 2.0*(iso_cs*iso_cs)*prim(IDN,k,j,i);
        }

        if (MAGNETIC_FIELDS_ENABLED) {
           m_kk += _sqr(bcc(IB3,k,j,i));
        }

        Real darea = 0.5 * (_sqr(x3f(k+1)) - _sqr(x3f(k)));
        Real dvol  = 1./3. * (_cub(x3f(k+1)) - _cub(x3f(k)));
        u(IM3,k,j,i) += dt * darea/dvol * m_kk;

        if (pmy_block->block_size.nx3 > 1) {
          // src_2 = -< M_{theta r} ><1/r> 
          u(IM2,k,j,i) -= dt * dx3f(k)/((x3f(k) + x3f(k+1)) * dvol) *
            (_sqr(x3f(k)) * flux[X3DIR](IM2,k,j,i)
           + _sqr(x3f(k+1)) * flux[X3DIR](IM2,k+1,j,i));

          // src_1 = -< M_{phi r} ><1/r> 
          u(IM1,k,j,i) -= dt * dx3f(k)/((x3f(k) + x3f(k+1)) * dvol) *
            (_sqr(x3f(k)) * flux[X3DIR](IM1,k,j,i)
           + _sqr(x3f(k+1)) * flux[X3DIR](IM1,k+1,j,i));
        }

        // src_2 = < M_{phi phi} ><cot theta/r>
        Real m_pp = prim(IDN,k,j,i)*_sqr(prim(IM1,k,j,i));

        if (EQUATION_OF_STATE == "adiabatic") {
          m_pp += prim(IPR,k,j,i);
        } else if (EQUATION_OF_STATE == "shallow_water") {
          m_pp += 0.5 * _sqr(prim(IDN,k,j,i));
        } else if (EQUATION_OF_STATE == "isothermal") {
          m_pp += (iso_cs*iso_cs)*prim(IDN,k,j,i);
        }

        if (MAGNETIC_FIELDS_ENABLED) {
          m_pp += 0.5*( _sqr(bcc(IB3,k,j,i)) + _sqr(bcc(IB2,k,j,i)) - _sqr(bcc(IB1,k,j,i)) );
        }
        Real tan_theta = (cos(x2f(j)) - cos(x2f(j+1))) / (sin(x2f(j+1)) - sin(x2f(j)));
        u(IM2,k,j,i) -= dt * darea/dvol * tan_theta * m_pp;

        // src_1 = -< M_{phi theta} ><cot theta/r> 
        if (use_x2_fluxes) {
          u(IM1,k,j,i) += dt * darea/dvol * tan_theta / (cos(x2f(j)) + cos(x2f(j+1)))*
              (cos(x2f(j))*flux[X2DIR](IM1,k,j,i)
              + cos(x2f(j+1))*flux[X2DIR](IM1,k,j+1,i));
        }
        else {
          Real m_ph = prim(IDN,k,j,i) * prim(IM1,k,j,i) * prim(IM2,k,j,i);
          if (MAGNETIC_FIELDS_ENABLED) {
            m_ph -= bcc(IB3,k,j,i) * bcc(IB2,k,j,i);
          }
          u(IM1,k,j,i) += dt * darea/dvol * tan_theta * m_ph;
        }
      }
    }
  }

  return;
}
