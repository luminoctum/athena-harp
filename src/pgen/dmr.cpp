//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file dmr.cpp
//  \brief Problem generator for double Mach reflection test.
//  Only works for genuinely 2D hydro problems in X1-X2 plane with adiabatic EOS.
//
// REFERENCE: P. Woodward & P. Colella, "The numerical simulation of two-dimensional
// fluid flow with strong shocks", JCP, 54, 115, sect. IVc.

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/mesh_refinement.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

#include <iostream>
#include <cmath>

// DMRInnerX1() - sets BCs on inner-x1 (left edge) of grid.  
// DMRInnerX2() - sets BCs on inner-x2 (bottom edge) of grid.  
// DMROuterX2() - sets BCs on outer-x2 (top edge) of grid.  
void DMRInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void DMRInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void DMROuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
int RefinementCondition(MeshBlock *pmb);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // Enroll user-defined boundary functions
  EnrollUserBoundaryFunction(INNER_X1, DMRInnerX1);
  EnrollUserBoundaryFunction(INNER_X2, DMRInnerX2);
  EnrollUserBoundaryFunction(OUTER_X2, DMROuterX2);
  // Enroll user-defined AMR criterion
  if(adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initialize DMR test
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::stringstream msg;

  if (block_size.nx3 > 1) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "nx3=" 
        << block_size.nx3 << " but this test only works for 2D" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Initialize shock using parameters defined in Woodward & Colella
  Real d0 = 8.0;
  Real e0 = 291.25;
  Real u0 =  8.25*sqrt(3.0)/2.0;
  Real v0 = -8.25*0.5;
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      Real shock_pos = 0.1666666666 + pcoord->x2v(j)/sqrt((double)3.0);
      // upstream conditions
      phydro->u(IDN,ks,j,i) = 1.4;
      phydro->u(IEN,ks,j,i) = 2.5;
      phydro->u(IM1,ks,j,i) = 0.0;
      phydro->u(IM2,ks,j,i) = 0.0;
      // downstream conditions
      if (pcoord->x1v(i) < shock_pos) {
        phydro->u(IDN,ks,j,i) = d0;
        phydro->u(IEN,ks,j,i) = e0 + 0.5*d0*(u0*u0+v0*v0);
        phydro->u(IM1,ks,j,i) = d0*u0;
        phydro->u(IM2,ks,j,i) = d0*v0;
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void DMRInnerX1()
//  \brief Sets boundary condition on left X boundary (iib) for dmr test
//  Quantities at this boundary are held fixed at the downstream state

void DMRInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Real d0 = 8.0;
  Real e0 = 291.25;
  Real u0 =  8.25*sqrt(3.0)/2.0;
  Real v0 = -8.25*0.5;
  Real gamma = pmb->peos->GetGamma();
  Real p0=e0*(gamma-1.0);

  for (int j=js; j<=je; ++j) {
    for (int i=1;  i<=(NGHOST); ++i) {
      prim(IDN,ks,j,is-i) = d0;
      prim(IVX,ks,j,is-i) = u0;
      prim(IVY,ks,j,is-i) = v0;
      prim(IPR,ks,j,is-i) = p0;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DMRInnerX2()
//  \brief  Sets boundary condition on lower Y boundary (ijb) for dmr test.
//  Quantaties at this boundary are held fixed at the downstream state for
//  x1 < 0.16666666, and are reflected for x1 > 0.16666666

void DMRInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Real d0 = 8.0;
  Real e0 = 291.25;
  Real u0 =  8.25*sqrt(3.0)/2.0;
  Real v0 = -8.25*0.5;
  Real gamma = pmb->peos->GetGamma();
  Real p0=e0*(gamma-1.0);

  for (int j=1;  j<=(NGHOST); ++j) {
    for (int i=is; i<=ie; ++i) {
      if (pco->x1v(i) < 0.1666666666) {
        // fixed at downstream state
        prim(IDN,ks,js-j,i) = d0;
        prim(IVX,ks,js-j,i) = u0;
        prim(IVY,ks,js-j,i) = v0;
        prim(IPR,ks,js-j,i) = p0;
      } else {
        // reflected
        prim(IDN,ks,js-j,i) = prim(IDN,ks,js+(j-1),i);
        prim(IVX,ks,js-j,i) = prim(IVX,ks,js+(j-1),i);
        prim(IVY,ks,js-j,i) = -prim(IVY,ks,js+(j-1),i);
        prim(IPR,ks,js-j,i) = prim(IPR,ks,js+(j-1),i);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DMROuterX2()
//  \brief Sets TIME-DEPENDENT boundary condition on upper Y boundary (ojb) for dmr test
//  Quantaties at this boundary are held fixed at the downstream state for
//  x1 < 0.16666666+v1_shock*time, and at the upstream state for
//  x1 > 0.16666666+v1_shock*time

void DMROuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Real d0 = 8.0;
  Real e0 = 291.25;
  Real u0 =  8.25*sqrt(3.0)/2.0;
  Real v0 = -8.25*0.5;
  Real shock_pos = 0.1666666666 + (1. + 20.*time)/sqrt(3.0);
  Real gamma = pmb->peos->GetGamma();
  Real p0=e0*(gamma-1.0);
  Real p1=2.5*(gamma-1.0);

  for (int j=1;  j<=(NGHOST); ++j) {
    for (int i=is; i<=ie; ++i) {
      if (pco->x1v(i) < shock_pos) {
        // fixed at downstream state
        prim(IDN,ks,je+j,i) = d0;
        prim(IVX,ks,je+j,i) = u0;
        prim(IVY,ks,je+j,i) = v0;
        prim(IPR,ks,je+j,i) = p0;
      } else {
        // fixed at upstream state
        prim(IDN,ks,je+j,i) = 1.4;
        prim(IVX,ks,je+j,i) = 0.0;
        prim(IVY,ks,je+j,i) = 0.0;
        prim(IPR,ks,je+j,i) = p1;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief refinement condition: density and pressure curvature

int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps=0.0;
  int k=pmb->ks;
  for(int j=pmb->js; j<=pmb->je; j++) {
    for(int i=pmb->is; i<=pmb->ie; i++) {
      Real epsr= (std::abs(w(IDN,k,j,i+1)-2.0*w(IDN,k,j,i)+w(IDN,k,j,i-1))
                 +std::abs(w(IDN,k,j+1,i)-2.0*w(IDN,k,j,i)+w(IDN,k,j-1,i)))/w(IDN,k,j,i);
      Real epsp= (std::abs(w(IPR,k,j,i+1)-2.0*w(IPR,k,j,i)+w(IPR,k,j,i-1))
                 +std::abs(w(IPR,k,j+1,i)-2.0*w(IPR,k,j,i)+w(IPR,k,j-1,i)))/w(IPR,k,j,i);
      Real eps = std::max(epsr, epsp);
      maxeps = std::max(maxeps, eps);
    }
  }
  // refine : curvature > 0.01
  if(maxeps > 0.01) return 1;
  // derefinement: curvature < 0.005
  if(maxeps < 0.005) return -1;
  // otherwise, stay
  return 0;
}
