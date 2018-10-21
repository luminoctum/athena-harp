//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hs94.cpp
//  \brief Problem generator for Held-Suarez-94 GCM bench mark.
//
// REFERENCE: I.M Held & M.J Suarez, "A Proposal for the Intercomparison of the
// Dynamical Cores of Atmospheric General Circulation Models"

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../athena_math.hpp" // _root
#include "../utils/utils.hpp"
#include "../physics/held_suarez_94.hpp"  // HeldSuarez94

// made global to share with BC and source functions
static Real grav_acc;
static HeldSuarez94 hs;

// functions for boundary conditions

void ProjectPressureInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ProjectPressureOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

// functions for Held-Suarez forcing
void HeldSuarezForcing(MeshBlock *pmb, Real const time, Real const dt,
    AthenaArray<Real> const &prim, AthenaArray<Real> const &bcc, AthenaArray<Real> &cons)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real sigma = prim(IPR,k,j,i) / prim(IPR,k,j,pmb->is);
        Real theta = M_PI/2. - pmb->pcoord->x2v(j);
        //std::cout << dt << " " << hs.Kv(sigma) << " " << hs.Kt(theta, sigma) << std::endl;
        // archive momentum
        Real m1 = cons(IM1,k,j,i),
             m2 = cons(IM2,k,j,i),
             m3 = cons(IM3,k,j,i);

        Real temp = prim(IPR,k,j,i) / (hs.GetRgas() * prim(IDN,k,j,i));
        Real cv = hs.GetRgas() / (pmb->peos->GetGamma() - 1.);

        // Rayleigh friction
        cons(IM1,k,j,i) += - dt * hs.Kv(sigma) * m1; 
        cons(IM2,k,j,i) += - dt * hs.Kv(sigma) * m2;
        cons(IM3,k,j,i) += - dt * hs.Kv(sigma) * m3;
        
        // Newtonian cooling
        //cons(IEN,k,j,i) += - dt * hs.Kt(theta, sigma) * (temp - hs.GetTempEq(theta, prim(IPR,k,j,i)))
        //  * cons(IDN,k,j,i) * cv;
        
        // debug
        //std::cout << prim(IPR,k,j,i) << " " << prim(IDN,k,j,i) << std::endl;
        //std::cout << temp << " " << hs.GetTempEq(theta, prim(IPR,k,j,i)) << std::endl;
        //std::cout << cv << " " << hs.Kt(theta, sigma) << " " << dt << " " << cons(IDN,k,j,i) << std::endl;
      }
}


//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // 3D problem
  // Enroll special BCs
  //EnrollUserBoundaryFunction(INNER_X1, ProjectPressureInnerX1);
  //EnrollUserBoundaryFunction(OUTER_X1, ProjectPressureOuterX1);

  // Enroll source term
  //EnrollUserExplicitSourceFunction(HeldSuarezForcing);
}

//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Held-Suarez problem generator
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // pertubation field
  long int iseed = -1;
  Real k3 = 10.*M_PI/(pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min);

  // setup global variable shared by boundary condition
  grav_acc = - phydro->psrc->GetG1();

  // setup initial pressure/temperature field
  hs.LoadInputFile(pin);
  Real p1, t1;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      Real theta = M_PI/2. - pcoord->x2v(j);
      hs.SetLatitude(theta);
      for (int i = is; i <= ie; ++i) {
        phydro->w(IVX,k,j,i) = 0.;
        phydro->w(IVY,k,j,i) = 0.;
        phydro->w(IVZ,k,j,i) = (ran2(&iseed) - 0.5)*(1. + cos(k3*pcoord->x3v(k)));
        if (i == is) {
          p1 = hs.GetSurfacePressure();
        } else {
          hs.SetDistance(pcoord->dx1f(i));
          int err = _root(p1, 0.5 * p1, 1., &p1, hs);

          if (err != 0)
            throw std::runtime_error("FATAL ERROR: HS94 hydrostatic integration not converge");
        }
        Real t1 = hs.GetTempEq(theta, p1);
        phydro->w(IDN,k,j,i) = p1 / (hs.GetRgas() * t1);
        phydro->w(IPR,k,j,i) = p1;
        hs.SetBottomPressure(p1);
      }
    }

  // transfer to conservative variables
  // bcc is cell-centered magnetic fields, it is only a place holder here
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}

//! \fn void ProjectPressureInnerX1()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm
void ProjectPressureInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int n = 0; n < NHYDRO; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          if (n == IVX) {
            prim(IVX,k,j,is - i) = - prim(IVX,k,j,is + i - 1);  // reflect 1-vel
          } else if (n == IPR) {
            prim(IPR,k,j,is - i) = prim(IPR,k,j,is + i - 1)
              + prim(IDN,k,j,is + i - 1) * grav_acc * (2 * i - 1) * pco->dx1f(i);
          } else {
            prim(n,k,j,is - i) = prim(n,k,j,is + i - 1);
          }
        }

  return;
}

//! \fn void ProjectPressureOuterX1()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int n = 0; n < NHYDRO; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          if (n == IVX) {
            prim(IVX,k,j,ie + i) = - prim(IVX,k,j,ie - i + 1);  // reflect 1-vel
          } else if (n == IPR) {
            prim(IPR,k,j,ie + i) = prim(IPR,k,j,ie - i + 1)
              - prim(IDN,k,j,ie - i + 1) * grav_acc * (2 * i - 1) * pco->dx1f(i);
          } else {
            prim(n,k,j,ie + i) = prim(n,k,j,ie - i + 1);
          }
        }

  return;
}
