//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_blast.cpp
//  \brief Problem generator for GRMHD spherical blast wave in flat spacetime.

// C++ headers
#include <algorithm>  // min()

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // macros, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro

// Configuration checking
#if not GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  // Get ratio of specific heats
  Real gamma_adi = peos->GetGamma();
  Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read problem parameters
  Real num_x = pin->GetReal("problem", "num_x");
  Real num_y = pin->GetReal("problem", "num_y");
  Real x_spacing = pin->GetReal("problem", "x_spacing");
  Real y_spacing = pin->GetReal("problem", "y_spacing");
  Real radius = pin->GetReal("problem", "radius");
  Real rho_inner = pin->GetReal("problem", "rho_inner");
  Real pgas_inner = pin->GetReal("problem", "pgas_inner");
  Real rho_outer = pin->GetReal("problem", "rho_outer");
  Real pgas_outer = pin->GetReal("problem", "pgas_outer");
  Real bx = 0.0, by = 0.0, bz = 0.0;
  if (MAGNETIC_FIELDS_ENABLED) {
    bx = pin->GetReal("problem", "bx");
    by = pin->GetReal("problem", "by");
    bz = pin->GetReal("problem", "bz");
  }

  // Prepare auxiliary array
  int ncells1 = block_size.nx1 + 2*NGHOST;
  int ncells2 = block_size.nx2;
  if (ncells2 > 1) {
    ncells2 += 2*NGHOST;
  }
  int ncells3 = block_size.nx3;
  if (ncells3 > 1) {
    ncells3 += 2*NGHOST;
  }
  AthenaArray<Real> b, g, gi;
  b.NewAthenaArray(3, ncells3, ncells2, ncells1);
  g.NewAthenaArray(NMETRIC, ncells1);
  gi.NewAthenaArray(NMETRIC, ncells1);

  // Initialize hydro variables
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i = il; i <= iu; ++i) {

        // Calculate distance to nearest blast center
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real min_separation = pcoord->DistanceBetweenPoints(x1, x2, x3, 0.0, 0.0, 0.0);
        for (int x_index = -num_x; x_index <= num_x; ++x_index) {
          Real center_x = x_index * x_spacing;
          for (int y_index = -num_y; y_index <= num_y; ++y_index) {
            if (x_index == 0 && y_index == 0) {
              continue;
            }
            Real center_y = y_index * y_spacing;
            Real separation = pcoord->DistanceBetweenPoints(x1, x2, x3, center_x,
                center_y, 0.0);
            min_separation = std::min(min_separation, separation);
          }
        }

        // Set pressure and density
        Real rho, pgas;
        if (min_separation < radius) {
          rho = rho_inner;
          pgas = pgas_inner;
        } else {
          rho = rho_outer;
          pgas = pgas_outer;
        }

        // Set velocity
        Real ut = 1.0;
        Real ux = 0.0;
        Real uy = 0.0;
        Real uz = 0.0;
        Real u0, u1, u2, u3;
        pcoord->TransformVectorCell(ut, ux, uy, uz, k, j, i, &u0, &u1, &u2, &u3);
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = u1 - gi(I01,i)/gi(I00,i) * u0;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = u2 - gi(I02,i)/gi(I00,i) * u0;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = u3 - gi(I03,i)/gi(I00,i) * u0;

        // Calculate cell-centered magnetic fields given Minkowski values
        Real bcont = 0.0;
        Real bconx = bx;
        Real bcony = by;
        Real bconz = bz;
        Real bcon0, bcon1, bcon2, bcon3;
        pcoord->TransformVectorCell(bcont, bconx, bcony, bconz, k, j, i, &bcon0, &bcon1,
            &bcon2, &bcon3);
        b(IB1,k,j,i) = bcon1 * u0 - bcon0 * u1;
        b(IB2,k,j,i) = bcon2 * u0 - bcon0 * u2;
        b(IB3,k,j,i) = bcon3 * u0 - bcon0 * u3;
      }
    }
  }
  peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, il, iu, jl, ju, kl, ku);

  // Delete auxiliary array
  b.DeleteAthenaArray();
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k = kl; k <= ku+1; ++k) {
      for (int j = jl; j <= ju+1; ++j) {
        for (int i = il; i <= iu+1; ++i) {
          Real ut = 1.0;
          Real ux = 0.0;
          Real uy = 0.0;
          Real uz = 0.0;
          Real bcont = 0.0;
          Real bconx = bx;
          Real bcony = by;
          Real bconz = bz;
          Real u0, u1, u2, u3;
          Real bcon0, bcon1, bcon2, bcon3;
          if (j != ju+1 && k != ku+1) {
            pcoord->TransformVectorFace1(ut, ux, uy, uz, k, j, i, &u0, &u1, &u2, &u3);
            pcoord->TransformVectorFace1(bcont, bconx, bcony, bconz, k, j, i, &bcon0,
                &bcon1, &bcon2, &bcon3);
            pfield->b.x1f(k,j,i) = bcon1 * u0 - bcon0 * u1;
          }
          if (i != iu+1 && k != ku+1) {
            pcoord->TransformVectorFace2(ut, ux, uy, uz, k, j, i, &u0, &u1, &u2, &u3);
            pcoord->TransformVectorFace2(bcont, bconx, bcony, bconz, k, j, i, &bcon0,
                &bcon1, &bcon2, &bcon3);
            pfield->b.x2f(k,j,i) = bcon2 * u0 - bcon0 * u2;
          }
          if (i != iu+1 && j != ju+1) {
            pcoord->TransformVectorFace3(ut, ux, uy, uz, k, j, i, &u0, &u1, &u2, &u3);
            pcoord->TransformVectorFace3(bcont, bconx, bcony, bconz, k, j, i, &bcon0,
                &bcon1, &bcon2, &bcon3);
            pfield->b.x3f(k,j,i) = bcon3 * u0 - bcon0 * u3;
          }
        }
      }
    }
  }
  return;
}
