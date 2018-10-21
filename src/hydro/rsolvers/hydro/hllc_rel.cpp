//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hllc_rel.cpp
//  \brief Implements HLLC Riemann solver for relativistic hydrodynamics.

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // abs(), sqrt()

// Athena++ headers
#include "../../hydro.hpp"
#include "../../../athena.hpp"                   // enums, macros
#include "../../../athena_arrays.hpp"            // AthenaArray
#include "../../../coordinates/coordinates.hpp"  // Coordinates
#include "../../../eos/eos.hpp"                  // EquationOfState
#include "../../../mesh/mesh.hpp"                // MeshBlock

// Declarations
static void HLLCTransforming(MeshBlock *pmb, const int k, const int j, const int il,
    const int iu, const int ivx, const AthenaArray<Real> &bb,
    AthenaArray<Real> &bb_normal, AthenaArray<Real> &g, AthenaArray<Real> &gi,
    AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r, AthenaArray<Real> &cons,
    AthenaArray<Real> &flux);
static void HLLENonTransforming(MeshBlock *pmb, const int k, const int j, const int il,
    const int iu, AthenaArray<Real> &g, AthenaArray<Real> &gi, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &flux);

//----------------------------------------------------------------------------------------
// Riemann solver
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   bb: 3D array of normal magnetic fields (not used)
//   prim_l, prim_r: left and right primitive states
// Outputs:
//   flux: fluxes across interface
// Notes:
//   prim_l, prim_r overwritten
//   tries to implement HLLC algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB2005)
//   otherwise implements HLLE algorithm similar to that of fluxcalc() in step_ch.c in
//       Harm

void Hydro::RiemannSolver(const int k, const int j, const int il, const int iu,
    const int ivx, const AthenaArray<Real> &bb, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &flux)
{
  if (GENERAL_RELATIVITY and ivx == IVY and pmy_block->pcoord->IsPole(j)) {
    HLLENonTransforming(pmy_block, k, j, il, iu, g_, gi_, prim_l, prim_r, flux);
  } else {
    HLLCTransforming(pmy_block, k, j, il, iu, ivx, bb, bb_normal_, g_, gi_, prim_l,
        prim_r, cons_, flux);
  }
  return;
}

//----------------------------------------------------------------------------------------
// Frame-transforming HLLC implementation
// Inputs:
//   pmb: pointer to MeshBlock object
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   bb: 3D array of normal magnetic fields (not used)
//   bb_normal: 1D scratch array for normal magnetic fields
//   g,gi: 1D scratch arrays for metric coefficients
//   prim_l, prim_r: left and right primitive states
//   cons: 1D scratch array for conserved quantities
// Outputs:
//   flux: fluxes across interface
// Notes:
//   prim_l, prim_r overwritten
//   implements HLLC algorithm from Mignone & Bodo 2005, MNRAS 364 126 (MB2005)
//   references Mignone & Bodo 2006, MNRAS 368 1040 (MB2006)

static void HLLCTransforming(MeshBlock *pmb, const int k, const int j, const int il,
    const int iu, const int ivx, const AthenaArray<Real> &bb,
    AthenaArray<Real> &bb_normal, AthenaArray<Real> &g, AthenaArray<Real> &gi,
    AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r, AthenaArray<Real> &cons,
    AthenaArray<Real> &flux)
{
  // Calculate metric if in GR
  int i01, i11;
  #if GENERAL_RELATIVITY
  {
    switch (ivx) {
      case IVX:
        pmb->pcoord->Face1Metric(k, j, il, iu, g, gi);
        i01 = I01;
        i11 = I11;
        break;
      case IVY:
        pmb->pcoord->Face2Metric(k, j, il, iu, g, gi);
        i01 = I02;
        i11 = I22;
        break;
      case IVZ:
        pmb->pcoord->Face3Metric(k, j, il, iu, g, gi);
        i01 = I03;
        i11 = I33;
        break;
    }
  }
  #endif  // GENERAL_RELATIVITY

  // Transform primitives to locally flat coordinates if in GR
  #if GENERAL_RELATIVITY
  {
    switch (ivx) {
      case IVX:
        pmb->pcoord->PrimToLocal1(k, j, il, iu, bb, prim_l, prim_r, bb_normal);
        break;
      case IVY:
        pmb->pcoord->PrimToLocal2(k, j, il, iu, bb, prim_l, prim_r, bb_normal);
        break;
      case IVZ:
        pmb->pcoord->PrimToLocal3(k, j, il, iu, bb, prim_l, prim_r, bb_normal);
        break;
    }
  }
  #endif  // GENERAL_RELATIVITY

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmb->peos->GetGamma();

  // Go through each interface
  #pragma simd
  for (int i = il; i <= iu; ++i) {

    // Extract left primitives
    const Real &rho_l = prim_l(IDN,i);
    const Real &pgas_l = prim_l(IPR,i);
    Real u_l[4];
    if (GENERAL_RELATIVITY) {
      u_l[1] = prim_l(ivx,i);
      u_l[2] = prim_l(ivy,i);
      u_l[3] = prim_l(ivz,i);
      u_l[0] = std::sqrt(1.0 + SQR(u_l[1]) + SQR(u_l[2]) + SQR(u_l[3]));
    } else {  // SR
      const Real &vx_l = prim_l(ivx,i);
      const Real &vy_l = prim_l(ivy,i);
      const Real &vz_l = prim_l(ivz,i);
      u_l[0] = std::sqrt(1.0 / (1.0 - SQR(vx_l) - SQR(vy_l) - SQR(vz_l)));
      u_l[1] = u_l[0] * vx_l;
      u_l[2] = u_l[0] * vy_l;
      u_l[3] = u_l[0] * vz_l;
    }

    // Extract right primitives
    const Real &rho_r = prim_r(IDN,i);
    const Real &pgas_r = prim_r(IPR,i);
    Real u_r[4];
    if (GENERAL_RELATIVITY) {
      u_r[1] = prim_r(ivx,i);
      u_r[2] = prim_r(ivy,i);
      u_r[3] = prim_r(ivz,i);
      u_r[0] = std::sqrt(1.0 + SQR(u_r[1]) + SQR(u_r[2]) + SQR(u_r[3]));
    } else {  // SR
      const Real &vx_r = prim_r(ivx,i);
      const Real &vy_r = prim_r(ivy,i);
      const Real &vz_r = prim_r(ivz,i);
      u_r[0] = std::sqrt(1.0 / (1.0 - SQR(vx_r) - SQR(vy_r) - SQR(vz_r)));
      u_r[1] = u_r[0] * vx_r;
      u_r[2] = u_r[0] * vy_r;
      u_r[3] = u_r[0] * vz_r;
    }

    // Calculate wavespeeds in left state (MB2005 23)
    Real lambda_p_l, lambda_m_l;
    Real wgas_l = rho_l + gamma_adi/(gamma_adi-1.0) * pgas_l;
    pmb->peos->SoundSpeedsSR(wgas_l, pgas_l, u_l[1]/u_l[0], SQR(u_l[0]), &lambda_p_l,
        &lambda_m_l);

    // Calculate wavespeeds in right state (MB2005 23)
    Real lambda_p_r, lambda_m_r;
    Real wgas_r = rho_r + gamma_adi/(gamma_adi-1.0) * pgas_r;
    pmb->peos->SoundSpeedsSR(wgas_r, pgas_r, u_r[1]/u_r[0], SQR(u_r[0]), &lambda_p_r,
        &lambda_m_r);

    // Calculate extremal wavespeeds
    Real lambda_l = std::min(lambda_m_l, lambda_m_r);
    Real lambda_r = std::max(lambda_p_l, lambda_p_r);

    // Calculate conserved quantities in L region (MB2005 3)
    Real cons_l[NWAVE];
    cons_l[IDN] = rho_l * u_l[0];
    cons_l[IEN] = wgas_l * u_l[0] * u_l[0] - pgas_l;
    cons_l[ivx] = wgas_l * u_l[1] * u_l[0];
    cons_l[ivy] = wgas_l * u_l[2] * u_l[0];
    cons_l[ivz] = wgas_l * u_l[3] * u_l[0];

    // Calculate fluxes in L region (MB2005 2,3)
    Real flux_l[NWAVE];
    flux_l[IDN] = rho_l * u_l[1];
    flux_l[IEN] = wgas_l * u_l[0] * u_l[1];
    flux_l[ivx] = wgas_l * u_l[1] * u_l[1] + pgas_l;
    flux_l[ivy] = wgas_l * u_l[2] * u_l[1];
    flux_l[ivz] = wgas_l * u_l[3] * u_l[1];

    // Calculate conserved quantities in R region (MB2005 3)
    Real cons_r[NWAVE];
    cons_r[IDN] = rho_r * u_r[0];
    cons_r[IEN] = wgas_r * u_r[0] * u_r[0] - pgas_r;
    cons_r[ivx] = wgas_r * u_r[1] * u_r[0];
    cons_r[ivy] = wgas_r * u_r[2] * u_r[0];
    cons_r[ivz] = wgas_r * u_r[3] * u_r[0];

    // Calculate fluxes in R region (MB2005 2,3)
    Real flux_r[NWAVE];
    flux_r[IDN] = rho_r * u_r[1];
    flux_r[IEN] = wgas_r * u_r[0] * u_r[1];
    flux_r[ivx] = wgas_r * u_r[1] * u_r[1] + pgas_r;
    flux_r[ivy] = wgas_r * u_r[2] * u_r[1];
    flux_r[ivz] = wgas_r * u_r[3] * u_r[1];

    // Calculate conserved quantities in HLL region in GR (MB2005 9)
    Real cons_hll[NWAVE];
    for (int n = 0; n < NWAVE; ++n) {
      cons_hll[n] = (lambda_r*cons_r[n] - lambda_l*cons_l[n] + flux_l[n] - flux_r[n])
          / (lambda_r-lambda_l);
    }

    // Calculate fluxes in HLL region (MB2005 11)
    Real flux_hll[NWAVE];
    for (int n = 0; n < NWAVE; ++n) {
      flux_hll[n] = (lambda_r*flux_l[n] - lambda_l*flux_r[n]
          + lambda_l*lambda_r * (cons_r[n] - cons_l[n])) / (lambda_r-lambda_l);
    }

    // Calculate contact wavespeed (MB2005 18)
    Real lambda_star;
    if (std::abs(flux_hll[IEN]) > TINY_NUMBER) {  // use quadratic formula
      // Follows algorithm in Numerical Recipes (section 5.6) for avoiding cancellations
      Real a = flux_hll[IEN];
      Real b = -(cons_hll[IEN] + flux_hll[ivx]);
      Real c = cons_hll[ivx];
      Real q = -0.5 * (b - std::sqrt(SQR(b) - 4.0*a*c));
      lambda_star = c / q;
    } else {  // no quadratic term
      lambda_star = cons_hll[ivx] / (cons_hll[IEN] + flux_hll[ivx]);
    }

    // Calculate contact pressure (MB2006 48)
    // Note: Could also use (MB2005 17), but note the first minus sign there is wrong.
    Real pgas_star = -flux_hll[IEN] * lambda_star + flux_hll[ivx];

    // Calculate conserved quantities in L* region (MB2005 16)
    Real cons_lstar[NWAVE];
    Real vx_l = u_l[1] / u_l[0];
    for (int n = 0; n < NWAVE; ++n) {
      cons_lstar[n] = cons_l[n] * (lambda_l-vx_l);
    }
    cons_lstar[IEN] += pgas_star*lambda_star - pgas_l*vx_l;
    cons_lstar[ivx] += pgas_star - pgas_l;
    for (int n = 0; n < NWAVE; ++n) {
      cons_lstar[n] /= lambda_l - lambda_star;
    }

    // Calculate fluxes in L* region (MB2005 14)
    Real flux_lstar[NWAVE];
    for (int n = 0; n < NWAVE; ++n) {
      flux_lstar[n] = flux_l[n] + lambda_l * (cons_lstar[n] - cons_l[n]);
    }

    // Calculate conserved quantities in R* region (MB2005 16)
    Real cons_rstar[NWAVE];
    Real vx_r = u_r[1] / u_r[0];
    for (int n = 0; n < NWAVE; ++n) {
      cons_rstar[n] = cons_r[n] * (lambda_r-vx_r);
    }
    cons_rstar[IEN] += pgas_star*lambda_star - pgas_r*vx_r;
    cons_rstar[ivx] += pgas_star - pgas_r;
    for (int n = 0; n < NWAVE; ++n) {
      cons_rstar[n] /= lambda_r - lambda_star;
    }

    // Calculate fluxes in R* region (MB2005 14)
    Real flux_rstar[NWAVE];
    for (int n = 0; n < NWAVE; ++n) {
      flux_rstar[n] = flux_r[n] + lambda_r * (cons_rstar[n] - cons_r[n]);
    }

    // Calculate interface velocity
    Real v_interface = 0.0;
    if (GENERAL_RELATIVITY) {
      v_interface = gi(i01,i) / std::sqrt(SQR(gi(i01,i)) - gi(I00,i)*gi(i11,i));
    }

    // Set conserved quantities in GR
    if (GENERAL_RELATIVITY) {
      for (int n = 0; n < NWAVE; ++n) {
        if (lambda_l >= v_interface) {  // L region
          cons(n,i) = cons_l[n];
        } else if (lambda_r <= v_interface) {  // R region
          cons(n,i) = cons_r[n];
        } else if (lambda_star >= v_interface) {  // L* region
          cons(n,i) = cons_lstar[n];
        } else {  // R* region
          cons(n,i) = cons_rstar[n];
        }
      }
    }

    // Set fluxes
    for (int n = 0; n < NWAVE; ++n) {
      if (lambda_l >= v_interface) {  // L region
        flux(n,i) = flux_l[n];
      } else if (lambda_r <= v_interface) {  // R region
        flux(n,i) = flux_r[n];
      } else if (lambda_star >= v_interface) {  // L* region
        flux(n,i) = flux_lstar[n];
      } else {  // R* region
        flux(n,i) = flux_rstar[n];
      }
    }
  }

  // Transform fluxes to global coordinates if in GR
  #if GENERAL_RELATIVITY
  {
    switch (ivx) {
      case IVX:
        pmb->pcoord->FluxToGlobal1(k, j, il, iu, cons, bb_normal, flux);
        break;
      case IVY:
        pmb->pcoord->FluxToGlobal2(k, j, il, iu, cons, bb_normal, flux);
        break;
      case IVZ:
        pmb->pcoord->FluxToGlobal3(k, j, il, iu, cons, bb_normal, flux);
        break;
    }
  }
  #endif  // GENERAL_RELATIVITY
  return;
}

//----------------------------------------------------------------------------------------
// Non-frame-transforming HLLE implementation
// Inputs:
//   pmb: pointer to MeshBlock object
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   g,gi: 1D scratch arrays for metric coefficients
//   prim_l, prim_r: left and right primitive states
// Outputs:
//   flux: fluxes across interface
// Notes:
//   implements HLLE algorithm similar to that of fluxcalc() in step_ch.c in Harm
//   derived from RiemannSolver() in hlle_rel_no_transform.cpp assuming ivx = IVY
//   same function as in hlle_rel.cpp

static void HLLENonTransforming(MeshBlock *pmb, const int k, const int j, const int il,
    const int iu, AthenaArray<Real> &g, AthenaArray<Real> &gi, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &flux)
#if GENERAL_RELATIVITY
{
  // Extract ratio of specific heats
  const Real gamma_adi = pmb->peos->GetGamma();

  // Get metric components
  pmb->pcoord->Face2Metric(k, j, il, iu, g, gi);

  // Go through each interface
  #pragma simd
  for (int i = il; i <= iu; ++i) {

    // Extract metric
    const Real &g_00 = g(I00,i), &g_01 = g(I01,i), &g_02 = g(I02,i), &g_03 = g(I03,i),
               &g_10 = g(I01,i), &g_11 = g(I11,i), &g_12 = g(I12,i), &g_13 = g(I13,i),
               &g_20 = g(I02,i), &g_21 = g(I12,i), &g_22 = g(I22,i), &g_23 = g(I23,i),
               &g_30 = g(I03,i), &g_31 = g(I13,i), &g_32 = g(I23,i), &g_33 = g(I33,i);
    const Real &g00 = gi(I00,i), &g01 = gi(I01,i), &g02 = gi(I02,i), &g03 = gi(I03,i),
               &g10 = gi(I01,i), &g11 = gi(I11,i), &g12 = gi(I12,i), &g13 = gi(I13,i),
               &g20 = gi(I02,i), &g21 = gi(I12,i), &g22 = gi(I22,i), &g23 = gi(I23,i),
               &g30 = gi(I03,i), &g31 = gi(I13,i), &g32 = gi(I23,i), &g33 = gi(I33,i);
    Real alpha = std::sqrt(-1.0/g00);

    // Extract left primitives
    const Real &rho_l = prim_l(IDN,i);
    const Real &pgas_l = prim_l(IPR,i);
    const Real &uu1_l = prim_l(IVX,i);
    const Real &uu2_l = prim_l(IVY,i);
    const Real &uu3_l = prim_l(IVZ,i);

    // Extract right primitives
    const Real &rho_r = prim_r(IDN,i);
    const Real &pgas_r = prim_r(IPR,i);
    const Real &uu1_r = prim_r(IVX,i);
    const Real &uu2_r = prim_r(IVY,i);
    const Real &uu3_r = prim_r(IVZ,i);

    // Calculate 4-velocity in left state
    Real ucon_l[4], ucov_l[4];
    Real tmp = g_11*SQR(uu1_l) + 2.0*g_12*uu1_l*uu2_l + 2.0*g_13*uu1_l*uu3_l
             + g_22*SQR(uu2_l) + 2.0*g_23*uu2_l*uu3_l
             + g_33*SQR(uu3_l);
    Real gamma_l = std::sqrt(1.0 + tmp);
    ucon_l[0] = gamma_l / alpha;
    ucon_l[1] = uu1_l - alpha * gamma_l * g01;
    ucon_l[2] = uu2_l - alpha * gamma_l * g02;
    ucon_l[3] = uu3_l - alpha * gamma_l * g03;
    ucov_l[0] = g_00*ucon_l[0] + g_01*ucon_l[1] + g_02*ucon_l[2] + g_03*ucon_l[3];
    ucov_l[1] = g_10*ucon_l[0] + g_11*ucon_l[1] + g_12*ucon_l[2] + g_13*ucon_l[3];
    ucov_l[2] = g_20*ucon_l[0] + g_21*ucon_l[1] + g_22*ucon_l[2] + g_23*ucon_l[3];
    ucov_l[3] = g_30*ucon_l[0] + g_31*ucon_l[1] + g_32*ucon_l[2] + g_33*ucon_l[3];

    // Calculate 4-velocity in right state
    Real ucon_r[4], ucov_r[4];
    tmp = g_11*SQR(uu1_r) + 2.0*g_12*uu1_r*uu2_r + 2.0*g_13*uu1_r*uu3_r
        + g_22*SQR(uu2_r) + 2.0*g_23*uu2_r*uu3_r
        + g_33*SQR(uu3_r);
    Real gamma_r = std::sqrt(1.0 + tmp);
    ucon_r[0] = gamma_r / alpha;
    ucon_r[1] = uu1_r - alpha * gamma_r * g01;
    ucon_r[2] = uu2_r - alpha * gamma_r * g02;
    ucon_r[3] = uu3_r - alpha * gamma_r * g03;
    ucov_r[0] = g_00*ucon_r[0] + g_01*ucon_r[1] + g_02*ucon_r[2] + g_03*ucon_r[3];
    ucov_r[1] = g_10*ucon_r[0] + g_11*ucon_r[1] + g_12*ucon_r[2] + g_13*ucon_r[3];
    ucov_r[2] = g_20*ucon_r[0] + g_21*ucon_r[1] + g_22*ucon_r[2] + g_23*ucon_r[3];
    ucov_r[3] = g_30*ucon_r[0] + g_31*ucon_r[1] + g_32*ucon_r[2] + g_33*ucon_r[3];

    // Calculate wavespeeds in left state
    Real lambda_p_l, lambda_m_l;
    Real wgas_l = rho_l + gamma_adi/(gamma_adi-1.0) * pgas_l;
    pmb->peos->SoundSpeedsGR(wgas_l, pgas_l, ucon_l[0], ucon_l[IVY], g00, g02, g22,
        &lambda_p_l, &lambda_m_l);

    // Calculate wavespeeds in right state
    Real lambda_p_r, lambda_m_r;
    Real wgas_r = rho_r + gamma_adi/(gamma_adi-1.0) * pgas_r;
    pmb->peos->SoundSpeedsGR(wgas_r, pgas_r, ucon_r[0], ucon_r[IVY], g00, g02, g22,
        &lambda_p_r, &lambda_m_r);

    // Calculate extremal wavespeeds
    Real lambda_l = std::min(lambda_m_l, lambda_m_r);
    Real lambda_r = std::max(lambda_p_l, lambda_p_r);

    // Calculate conserved quantities in L region (rho u^0 and T^0_\mu)
    Real cons_l[NWAVE];
    cons_l[IDN] = rho_l * ucon_l[0];
    cons_l[IEN] = wgas_l * ucon_l[0] * ucov_l[0] + pgas_l;
    cons_l[IVX] = wgas_l * ucon_l[0] * ucov_l[1];
    cons_l[IVY] = wgas_l * ucon_l[0] * ucov_l[2];
    cons_l[IVZ] = wgas_l * ucon_l[0] * ucov_l[3];

    // Calculate fluxes in L region (rho u^i and T^i_\mu, where i = IVY)
    Real flux_l[NWAVE];
    flux_l[IDN] = rho_l * ucon_l[IVY];
    flux_l[IEN] = wgas_l * ucon_l[IVY] * ucov_l[0];
    flux_l[IVX] = wgas_l * ucon_l[IVY] * ucov_l[1];
    flux_l[IVY] = wgas_l * ucon_l[IVY] * ucov_l[2];
    flux_l[IVZ] = wgas_l * ucon_l[IVY] * ucov_l[3];
    flux_l[IVY] += pgas_l;

    // Calculate conserved quantities in R region (rho u^0 and T^0_\mu)
    Real cons_r[NWAVE];
    cons_r[IDN] = rho_r * ucon_r[0];
    cons_r[IEN] = wgas_r * ucon_r[0] * ucov_r[0] + pgas_r;
    cons_r[IVX] = wgas_r * ucon_r[0] * ucov_r[1];
    cons_r[IVY] = wgas_r * ucon_r[0] * ucov_r[2];
    cons_r[IVZ] = wgas_r * ucon_r[0] * ucov_r[3];

    // Calculate fluxes in R region (rho u^i and T^i_\mu, where i = IVY)
    Real flux_r[NWAVE];
    flux_r[IDN] = rho_r * ucon_r[IVY];
    flux_r[IEN] = wgas_r * ucon_r[IVY] * ucov_r[0];
    flux_r[IVX] = wgas_r * ucon_r[IVY] * ucov_r[1];
    flux_r[IVY] = wgas_r * ucon_r[IVY] * ucov_r[2];
    flux_r[IVZ] = wgas_r * ucon_r[IVY] * ucov_r[3];
    flux_r[IVY] += pgas_r;

    // Calculate fluxes in HLL region
    Real flux_hll[NWAVE];
    for (int n = 0; n < NWAVE; ++n) {
      flux_hll[n] = (lambda_r*flux_l[n] - lambda_l*flux_r[n]
          + lambda_r*lambda_l * (cons_r[n] - cons_l[n])) / (lambda_r-lambda_l);
    }

    // Set fluxes
    for (int n = 0; n < NWAVE; ++n) {
      if (lambda_l >= 0.0) {  // L region
        flux(n,i) = flux_l[n];
      } else if (lambda_r <= 0.0) {  // R region
        flux(n,i) = flux_r[n];
      } else {  // HLL region
        flux(n,i) = flux_hll[n];
      }
    }
  }
  return;
}

#else
{
  return;
}
#endif  // GENERAL_RELATIVITY
