//! \file  roe_shallow_water.cpp
//  \brief Roe's linearized Riemann solver for shallow water model

// C/C++ headers
#include <cmath>      // sqrt()
#include <iostream>
#include <algorithm>  // min()

// Athena++ headers
#include "../hydro.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../eos/eos.hpp"
#include "../../math_funcs.hpp"
#include "../../globals.hpp"

inline Real Enthalpy(Real temp, Real const rhon[], Real const mu[], Real const cv[], Real const latent[])
{
  Real cp = 0, lt = 0., rho = 0.;
  for (int n = 0; n < NGAS; ++n) {
    cp += (cv[n] + Globals::Rgas)/mu[n]*rhon[n];
    rho += rhon[n];
  }
  for (int n = NGAS; n < NCOMP; ++n) {
    cp += cv[n]/mu[n]*rhon[n];
    lt += latent[n]/mu[n]*rhon[n];
    rho += rhon[n];
  }
  return (cp*temp + lt)/rho;
}

inline Real TotalDensity(AthenaArray<Real> &w, int i, Real const mu[], Real rhon[])
{
  Real xt = 1., rho = 0.;
  for (int n = NGAS; n < NCOMP; ++n)
    xt -= w(n,i);
  Real x1 = xt;
  for (int n = 1; n < NGAS; ++n) {
    rhon[n] = w(n,i)/xt*w(IPR,i)*mu[n]/(Globals::Rgas*w(IT,i));
    x1 -= w(n,i);
    rho += rhon[n];
  }
  rhon[0] = x1/xt*w(IPR,i)*mu[0]/(Globals::Rgas*w(IT,i));
  rho += rhon[0];
  for (int n = NGAS; n < NCOMP; ++n) {
    rhon[n] = (w(n,i)*mu[n])/(x1*mu[0])*rhon[0];
    rho += rhon[n];
  }
  return rho;
}

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
    int const ivx, AthenaArray<Real> const& bx, AthenaArray<Real> &wl,
    AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  EquationOfState *peos = pmy_block->peos;
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  Real wave[3][NHYDRO], speed[3];
  Real ubar, vbar, wbar, cbar, hbar, hl, hr, rhobar;
  Real qbar[NCOMP], rholn[NCOMP], rhorn[NCOMP];
  Real alpha[NHYDRO];
  Real rhol, rhor, sqrhol, sqrhor, isqrho;
  Real du, dv, dw, dp;

  for (int i = il; i <= iu; ++i) {
    rhol = TotalDensity(wl, i, peos->mu_, rholn);
    rhor = TotalDensity(wr, i, peos->mu_, rhorn);
    sqrhol = sqrt(rhol);
    sqrhor = sqrt(rhor);
    isqrho = 1./(sqrhol + sqrhor);
    rhobar = sqrhol*sqrhor;

    // u,v,w
    ubar = isqrho*(wl(ivx,i)*sqrhol + wr(ivx,i)*sqrhor);
    vbar = isqrho*(wl(ivy,i)*sqrhol + wr(ivy,i)*sqrhor);
    wbar = isqrho*(wl(ivz,i)*sqrhol + wr(ivz,i)*sqrhor);

    hl = Enthalpy(wl(IT,i), rholn, peos->mu_, peos->cv_, peos->latent_)
        + 0.5*(_sqr(wl(IVX,i)) + _sqr(wl(IVY,i)) + _sqr(wl(IVZ,i)));
    hr = Enthalpy(wr(IT,i), rhorn, peos->mu_, peos->cv_, peos->latent_)
        + 0.5*(_sqr(wr(IVX,i)) + _sqr(wr(IVY,i)) + _sqr(wr(IVZ,i)));
    hbar = isqrho*(hl*sqrhol + hr*sqrhor);

    // qbar, mass mixing ratio
    for (int n = 0; n < NCOMP; ++n)
      qbar[n] = isqrho*(rholn[n]/sqrhol + rhorn[n]/sqrhor);

    // kbar = rck/rc
    Real rc = 0., rck = 0.;
    for (int n = 0; n < NGAS; ++n) {
      rc += peos->cv_[n]/peos->mu_[n]*qbar[n];
      rck += Globals::Rgas/peos->mu_[n]*qbar[n];
    }
    for (int n = NGAS; n < NCOMP; ++n)
      rc += peos->cv_[n]/peos->mu_[n]*qbar[n];

    // cbar
    Real latent = 0.;
    for (int n = NGAS; n < NCOMP; ++n)
      latent += peos->latent_[n]/peos->mu_[n]*qbar[n];
    cbar = sqrt(rck/rc*(hbar - 0.5*(_sqr(ubar) + _sqr(vbar) + _sqr(wbar)) - latent));

    // primitive variable difference
    // using low Mach number fix (Rieper, 2011)
    //du = (wr(ivx,i) - wl(ivx,i))*
    //  std::min((fabs(ubar) + fabs(vbar) + fabs(wbar))/cbar, 1.);
    du = wr(ivx,i) - wl(ivx,i);
    dv = wr(ivy,i) - wl(ivy,i);
    dw = wr(ivz,i) - wl(ivz,i);
    dp = wr(IPR,i) - wl(IPR,i);

    // coefficients of eigenvectors
    alpha[ivx] = - 0.5*rhobar/cbar*du + 0.5*dp/_sqr(cbar);
    alpha[IPR] = + 0.5*rhobar/cbar*du + 0.5*dp/_sqr(cbar);
    alpha[ivy] = rhobar*dv;
    alpha[ivz] = rhobar*dw;
    for (int n = 0; n < NCOMP; ++n)
      alpha[n] = rhorn[n] - rholn[n] - dp/_sqr(cbar)*qbar[n];

    // wave 1, u-c
    for (int n = 0; n < NCOMP; ++n)
      wave[0][n] = alpha[ivx]*qbar[n];
    wave[0][ivx] = alpha[ivx]*(ubar - cbar);
    wave[0][ivy] = alpha[ivx]*vbar;
    wave[0][ivz] = alpha[ivx]*wbar;
    wave[0][IPR] = alpha[ivx]*(hbar - ubar*cbar);

    // wave 3, u+c
    for (int n = 0; n < NCOMP; ++n)
      wave[2][n] = alpha[IPR]*qbar[n];
    wave[2][ivx] = alpha[IPR]*(ubar + cbar);
    wave[2][ivy] = alpha[IPR]*vbar;
    wave[2][ivz] = alpha[IPR]*wbar;
    wave[2][IPR] = alpha[IPR]*(hbar + ubar*cbar);

    // wave 2, u
    Real alpha0 = 0.;
    for (int n = 0; n < NCOMP; ++n) {
      wave[1][n] = alpha[n];
      alpha0 += alpha[n];
    }
    wave[1][ivx] = ubar*alpha0;
    wave[1][ivy] = vbar*alpha0 + alpha[ivy];
    wave[1][ivz] = wbar*alpha0 + alpha[ivz];
    wave[1][IPR] = rhobar*(hr - hl - ubar*du) + alpha0*hbar - dp;

    // speed
    speed[0] = fabs(ubar - cbar);
    speed[1] = fabs(ubar);
    speed[2] = fabs(ubar + cbar);

    // flux
    for (int n = 0; n < NCOMP; ++n)
      flx(n,i) = 0.5*(rholn[n]*wl(ivx,i) + rhorn[n]*wr(ivx,i));
    flx(ivx,i) = 0.5*(rhol*_sqr(wl(ivx,i)) + wl(IPR,i)
                    + rhor*_sqr(wr(ivx,i)) + wr(IPR,i));
    flx(ivy,i) = 0.5*(rhol*wl(ivx,i)*wl(ivy,i)
                    + rhor*wr(ivx,i)*wr(ivy,i));
    flx(ivz,i) = 0.5*(rhol*wl(ivx,i)*wl(ivz,i)
                    + rhor*wr(ivx,i)*wr(ivz,i));
    flx(IPR,i) = 0.5*(rhol*hl*wl(ivx,i)
                    + rhor*hr*wr(ivx,i));
    for (int r = 0; r < 3; ++r)
      for (int n = 0; n < NHYDRO; ++n)
        flx(n,i) -= 0.5*speed[r]*wave[r][n];
  }
}
