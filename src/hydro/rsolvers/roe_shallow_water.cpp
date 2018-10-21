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

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
    int const ivx, AthenaArray<Real> const& bx, AthenaArray<Real> &wl,
    AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = ivx == 1 ? 2 : 1;

  Real wli[3], wri[3], wave[3][3], speed[3];

  Real ubar, vbar, cbar, delh, delu, delv, hbar, a1, a2, a3;

  for (int i = il; i <= iu; ++i) {
    wli[0] = wl(IDN, i);
    wli[1] = wl(ivx, i);
    wli[2] = wl(ivy, i);

    wri[0] = wr(IDN, i);
    wri[1] = wr(ivx, i);
    wri[2] = wr(ivy, i);

    ubar = (wli[1] * sqrt(wli[0]) + wri[1] * sqrt(wri[0])) / (sqrt(wli[0]) + sqrt(wri[0]));
    vbar = (wli[2] * sqrt(wli[0]) + wri[2] * sqrt(wri[0])) / (sqrt(wli[0]) + sqrt(wri[0]));
    cbar = sqrt(0.5 * (wli[0] + wri[0]));

    delh = wri[0] - wli[0];
    //delu = wri[1] - wli[1];
    // low Mach number fix (Rieper, 2011)
    delu = (wri[1] - wli[1]) * std::min((fabs(ubar) + fabs(vbar))/cbar, 1.);
    delv = wri[2] - wli[2];
    hbar = sqrt(wli[0] * wri[0]);

    a1 = 0.5 * (cbar * delh - hbar * delu) / cbar;
    a2 = hbar * delv;
    a3 = 0.5 * (cbar * delh + hbar * delu) / cbar;

    wave[0][0] = a1; wave[0][1] = a1 * (ubar - cbar); wave[0][2] = a1 * vbar;
    wave[1][0] = 0.; wave[1][1] = 0.; wave[1][2] = a2;
    wave[2][0] = a3; wave[2][1] = a3 * (ubar + cbar); wave[2][2] = a3 * vbar;

    speed[0] = fabs(ubar - cbar);
    speed[1] = fabs(ubar);
    speed[2] = fabs(ubar + cbar);

    flx(IDN, i) = 0.5 * (wli[0] * wli[1] + wri[0] * wri[1]);
    flx(ivx, i) = 0.5 * (wli[0] * wli[1] * wli[1] + 0.5 * wli[0] * wli[0]
                      +  wri[0] * wri[1] * wri[1] + 0.5 * wri[0] * wri[0]);
    flx(ivy, i) = 0.5 * (wli[0] * wli[1] * wli[2]
                      +  wri[0] * wri[1] * wri[2]);
    flx(IVZ, i) = 0.;

    for (int r = 0; r < 3; ++r) {
      flx(IDN, i) -= 0.5 * speed[r] * wave[r][0];
      flx(ivx, i) -= 0.5 * speed[r] * wave[r][1];
      flx(ivy, i) -= 0.5 * speed[r] * wave[r][2];
    }
  }
}
