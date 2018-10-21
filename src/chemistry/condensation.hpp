// C++ header
#include <cmath>

// Athena++ headers
#include "../math_funcs.hpp"
#include "reaction.hpp"
#include "sat_vapor_pres.hpp"

inline Real GasGasSolidNH4SH(Reaction const& rc, Real const prim[], Real alpha = 0.)
{
  Real rate;
  Real xnh3 = prim[rc.reactor[0]];
  Real xh2s = prim[rc.reactor[1]];
  Real xnh4sh = prim[rc.reactor[2]];
  Real s = SatVaporPresNH4SHLewis(prim[IT])/_sqr(prim[IPR]);
  if (s > 1.) return -xnh4sh;
  Real xt = 1.;
  for (int g = NGAS; g < NCOMP; ++g) xt -= prim[g];
  rate = (xnh3+xh2s-4.*s*xt-sqrt(_sqr(xnh3-xh2s)+4.*s*(xt-2.*xnh3)*(xt-2.*xh2s)))/(2.*(1.-4.*s));
  if (rate > 0.) rate = std::min(std::min(rate, xnh3), xh2s);
  if (rate < 0.) rate = std::max(rate, -xnh4sh);
  return rate;
}

// alpha = L/cv
inline Real GasCloudIdeal(Reaction const& rc, Real const prim[], Real alpha = 0.)
{
  Real rate;
  Real tr = rc.coeff[0],
       pr = rc.coeff[1],
       beta = rc.coeff[2],
       gamma = rc.coeff[3];
  Real x1 = prim[rc.reactor[0]];
  Real xc = prim[rc.reactor[1]];
  Real s = SatVaporPresIdeal(prim[IT]/tr,pr,beta,gamma)/prim[IPR];
  if (s > 1.) return -xc;
  Real xt = 1.;
  for (int g = NGAS; g < NCOMP; ++g) xt -= prim[g];
  Real dsdt = s/prim[IT]*(beta/prim[IT]*tr - gamma);
  rate = (x1 - s*xt)/((1. - s)*(1. - alpha*(xt - x1)/_sqr(1. - s)*dsdt));
  //std::cout << rc.reactor[0] << " " << rc.reactor[1] << " " << rate 
  //          << " " << x1 << " " << s << std::endl;
  if (rate > 0.) rate = std::min(rate, x1);
  if (rate < 0.) rate = std::max(rate, -xc);
  return rate;
}

inline Real LiquidSolidIdeal(Reaction const& rc, Real const prim[], Real alpha = 0.)
{
  Real rate;
  Real tr = rc.coeff[0];
  if (prim[IT] < tr)
    rate = prim[rc.reactor[0]];
  else if (prim[IT] > tr)
    rate = - prim[rc.reactor[1]];
  return rate;
}
