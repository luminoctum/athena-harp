#include <cmath>  // sin, cos, pow
#include "../parameter_input.hpp"
#include "../physics/held_suarez_94.hpp"
#include "../athena_math.hpp" // _sqr
#include "../globals.hpp" // Rgas

//! \file held_suarez_94.cpp
//  \brief implementation of held_suarez_94.hpp

HeldSuarez94::HeldSuarez94() :
  tdy(60.), tdz(10.), psrf(1.E5), tsrf(315.), tmin(200.),
  kappa(2./7.), rgas(286.5), grav(9.8),
  sigma_b(0.7), kf(1./86400.), ka(1./40.), ks(1./4.),
  lat(0.), pbot(1.E5), dz(0.)
{}

void HeldSuarez94::LoadInputFile(ParameterInput *pin)
{
  tdy = pin->GetReal("problem", "tdy");
  tdz = pin->GetReal("problem", "tdz");
  psrf = pin->GetReal("problem", "psrf");
  tsrf = pin->GetReal("problem", "tsrf");
  tmin = pin->GetReal("problem", "tmin");

  kappa = (pin->GetReal("hydro", "gamma") - 1.) / pin->GetReal("hydro", "gamma");
  rgas = Globals::Rgas / pin->GetReal("hydro", "mu");
  grav = - pin->GetReal("hydro", "grav_acc1");

  sigma_b = pin->GetReal("problem", "sigma_b");
  kf = pin->GetReal("problem", "kf");
  ka = pin->GetReal("problem", "ka");
  ks = pin->GetReal("problem", "ks");
}

Real HeldSuarez94::GetTempEq(Real theta, Real pres) const
{
  Real temp = (tsrf - tdy * _sqr(sin(theta)) - tdz * log(pres/psrf) * _sqr(cos(theta))) 
    * pow(pres/psrf, kappa);
  return _max(tmin, temp);
}

Real HeldSuarez94::operator()(Real ptop) const
{
  Real temp1 = GetTempEq(lat, pbot),
       temp2 = GetTempEq(lat, ptop);

  Real rho1 = pbot / (rgas * temp1),
       rho2 = ptop / (rgas * temp2);

  return (ptop - pbot)/dz + sqrt(rho1 * rho2) * grav;
}

Real HeldSuarez94::Kv(Real sigma) const {
  return kf * _max(0., (sigma - sigma_b) / (1. - sigma_b));
}

Real HeldSuarez94::Kt(Real theta, Real sigma) const {
  return ka + (ks - ka) * _max(0., (sigma - sigma_b) / (1. - sigma_b))
    * _sqr(_sqr(cos(theta)));
}
