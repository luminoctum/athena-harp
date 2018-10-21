#ifndef SAT_VAPOR_PRES_HPP_
#define SAT_VAPOR_PRES_HPP_

// C++ header
#include <cmath>

// Athena++ headers
#include "../athena.hpp"

inline Real SatVaporPresIdeal(Real t, Real p, Real beta, Real gamma) {
  return p*exp((1. - 1./t)*beta - gamma*log(t));
}

inline Real DlnpDlnTIdeal(Real t, Real beta, Real gamma) {
  return beta/t - gamma;
}

inline Real SatVaporPresH2OUMich(Real T) {
  Real x =
    -2445.5646/T +
    8.2312*log10(T) -
    0.01677006*T +
    1.20514e-05*T*T -
    6.757169;
  return pow(10.0,x) * 1.E5 / 760.;
}

inline Real SatVaporPresH2OAntoine(Real T)
{
  Real result;
  
  if (T < 303.)
    result = pow(10., 5.40221 - (1838.675 / (T - 31.737)));
  else if (T < 333.)
    result = pow(10., 5.20389 - (1733.926 / (T - 39.485)));
  else if (T < 363.)
    result = pow(10., 5.0768 - (1659.793 / (T - 45.854)));
  else
    result = pow(10., 5.08354 - (1663.125 / (T - 45.622)));

  return 1.E5 * result;
}

inline Real SatVaporPresH2OHubner(Real T)
{
  Real A, B ,C, D;
  A = 4.07023;
  B = -2484.98;
  C = 3.56654;
  D = -0.00320981;
  Real x = A + B/T + C*log10(T) + D*T; // best fit in [100K; 273.16K]

  return pow(10.0, x);
}

inline Real SatVaporPresH2OFray(Real T)
{
  Real a[7], x, pt, Tt, tr;
  pt = 6.11657e-03;
  Tt = 273.16;
  tr = T/Tt;
  x = 0.;
  a[0] = 20.9969665107897;
  a[1] = 3.72437478271362;
  a[2] = -13.9205483215524;
  a[3] = 29.6988765013566;
  a[4] = -40.1972392635944;
  a[5] = 29.7880481050215;
  a[6] = -9.13050963547721;

  for(int i = 0; i < 7; i++)
    x = x + a[i]*pow(tr, i);  // best fit in [0K; 273.16K]

  return 1.E5 * pt * exp(1.5*log(tr) + (1 - 1/tr)*x);
}

inline Real SatVaporPresH2OBriggsS(Real T)
{
  Real a[6], x;
  if ( T < 273.16 ) {
    a[1] = -5631.1206;
    a[2] = -8.363602;
    a[3] = 8.2312;
    a[4] = -3.861449e-02;
    a[5] = 2.77494e-05;
  }
  else {
    a[1] = -2313.0338;
    a[2] = -164.03307;
    a[3] = 38.053682;
    a[4] = -1.3844344e-01;
    a[5] = 7.4465367e-05;
  }
  x = a[1]/T + a[2] + a[3]*log(T) + a[4]*T + a[5]*pow(T, 2);
  return exp(x) / 10.;
}

inline Real SatVaporPresH2OBolton(Real T)
{
  Real result;
  result = 612.2 * exp(17.67 * (T - 273.15) / (T - 29.65));

  return result;
}

inline Real SatVaporPresH2OSmithsonian(Real T)
{
  Real result;
  result = 100. * exp(23.33086 - 6111.72784 / T + 0.15215 * log(T));

  return result;
}

inline Real SatVaporPresH2OIdeal(Real T)
{
  Real betal = 24.88,
         gammal = 5.06,
         betas = 22.98,
         gammas = 0.52,
         tr = 273.16,
         pr = 611.7;

  return T > tr ? SatVaporPresIdeal(T / tr, pr, betal, gammal)
    : SatVaporPresIdeal(T / tr, pr, betas, gammas);
}

inline Real SatVaporPresNH3UMich(Real T)
{
  Real x = -1790.00/T - 1.81630*log10(T) + 14.97593;
  return pow(10.0,x) * 1.E5 / 760.;
}

inline Real SatVaporPresNH3Antoine(Real T)
{
  Real result;
  if (T < 239.6)
    result = pow(10., 3.18757 - (506.713 / (T - 80.78)));
  else
    result = pow(10., 4.86886 - (1113.928 / (T - 10.409)));

  return 1.E5 * result;
}

inline Real SatVaporPresNH3Hubner(Real T)
{

  Real A = 24.3037,
       B = -1766.28,
       C = -5.64472,
       D = 0.00740241;

  Real x = A + B/T + C*log10(T) + D*T; // best fit in [130K; 200K]

  return pow(10.0, x);
}

inline Real SatVaporPresNH3BriggsS(Real T)
{
  Real a[6], x;
  if ( T < 195 ) {
    a[1] = -4122.;
    a[2] = 41.67871;
    a[3] = -1.81630;
    a[4] = 0.;
    a[5] = 0.;
  }
  else {
    a[1] = -4409.3512;
    a[2] = 76.864252;
    a[3] = -8.4598340;
    a[4] = 5.51029e-03;
    a[5] = 6.80463e-06;
  }
  x = a[1]/T + a[2] + a[3]*log(T) + a[4]*T + a[5]*pow(T, 2);
  return exp(x) / 10.;
}

inline Real SatVaporPresNH3Fray(Real T)
{
  Real a[7], x = 0;
  a[0] = 1.596e+01;
  a[1] = -3.537e+03;
  a[2] = -3.310e+04;
  a[3] = 1.742e+06;
  a[4] = -2.995e+07;
  a[5] = 0.;
  a[6] = 0.;

  for(int i = 1; i < 7; i++)
    x = x + a[i]/pow(T, i);  // best fit in [15K; 195.41K]

   return 1.E5 * exp(x + a[0]);
}

inline Real SatVaporPresNH3Ideal(Real T)
{
  Real betal = 20.08,
       gammal = 5.62,
       betas = 20.64,
       gammas = 1.43,
       tr = 195.4,
       pr = 6060.;

  return T > tr ? SatVaporPresIdeal(T / tr, pr, betal, gammal)
    : SatVaporPresIdeal(T / tr, pr, betas, gammas);
}

inline Real SatVaporPresH2SAntoine(Real T)
{
  Real result;
  if (T < 212.8)
    result = pow(10., 4.43681 - (829.439 / (T - 25.412)));
  else
    result = pow(10., 4.52887 - (958.587 / (T - 0.539)));

  return 1.E5 * result;
}

inline Real SatVaporPresNH4SHLewis(Real T)
{
  return pow(10., 14.82 - 4705. / T) * 101325. * 101325.;
}

#endif
