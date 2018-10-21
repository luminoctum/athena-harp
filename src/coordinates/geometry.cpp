#include <cmath>
#include "geometry.hpp"

void RotatePoleToEquator(Real *theta1, Real *phi1, Real theta0, Real phi0)
{
  *theta1 = asin(cos(theta0)*sin(phi0));
  *phi1   = asin(cos(theta0)*cos(phi0)/cos(*theta1));

  if (theta0 < 0. && (phi0 > -M_PI/2. && phi0 < M_PI/2.))
    *phi1 = M_PI - *phi1;
  if (theta0 < 0. && (phi0 < -M_PI/2. || phi0 > M_PI/2.))
    *phi1 = - M_PI - *phi1;
}

void RotateEquatorToPole(Real *theta1, Real *phi1, Real theta0, Real phi0)
{
  *theta1 = asin(cos(theta0)*cos(phi0));
  *phi1   = asin(sin(theta0)/cos(*theta1));

  if (phi0 < 0. && theta0 > 0.)
    *phi1 = M_PI - *phi1;
  if (phi0 < 0. && theta0 < 0.)
    *phi1 = - M_PI - *phi1;
}

void SphericalLatlonToCartesian(
    Real *x, Real *y, Real *z, 
    Real a, Real b, Real c,
    Real phi, Real theta)
{
  *x = -a*sin(phi) - b*sin(theta)*cos(phi) + c*cos(theta)*cos(phi);
  *y = a*cos(phi) - b*sin(theta)*sin(phi) + c*cos(theta)*sin(phi);
  *z = b*cos(theta) + c*sin(theta);
}

void CartesianToSphericalLatlon(
    Real *a, Real *b, Real *c, 
    Real x, Real y, Real z, 
    Real phi, Real theta)
{
  *a = -x*sin(phi) + y*cos(phi);
  *b = -x*sin(theta)*cos(phi) - y*sin(theta)*sin(phi) + z*cos(theta);
  *c = x*cos(theta)*cos(phi) + y*cos(theta)*sin(phi) + z*sin(theta);
}
