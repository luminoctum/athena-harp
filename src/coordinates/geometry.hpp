#ifndef GEOMETRY_HPP_
#define GEOMETRY_HPP_
#include "../athena.hpp"

void RotatePoleToEquator(Real *theta1, Real *phi1, Real theta0, Real phi0);
void RotateEquatorToPole(Real *theta1, Real *phi1, Real theta0, Real phi0);
void SphericalLatlonToCartesian(Real *x, Real *y, Real *z, Real a, Real b, Real c, Real phi, Real theta);
void CartesianToSphericalLatlon(Real *a, Real *b, Real *c, Real x, Real y, Real z, Real phi, Real theta);

#endif
