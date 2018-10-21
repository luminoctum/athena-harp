//! \file jupiter_pole.cpp
//  \brief jupiter polar model

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../math_funcs.hpp"     // _root, _sqr
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../coordinates/geometry.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../particle/particle.hpp"

// Coriolis parameter
Real omegax, omegay, omegaz;

Real ZonalWind(Real lat, Real lat1, Real lat2, Real umax)
{
  if (lat <= lat1 || lat >= lat2) return 0.;
  Real en = exp(-4./_sqr(lat2 - lat1));
  return umax/en * exp(1./((lat - lat1)*(lat - lat2)));
}

// integrate using simpson's 3/8 rule
Real Geopotential(Real lat, Real lat1, Real lat2, Real umax, Real omega, Real radius)
{
  if (lat <= lat1) return 0.;
  int num = 100;
  Real ds = (lat - lat1)/num, gh = 0.;
  for (int i = 0; i < num; ++i) {
    Real s1 = lat1 + ds*i;
    Real s2 = lat1 + ds*(i + 1./3.);
    Real s3 = lat1 + ds*(i + 2./3.);
    Real s4 = lat1 + ds*(i + 1);

    Real u1 = ZonalWind(s1, lat1, lat2, umax);
    Real u2 = ZonalWind(s2, lat1, lat2, umax);
    Real u3 = ZonalWind(s3, lat1, lat2, umax);
    Real u4 = ZonalWind(s4, lat1, lat2, umax);

    Real f1 = - 2.*omega*sin(s1)*radius*u1 - u1*u1*tan(s1);
    Real f2 = - 2.*omega*sin(s2)*radius*u2 - u2*u2*tan(s2);
    Real f3 = - 2.*omega*sin(s3)*radius*u3 - u3*u3*tan(s3);
    Real f4 = - 2.*omega*sin(s4)*radius*u4 - u4*u4*tan(s4);

    gh += ds/8.*(f1 + 3.*f2 + 3.*f3 + f4);
  }
  return gh;
}

struct FarFieldGeopotentialSolver {
  Real lat1, lat2, gheq, omega, num, umax, radius;
  Real operator()(Real gh) {
    Real ds = M_PI/num, total_gh = 0.;
    for (int i = 0; i < num; ++i) {
      Real s1 = -M_PI/2. + ds*i;
      Real s2 = -M_PI/2. + ds*(i + 1./3.);
      Real s3 = -M_PI/2. + ds*(i + 2./3.);
      Real s4 = -M_PI/2. + ds*(i + 1);

      Real u1 = ZonalWind(s1, lat1, lat2, umax);
      Real u2 = ZonalWind(s2, lat1, lat2, umax);
      Real u3 = ZonalWind(s3, lat1, lat2, umax);
      Real u4 = ZonalWind(s4, lat1, lat2, umax);

      Real f1 = - 2.*omega*sin(s1)*radius*u1 - u1*u1*tan(s1);
      Real f2 = - 2.*omega*sin(s2)*radius*u2 - u2*u2*tan(s2);
      Real f3 = - 2.*omega*sin(s3)*radius*u3 - u3*u3*tan(s3);
      Real f4 = - 2.*omega*sin(s4)*radius*u4 - u4*u4*tan(s4);

      gh += ds/8.*(f1 + 3.*f2 + 3.*f3 + f4);
      total_gh += gh*cos(-M_PI/2. + 0.5*ds*i)*ds;
    }
    return total_gh/2. - gheq;
  }
};

// particle functions
bool ParticleTranslate(MeshBlock *pmb, Real const time, Real const dt,
  Particle &pt, AthenaArray<Real> const& prim, AthenaArray<Real> const& cons,
  int kji[3], AthenaArray<Real>& cons_out)
{
  Real radius = pmb->pcoord->x3v(0),
       lat = pmb->pcoord->x2v(kji[1]);
  pt.x1 += pt.v1 * dt / (radius * cos(lat));
  pt.x2 += pt.v2 * dt / radius;
  return true;
}

// Coriolis force for spherical latlon grid
void Hydro::SourceTerm(const Real time, const Real dt, const int step,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &cons, 
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons_out)
{
  MeshBlock *pmb = pmy_block;
  Coordinates *pcoord = pmb->pcoord;
  bool nx3 = pmb->pmy_mesh->mesh_size.nx3;

  Real phi, theta, omega1, omega2, omega3;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) 
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        phi = pmb->pcoord->x1v(i);
        theta = pmb->pcoord->x2v(j);

        omega1 = - sin(phi)*omegax + cos(phi)*omegay;
        omega2 = - sin(theta)*cos(phi)*omegax - sin(theta)*sin(phi)*omegay + cos(theta)*omegaz;
        omega3 = cos(theta)*cos(phi)*omegax + cos(theta)*sin(phi)*omegay + sin(theta)*omegaz;

        cons_out(IM1,k,j,i) += 2.*dt*(omega3*cons(IM2,k,j,i) - omega2*cons(IM3,k,j,i));
        cons_out(IM2,k,j,i) += 2.*dt*(omega1*cons(IM3,k,j,i) - omega3*cons(IM1,k,j,i));
        if (nx3 > 1) // 3D
          cons_out(IM3,k,j,i) += 2.*dt*(omega2*cons(IM1,k,j,i) - omega1*cons(IM2,k,j,i));
      }
}

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // read and save problem parameter
  Real jlat1  = pin->GetReal("problem", "jlat1")/180.*M_PI;
  Real jlat2  = pin->GetReal("problem", "jlat2")/180.*M_PI;
  Real umax   = pin->GetReal("problem", "umax");
  Real gh0    = pin->GetReal("problem", "gh0");
  int  vnum   = pin->GetInteger("problem", "vnum");
  Real vrad = pin->GetReal("problem", "vrad");
  Real vlat = pin->GetReal("problem", "vlat")/180.*M_PI;
  Real vgh  = pin->GetReal("problem", "vgh");
  int  ntracers = pin->GetInteger("problem", "ntracers");

  // coriolis parameter
  omegax = pin->GetOrAddReal("problem", "omegax", 0.);
  omegay = pin->GetOrAddReal("problem", "omegay", 0.);
  omegaz = pin->GetOrAddReal("problem", "omegaz", 0.);

  // planet radius
  Real radius = pcoord->x3v(0);

  // setup vortex longitude
  Real *vlon = new Real [vnum];
  for (int n = 0; n < vnum; ++n)
    vlon[n] = 2.*M_PI*n/vnum;

  // particle
  ppg = new ParticleGroup(this, "tracer", ParticleTranslate);

  Real x1min = block_size.x1min*0.99;
  Real x1max = block_size.x1max*0.99;
  Real x2min = block_size.x2min*0.99;
  Real x2max = block_size.x2max*0.99;
  Real lat, lon;
  for (int n2 = 0; n2 < sqrt(ntracers); ++n2)
    for (int n1 = 0; n1 < sqrt(ntracers); ++n1) {
      Particle tracer;
      tracer.x1 = x1min + (x1max - x1min)*(n1/sqrt(ntracers));
      tracer.x2 = asin(sin(x2min) + (sin(x2max) - sin(x2min))*(n2/sqrt(ntracers)));
      tracer.x3 = 0.;
      RotateEquatorToPole(&lat, &lon, tracer.x2, tracer.x1);
      tracer.time = lat;
      ppg->q.push_back(tracer);
  }

  // setup initial condition
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real theta, phi, xx, yy, zz, u1, u2, u3;
        RotateEquatorToPole(&theta, &phi, pcoord->x2v(j), pcoord->x1v(i));

        // setup mean flow with jet
        phydro->w(IDN,k,j,i) = gh0 + Geopotential(theta, jlat1, jlat2, umax, omegax, radius);
        SphericalLatlonToCartesian(&xx, &yy, &zz, ZonalWind(theta, jlat1, jlat2, umax), 0., 0.,
          phi, theta);
        CartesianToSphericalLatlon(&phydro->w(IVX,k,j,i), &phydro->w(IVY,k,j,i), &u3, zz, xx, yy,
          pcoord->x1v(i), pcoord->x2v(j));

        // add vortices around the pole
        // vortex center position is C = [vlat, 2*pi*n/vnum]
        // vortex height is gh_vortex * exp(-0.5*sqr(X - C)/sqr(vrad))
        for (int n = 0; n < vnum; ++n) {
          Real dist = radius*acos(sin(vlat)*sin(theta) + cos(vlat)*cos(theta)*cos(vlon[n] - phi));
          phydro->w(IDN,k,j,i) += vgh*exp(-0.5*_sqr(dist)/_sqr(vrad));
        }

        // another vortex at the pole
        Real dist = radius*(M_PI/2. - theta);
        phydro->w(IDN,k,j,i) += vgh*exp(-0.5*_sqr(dist)/_sqr(vrad));

        // add small pertubation
        //phydro->w(IDN,k,j,i) += 240.*cos(theta)*exp(-_sqr(10.*phi))
        //  *exp(-_sqr(30.*((jlat1 + jlat2)/2. - theta)));
      }

  // convert to conserved variables
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  delete[] vlon;
}
