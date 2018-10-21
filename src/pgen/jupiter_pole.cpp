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
#include "../utils/utils.hpp" // ran2
#include "../particle/particle.hpp"

// Coriolis parameter
Real omegax, omegay, omegaz;

// History output total absolute angular momentum and energy
Real TotalAbsoluteAngularMomentum(MeshBlock *pm, int iout);
Real TotalEnergy(MeshBlock *pm, int iout);

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

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, TotalAbsoluteAngularMomentum, "AM");
  EnrollUserHistoryOutput(1, TotalEnergy, "EN");
}

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // read and save problem parameter
  int  vnum = pin->GetInteger("problem", "vnum");
  Real vrad = pin->GetReal("problem", "vrad");
  Real vlat = pin->GetReal("problem", "vlat")/180.*M_PI;
  Real vgh  = pin->GetReal("problem", "vgh");
  Real gh0  = pin->GetReal("problem", "gh0");
  int  ntracers = pin->GetInteger("problem", "ntracers");

  // coriolis parameter
  omegax = pin->GetOrAddReal("problem", "omegax", 0.);
  omegay = pin->GetOrAddReal("problem", "omegay", 0.);
  omegaz = pin->GetOrAddReal("problem", "omegaz", 0.);

  Real radius = pcoord->x3v(0);

  // setup vortex longitude
  Real *vlon = new Real [vnum];
  for (int n = 0; n < vnum; ++n)
    vlon[n] = 2.*M_PI*n/vnum;

  // setup initial condition
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real theta, phi, xx, yy, zz, u1, u2, u3;
        RotateEquatorToPole(&theta, &phi, pcoord->x2v(j), pcoord->x1v(i));

        // add vortices around the pole
        // vortex center position is C = [vlat, 2*pi*n/vnum]
        // vortex height is gh_vortex * exp(-0.5*sqr(X - C)/sqr(vrad))
        phydro->u(IDN,k,j,i) = gh0;
        for (int n = 0; n < vnum; ++n) {
          Real dist = radius * acos(sin(vlat)*sin(theta) + cos(vlat)*cos(theta)*cos(vlon[n]-phi));
          phydro->u(IDN,k,j,i) += vgh * exp(-0.5*_sqr(dist)/_sqr(vrad));
        }

        // another vortex at the pole
        Real dist = radius * (M_PI/2. - theta);
        phydro->u(IDN,k,j,i) += vgh * exp(-0.5*_sqr(dist)/_sqr(vrad));
      }

  // randomly distribute tracers over the domain
  ppg = new ParticleGroup(this, "tracer", ParticleTranslate);
  long int iseed = -1 - Globals::my_rank;

  Real x1min = block_size.x1min*0.99;
  Real x1max = block_size.x1max*0.99;
  Real x2min = block_size.x2min*0.99;
  Real x2max = block_size.x2max*0.99;

  /*for (int n = 0; n < ntracers; ++n) {
    Particle tracer;
    tracer.time = 0.;
    tracer.x3 = 0.;
    tracer.x2 = asin(sin(x2min) + (sin(x2max) - sin(x2min))*ran2(&iseed));
    tracer.x1 = x1min + (x1max - x1min) * ran2(&iseed);
    ppg->q.push_back(tracer);
  }*/

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

  delete[] vlon;
}


Real TotalAbsoluteAngularMomentum(MeshBlock *pmb, int iout)
{
  AthenaArray<Real> vol;
  int ncells1 = pmb->block_size.nx1 + 2*NGHOST;
  vol.NewAthenaArray(ncells1);
  Hydro *phyd = pmb->phydro;
  Real am = 0.;
  Real omega = phyd->psrc->GetOmegaZ();

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real phi = phyd->u(IDN,k,j,i);
        Real r = pmb->pcoord->x3v(k);
        Real theta = pmb->pcoord->x2v(j);
        Real u = phyd->u(IM1,k,j,i);
        am += vol(i) * (omega*phi*r*cos(theta) + u) * r*cos(theta);
      }
    }

  return am;
}

Real TotalEnergy(MeshBlock *pmb, int iout)
{
  AthenaArray<Real> vol;
  int ncells1 = pmb->block_size.nx1 + 2*NGHOST;
  vol.NewAthenaArray(ncells1);
  Hydro *phyd = pmb->phydro;
  Real en = 0.;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real phi = phyd->u(IDN,k,j,i);
        Real ke1 = 0.5 * _sqr(phyd->u(IM1,k,j,i)) / phi;
        Real ke2 = 0.5 * _sqr(phyd->u(IM2,k,j,i)) / phi;
        Real pe = 0.5 * _sqr(phi);
        en += vol(i) * (ke1 + ke2 + pe);
      }
    }

  return en;
}
