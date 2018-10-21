//! \file sw.cpp
//  \brief Galewsky 2004 test of global shallow water model

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../athena_math.hpp"     // _root, _sqr
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"

// model parameters shared by all subroutine
static Real theta0;
static Real theta1;
static Real umax;
static Real phim;
static Real omega;
static Real radius;

// support functions
Real get_zonal_wind(Real theta);
Real get_geopotential(Real phi0, Real theta);
Real solve_phi0(Real phi0);

// History output total absolute angular momentum and energy
Real TotalAbsoluteAngularMomentum(MeshBlock *pm, int iout);
Real TotalEnergy(MeshBlock *pm, int iout);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, TotalAbsoluteAngularMomentum, "AM");
  EnrollUserHistoryOutput(1, TotalEnergy, "EN");
}

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  theta0 = pin->GetOrAddReal("problem", "theta0", M_PI/7.);
  theta1 = pin->GetOrAddReal("problem", "theta1", M_PI/2. - theta0);
  umax = pin->GetOrAddReal("problem", "umax", 80.);
  phim = pin->GetOrAddReal("problem", "phim", 98061.6);
  omega = phydro->psrc->GetOmegaZ();
  radius = pcoord->x3v(0);

  // sovle for phi0
  Real phi0;
  int err = _root(phim/2., phim*2., 1., &phi0, solve_phi0);

  if (err != 0)
    throw std::runtime_error("FATAL ERROR: Galewsky04 solve phi0 does not converge");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real lambda = pcoord->x1v(i);
        Real theta = pcoord->x2v(j);
        // setup mean flow
        phydro->u(IDN,k,j,i) = get_geopotential(phi0, theta);
        phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i) * get_zonal_wind(theta);

        // add pertubation
        phydro->u(IDN,k,j,i) += 1176.74*cos(theta)*exp(-_sqr(3.*(lambda-M_PI)))
          * exp(-_sqr(15.*(M_PI/4.-theta)));
      }
}

Real get_zonal_wind(Real theta)
{
  if (theta <= theta0 || theta >= theta1)
    return 0.;

  Real en = exp(-4./_sqr(theta1 - theta0));

  return umax/en * exp(1./((theta-theta0)*(theta-theta1)));
}

Real get_geopotential(Real phi0, Real theta)
{
  Real ds = 0.01;
  for (Real s = - M_PI/2. + ds; s < theta; s += ds) {
    Real f1r = 2.*omega*sin(s-ds)*radius,
         f2r = 2.*omega*sin(s)*radius;

    Real u1 = get_zonal_wind(s-ds),
         u2 = get_zonal_wind(s);

    phi0 -= 0.5 * ((f1r*u1 + u1*u1*tan(s-ds)) + (f2r*u2 + u2*u2*tan(s))) * ds;
  }

  return phi0;
}

Real solve_phi0(Real phi)
{
  Real total_phi = 0.;
  Real ds = 0.01, theta = - M_PI/2. + ds;

  for (; theta < M_PI/2.; theta += ds) {
    Real f1r = 2.*omega*sin(theta-ds)*radius,
         f2r = 2.*omega*sin(theta)*radius;

    Real u1 = get_zonal_wind(theta-ds),
         u2 = get_zonal_wind(theta);

    phi -= 0.5 * ((f1r*u1 + u1*u1*tan(theta-ds)) + (f2r*u2 + u2*u2*tan(theta))) * ds;
    total_phi += phi * cos(theta) * ds;
  }

  return total_phi/2. - phim;
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
