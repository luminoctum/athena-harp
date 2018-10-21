//! \file swgv.cpp
//  \brief motion of a single vortex in a shallow water model

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <vector>

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
static Real r_storm;
static Real smax;
static Real gheq;
static Real x1;
static Real x2;

static Real f0;
static Real beta;

// History output total energy
Real ShallowWaterTotalEnergy(MeshBlock *pm, int iout);

// Forcing
void BetaPlaneForcing(MeshBlock *pmb, Real const time, Real const dt,
    int const step, AthenaArray<Real> const& prim, AthenaArray<Real> const& bcc,
    AthenaArray<Real>& cons);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, ShallowWaterTotalEnergy, "EN");
  EnrollUserExplicitSourceFunction(BetaPlaneForcing);
}

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  r_storm = pin->GetReal("problem", "r_storm");
  smax    = pin->GetReal("problem", "smax");
  gheq    = pin->GetReal("problem", "gheq");

  x1      = pin->GetReal("problem", "x1");
  x2      = pin->GetReal("problem", "x2");    

  f0      = pin->GetReal("problem", "f0");
  beta    = pin->GetReal("problem", "beta");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        // setup initial thickness and heigh anomaly
        Real x1p = pcoord->x1v(i) - x1;
        Real x2p = pcoord->x2v(j) - x2;
        phydro->u(IDN,k,j,i) = gheq + smax * exp(-(_sqr(x1p)+_sqr(x2p))/(2.*_sqr(r_storm)));
      }
}

void BetaPlaneForcing(MeshBlock *pmb, Real const time, Real const dt,
    int const step, AthenaArray<Real> const& prim, AthenaArray<Real> const& bcc,
    AthenaArray<Real>& cons)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      Real x2 = pmb->pcoord->x2v(j);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        cons(IM1,k,j,i) += dt * (f0 + beta * x2) * prim(IDN,k,j,i) * prim(IM2,k,j,i);
        cons(IM2,k,j,i) += - dt * (f0 + beta * x2) * prim(IDN,k,j,i) * prim(IM1,k,j,i);
      }
    }
}

Real ShallowWaterTotalEnergy(MeshBlock *pmb, int iout)
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
