//! \file sw.cpp
//  \brief Showman 2006 shallow water model

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
static Real tau_interval;
static Real tau_storm;
static Real smax;
static Real tau_mass;
static Real tau_ape;
static Real gheq;

static Real f0;
static Real beta;

// mass pulse
struct MassPulse {
  Real tpeak;
  Real x1;
  Real x2;
};

// all mass pulses
static std::vector<MassPulse> storms;

// last time of mass pulse in this MeshBlock
static Real t_last_storm = -1.E-8;

// seed of random number
static long int iseed = -1;

// external forcing
void MassPulseRadiationSink(MeshBlock *pmb, const Real time, const Real dt, const int step,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

// support functions
Real get_zonal_wind(Real theta);
Real get_geopotential(Real phi0, Real theta);
Real solve_phi0(Real phi0);

// History output total energy
Real TotalEnergy(MeshBlock *pm, int iout);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, TotalEnergy, "EN");
  EnrollUserExplicitSourceFunction(MassPulseRadiationSink);
}

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  r_storm       = pin->GetReal("problem", "r_storm");
  tau_interval  = pin->GetReal("problem", "tau_interval");
  tau_storm     = pin->GetReal("problem", "tau_storm");
  smax          = pin->GetReal("problem", "smax");
  tau_mass      = pin->GetReal("problem", "tau_mass");
  tau_ape       = pin->GetReal("problem", "tau_ape");
  gheq          = pin->GetReal("problem", "gheq");
  f0            = pin->GetReal("problem", "f0");
  beta          = pin->GetReal("problem", "beta");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        // setup initial thickness
        phydro->u(IDN,k,j,i) = gheq;
      }
}

Real get_zonal_wind(Real theta)
{
}

Real get_geopotential(Real phi0, Real theta)
{
}

Real solve_phi0(Real phi)
{
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

void MassPulseRadiationSink(MeshBlock *pmb, const Real time, const Real dt, const int step,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  int current_storms = storms.size();
  Real x1min = pmb->pmy_mesh->mesh_size.x1min;
  Real x1max = pmb->pmy_mesh->mesh_size.x1max;
  Real x2min = pmb->pmy_mesh->mesh_size.x2min;
  Real x2max = pmb->pmy_mesh->mesh_size.x2max;

  // insert new storm
  if ((step == 1) && (ran2(&iseed) < exp(-tau_interval/(time - t_last_storm)))) {
    MassPulse storm;
    storm.tpeak = time + tau_storm/2.;
    storm.x1 = x1min + (x1max - x1min) * ran2(&iseed);
    storm.x2 = x2min + (x2max - x2min) * ran2(&iseed);
    storms.push_back(storm);
    t_last_storm = time;
  }

  // calculate average height for this instant
  Real gh_avg[2];
  AthenaArray<Real> vol;
  int ncells1 = pmb->block_size.nx1 + 2*NGHOST;
  vol.NewAthenaArray(ncells1);

  gh_avg[0] = 0.;
  gh_avg[1] = 0.;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        gh_avg[0] += prim(IDN,k,j,i) * vol(i);
        gh_avg[1] += vol(i);
      }
    }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, gh_avg, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  gh_avg[0] /= gh_avg[1];

  //if (Globals::my_rank == 0)
  //  std::cout << gh_avg[0] << std::endl;

  // add storm source, radiation, and beta plane
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real x1a = pmb->pcoord->x1v(i);
        Real x2a = pmb->pcoord->x2v(j);
        for (size_t n = 0; n < current_storms; ++n) {
          Real x1b = storms[n].x1;
          Real x2b = storms[n].x2;
          Real dr = sqrt(_sqr(x1a - x1b) + _sqr(x2a - x2b));
          if (dr > 2.2 * r_storm) continue;

          // mass pulse
          cons(IDN,k,j,i) += dt * smax * exp(-_sqr(dr)/_sqr(r_storm) - _sqr(time - storms[n].tpeak)/_sqr(tau_storm));
        }
        // radiation
        cons(IDN,k,j,i) += - dt * (gh_avg[0] - gheq)/tau_mass;

        // beta plane
        cons(IM1,k,j,i) += dt * (f0 + beta * x2a) * prim(IDN,k,j,i) * prim(IM2,k,j,i);
        cons(IM2,k,j,i) += - dt * (f0 + beta * x2a) * prim(IDN,k,j,i) * prim(IM1,k,j,i);
      }

  // remove deceased storms
  if ((step == 1) && (storms.size() > 0)) {
    size_t n = 0;
    for (; n < storms.size(); ++n)
      if (time - storms[n].tpeak < 2.2*tau_storm)
        break;
    storms.erase(storms.begin(), storms.begin() + n);
  }
}
