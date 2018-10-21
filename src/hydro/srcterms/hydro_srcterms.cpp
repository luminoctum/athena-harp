//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief Class to implement source terms in the hydro equations

// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "hydro_srcterms.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../parameter_input.hpp"

// HydroSourceTerms constructor 

HydroSourceTerms::HydroSourceTerms(Hydro *phyd, ParameterInput *pin)
{
  pmy_hydro_ = phyd;
  hydro_sourceterms_defined = false;

  // read point mass or constant acceleration parameters from input block

  // set the point source only when the coordinate is spherical or 2D cylindrical.
  gm_ = pin->GetOrAddReal("problem","GM",0.0);
  if (gm_ != 0.0) {
    if (COORDINATE_SYSTEM == "spherical_polar"
    || (COORDINATE_SYSTEM == "cylindrical" && phyd->pmy_block->block_size.nx3==1)) {
      hydro_sourceterms_defined = true;
    }
    else {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms constructor" << std::endl
          << "The point source gravity works only in spherical polar coordinates " 
          << "or in 2D cylindrical coordinates." << std::endl
          << "Check <problem> GM parameter in the input file." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
  g1_ = pin->GetOrAddReal("hydro","grav_acc1",0.0);
  if (g1_ != 0.0) hydro_sourceterms_defined = true;

  g2_ = pin->GetOrAddReal("hydro","grav_acc2",0.0);
  if (g2_ != 0.0) hydro_sourceterms_defined = true;

  g3_ = pin->GetOrAddReal("hydro","grav_acc3",0.0);
  if (g3_ != 0.0) hydro_sourceterms_defined = true;

  omega1_ = pin->GetOrAddReal("hydro","coriolis_acc1",0.0);
  if (omega1_ != 0.0) hydro_sourceterms_defined = true;

  omega2_ = pin->GetOrAddReal("hydro","coriolis_acc2",0.0);
  if (omega2_ != 0.0) hydro_sourceterms_defined = true;

  omega3_ = pin->GetOrAddReal("hydro","coriolis_acc3",0.0);
  if (omega3_ != 0.0) hydro_sourceterms_defined = true;

  omegax_ = pin->GetOrAddReal("hydro","coriolis_accx",0.0);
  if (omegax_ != 0.0) hydro_sourceterms_defined = true;

  omegay_ = pin->GetOrAddReal("hydro","coriolis_accy",0.0);
  if (omegay_ != 0.0) hydro_sourceterms_defined = true;

  omegaz_ = pin->GetOrAddReal("hydro","coriolis_accz",0.0);
  if (omegaz_ != 0.0) hydro_sourceterms_defined = true;

  UserSourceTerm = phyd->pmy_block->pmy_mesh->UserSourceTerm_;
  if(UserSourceTerm != NULL) hydro_sourceterms_defined = true;
}

// destructor

HydroSourceTerms::~HydroSourceTerms()
{
}

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::AddHydroSourceTerms
//  \brief Adds source terms to conserved variables

void HydroSourceTerms::AddHydroSourceTerms(const Real time, const Real dt, const int step,
     const AthenaArray<Real> *flux, const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  // accleration due to point mass (MUST BE AT ORIGIN)
  if (gm_ != 0.0) PointMass(dt, flux, prim, cons);

  // constant acceleration (e.g. for RT instability)
  if (g1_ != 0.0 || g2_ != 0.0 || g3_ != 0.0)
    ConstantAcceleration(dt, flux, prim,cons);

  // coriolis acceleration in the axial direction
  if (omega1_ != 0.0 || omega2_ != 0.0 || omega3_ != 0.0)
    Coriolis123(dt, flux, prim, cons);

  // coriolis acceleration in the cartesian direction
  if (omegax_ != 0.0 || omegay_ != 0.0 || omegaz_ != 0.0)
    CoriolisXYZ(dt, flux, prim, cons);

  //  user-defined source terms
  if (UserSourceTerm != NULL)
    UserSourceTerm(pmb,time,dt,step,prim,bcc,cons);

  // Add new source terms here
  // MyNewSourceTerms()

  return;
}

