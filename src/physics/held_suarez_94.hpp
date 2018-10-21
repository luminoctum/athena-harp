#ifndef HELD_SUAREZ_94_HPP
#define HELD_SUAREZ_94_HPP

#include "../athena.hpp"  // Real, MeshBlock

//! \file held_suarez_94.hpp
//  \brief classes and functions related to Held & Suarez, 94

// Forward declarations
class ParameterInput;

//! \class HeldSuarez94
//  \brief contains all Held & Suarez 94 parameters
class HeldSuarez94
{
  //friend class MeshBlock;
public:
  //! default constructor
  HeldSuarez94();

  //! read parameters from input file
  void LoadInputFile(ParameterInput *pin);

  //! \return equilibrium temperature at pressure and latitude
  Real GetTempEq(Real theta, Real pres) const;

  //! used in solving for hydrostatic equilibrium
  Real operator()(Real ptop) const;

  //! \return Rayleigh friction coefficient
  Real Kv(Real sigma) const;

  //! \return Newtonian cooling coefficient
  Real Kt(Real theta, Real sigma) const;

  // acess functions
  void SetLatitude(Real lat_) { lat = lat_; }
  void SetBottomPressure(Real pbot_) { pbot = pbot_; }
  void SetDistance(Real dz_) { dz = dz_; }
  Real GetSurfacePressure() const { return psrf; }
  Real GetRgas() const { return rgas; }

protected:
  Real tdy;
  Real tdz;
  Real psrf;
  Real tsrf;
  Real tmin;

  Real kappa;
  Real rgas;
  Real grav;

  Real sigma_b;
  Real kf;
  Real ka;
  Real ks;

  Real lat;   //!< latitude
  Real pbot;  //!< pressure of the bottom cell
  Real dz;    //!< distance between two cells
};


#endif
