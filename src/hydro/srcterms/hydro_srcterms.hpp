#ifndef HYDRO_SRCTERMS_HPP
#define HYDRO_SRCTERMS_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_srcterms.hpp
//  \brief defines class HydroSourceTerms
//  Contains data and functions that implement physical (not coordinate) source terms

// Athena headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

// Forward declarations
class Hydro;
class ParameterInput;

//! \class HydroSourceTerms
//  \brief data and functions for physical source terms in the hydro

class HydroSourceTerms {
public:
  HydroSourceTerms(Hydro *phyd, ParameterInput *pin);
  ~HydroSourceTerms();

  // accessors
  Real GetGM() const {return gm_;}
  Real GetG1() const {return g1_;}
  Real GetG2() const {return g2_;}
  Real GetG3() const {return g3_;}
  Real GetOmega1() const {return omega1_;}
  Real GetOmega2() const {return omega2_;}
  Real GetOmega3() const {return omega3_;}
  Real GetOmegaX() const {return omegax_;}
  Real GetOmegaY() const {return omegay_;}
  Real GetOmegaZ() const {return omegaz_;}

  // data
  bool hydro_sourceterms_defined;

  // functions
  void AddHydroSourceTerms(const Real time, const Real dt, const int step, const AthenaArray<Real> *flx,
    const AthenaArray<Real> &p, const AthenaArray<Real> &b, AthenaArray<Real> &c);
  void PointMass(const Real dt, const AthenaArray<Real> *flx,const AthenaArray<Real> &p,
    AthenaArray<Real> &c);
  void ConstantAcceleration(const Real dt, const AthenaArray<Real> *flx,
    const AthenaArray<Real> &p, AthenaArray<Real> &c);
  void Coriolis123(const Real dt, const AthenaArray<Real> *flx,
    const AthenaArray<Real> &p, AthenaArray<Real> &c);
  void CoriolisXYZ(const Real dt, const AthenaArray<Real> *flx,
    const AthenaArray<Real> &p, AthenaArray<Real> &c);
  void EnrollSrcTermFunction(SrcTermFunc_t my_func);
  SrcTermFunc_t UserSourceTerm;

private:
  Hydro *pmy_hydro_;  // ptr to Hydro containing this HydroSourceTerms

  Real gm_;           // GM for point mass MUST BE LOCATED AT ORIGIN
  Real g1_, g2_, g3_; // constant acc'n in each direction

  Real omega1_, omega2_, omega3_; // constant coriolis acc'n in each direction
  Real omegax_, omegay_, omegaz_; // coriolis acc'n in cartesian direction
};
#endif
