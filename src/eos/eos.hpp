#ifndef EOS_HPP
#define EOS_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos.hpp
//  \brief defines class EquationOfState
//  Contains data and functions that implement the equation of state

// Athena headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../coordinates/coordinates.hpp" // Coordinates
#include "../globals.hpp"					// Rgas

// Declarations
class Hydro;
class HydroSourceTerms;
class ParameterInput;
struct FaceField;

//! \class EquationOfState
//  \brief data and functions that implement EoS

class EquationOfState {
  friend class Hydro;
  friend class HydroSourceTerms;
  friend class ReactionGroup;
public:
  EquationOfState(MeshBlock *pmb) :
    pmy_block_(pmb)
  {}
  virtual ~EquationOfState() {}

  virtual void ConservedToPrimitive(AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old,
    const FaceField &b, AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int is, int ie, int js, int je, int ks, int ke) = 0;
  virtual void PrimitiveToConserved(const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
    AthenaArray<Real> &cons, Coordinates *pco,
    int is, int ie, int js, int je, int ks, int ke) = 0;

  virtual Real SoundSpeed(Real const prim[]) const {return 0.;}
  virtual Real Cp(Real const prim[]) const {return 0.;}
  virtual Real Cv(Real const prim[]) const {return 0.;}
  virtual Real Mass(Real const prim[]) const {return 0.;}
  virtual Real Entropy(Real const prim[]) const {return 0.;}
  virtual Real Enthalpy(Real const prim[]) const {return 0.;}
  Real Energy(Real const prim[]) const {
		return Enthalpy(prim) - Globals::Rgas;
	}

  Real GetDensityFloor() const {return density_floor_;}
  Real GetPressureFloor() const {return pressure_floor_;}
  Real GetGamma() const {return gamma_;}
  Real GetIsoSoundSpeed() const {return iso_sound_speed_;}
  Real GetCv(int n) const {return cv_[n];}
  //Real GetMu(int n) const {return mu_[n];}

protected:
  MeshBlock *pmy_block_;                        // ptr to MeshBlock containing this EOS
  Real density_floor_, pressure_floor_;         // density and pressure floors
  Real gamma_;                                  // ratio of specific heats
  Real iso_sound_speed_;                        // isothermal sound speed
  Real mu_[NCOMP], cv_[NCOMP], latent_[NCOMP];  // molecular weight, molar heat capacity and molar latent heat
};

class ShallowWaterHydro: public EquationOfState {
public:
  ShallowWaterHydro(MeshBlock *pmb, ParameterInput *pin);
  virtual ~ShallowWaterHydro();

  void ConservedToPrimitive(AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old,
    const FaceField &b, AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int is, int ie, int js, int je, int ks, int ke);
  void PrimitiveToConserved(const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
    AthenaArray<Real> &cons, Coordinates *pco,
    int is, int ie, int js, int je, int ks, int ke);

  Real SoundSpeed(Real const prim[]) const {
    return sqrt(prim[IDN]);
  }
};

class AdiabaticHydro: public EquationOfState {
public:
  AdiabaticHydro(MeshBlock *pmb, ParameterInput *pin);
  virtual ~AdiabaticHydro();

  void ConservedToPrimitive(AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old,
    const FaceField &b, AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int is, int ie, int js, int je, int ks, int ke);
  void PrimitiveToConserved(const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
    AthenaArray<Real> &cons, Coordinates *pco,
    int is, int ie, int js, int je, int ks, int ke);

  Real SoundSpeed(Real const prim[]) const {
    return gamma_*prim[IPR]/prim[IDN];
  }
};

class HeterogeneousHydro: public EquationOfState {
public:
  HeterogeneousHydro(MeshBlock *pmb, ParameterInput *pin);
  virtual ~HeterogeneousHydro();

  void ConservedToPrimitive(AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old,
    const FaceField &b, AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int is, int ie, int js, int je, int ks, int ke);
  void PrimitiveToConserved(const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
    AthenaArray<Real> &cons, Coordinates *pco,
    int is, int ie, int js, int je, int ks, int ke);

  Real SoundSpeed(Real const prim[]) const;
  Real Cp(Real const prim[]) const;
  Real Cv(Real const prim[]) const;
  Real Mass(Real const prim[]) const;
  Real Entropy(Real const prim[]) const;
  //Real Energy(Real const prim[]) const;
  Real Enthalpy(Real const prim[]) const;
};

class IsothermalHydro : public EquationOfState {
public:
  IsothermalHydro(MeshBlock *pmb, ParameterInput *pin);
  virtual ~IsothermalHydro();

  void ConservedToPrimitive(AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old,
    const FaceField &b, AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int is, int ie, int js, int je, int ks, int ke);
  void PrimitiveToConserved(const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
    AthenaArray<Real> &cons, Coordinates *pco,
    int is, int ie, int js, int je, int ks, int ke);

  Real SoundSpeed(Real const prim[]) {
    return iso_sound_speed_;
  }
};

class MagnetoAdiabaticHydro: public AdiabaticHydro {
public:
  MagnetoAdiabaticHydro(MeshBlock *pmb, ParameterInput *pin);
  Real FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx);
};

class MagnetoIsothermalHydro: public IsothermalHydro {
public:
  MagnetoIsothermalHydro(MeshBlock *pmb, ParameterInput *pin);
  Real FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx);
};

class SrHydro: public EquationOfState {
public:
  SrHydro(MeshBlock *pmb, ParameterInput *pin);
  void SoundSpeeds(Real rho_h, Real pgas, Real vx, Real gamma_lorentz_sq,
      Real *plambda_plus, Real *plambda_minus);
protected:
  AthenaArray<Real> g_, g_inv_;          // metric and its inverse, used in GR
  Real sigma_max_, beta_min_;            // limits on ratios of gas quantities to pmag
  Real gamma_max_;                       // maximum Lorentz factor
};

class SrMagnetoHydro: public SrHydro {
public:
  SrMagnetoHydro(MeshBlock *pmb, ParameterInput *pin);
  void FastMagnetosonicSpeeds(const AthenaArray<Real> &prim,
      const AthenaArray<Real> &bbx_vals, int il, int iu, int ivx,
      AthenaArray<Real> &lambdas_p, AthenaArray<Real> &lambdas_m);
};

class GrHydro: public EquationOfState {
public:
  GrHydro(MeshBlock *pmb, ParameterInput *pin);
  void SoundSpeeds(Real rho_h, Real pgas, Real u0, Real u1,
      Real g00, Real g01, Real g11,
      Real *plambda_plus, Real *plambda_minus);
protected:
  AthenaArray<Real> g_, g_inv_;          // metric and its inverse, used in GR
  Real sigma_max_, beta_min_;            // limits on ratios of gas quantities to pmag
  Real gamma_max_;                       // maximum Lorentz factor
  Real rho_min_, rho_pow_;               // variables to control power-law denity floor
  Real pgas_min_, pgas_pow_;             // variables to control power-law pressure floor
  AthenaArray<Real> fixed_;              // cells with problems, used in GR hydro
};

class GrMagnetoHydro: public GrHydro {
public:
  GrMagnetoHydro(MeshBlock *pmb, ParameterInput *pin);
  void FastMagnetosonicSpeeds(Real rho_h, Real pgas, Real u0, Real u1, Real b_sq,
      Real g00, Real g01, Real g11,
      Real *plambda_plus, Real *plambda_minus);
protected:
  Real sigma_max_, beta_min_;            // limits on ratios of gas quantities to pmag
  AthenaArray<Real> normal_dd_;          // normal-frame densities, used in GR MHD
  AthenaArray<Real> normal_ee_;          // normal-frame energies, used in GR MHD
  AthenaArray<Real> normal_mm_;          // normal-frame momenta, used in GR MHD
  AthenaArray<Real> normal_bb_;          // normal-frame fields, used in GR MHD
  AthenaArray<Real> normal_tt_;          // normal-frame M.B, used in GR MHD
};

#endif
