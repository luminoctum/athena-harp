#ifndef ATHENA_HPP
#define ATHENA_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena.hpp
//  \brief contains Athena++ general purpose types, structures, enums, etc.

#include "defs.hpp"
#include "athena_arrays.hpp"
#include <math.h>

// typedefs that allow code to run with either floats or doubles
typedef double Real;
#ifdef MPI_PARALLEL
#define MPI_ATHENA_REAL MPI_DOUBLE
#endif

class MeshBlock;
class Coordinates;
struct RegionSize;
struct Particle;
struct Reaction;

//---------------------------------------------------------------------------------------
//! \struct FaceField
//  \brief container for face-centered fields

typedef struct FaceField {
  AthenaArray<Real> x1f,x2f,x3f;
} FaceField;

//----------------------------------------------------------------------------------------
//! \struct EdgeField
//  \brief container for edge-centered fields

typedef struct EdgeField {
  AthenaArray<Real> x1e,x2e,x3e;
} EdgeField;

//----------------------------------------------------------------------------------------
// enums used everywhere

// array indices for conserved: density, momemtum, total energy, face-centered field 
enum {IDN=0, IM1=NCOMP, IM2=NCOMP+1, IM3=NCOMP+2, IEN=NCOMP+3};
enum {IB1=0, IB2=1, IB3=2};

// array indices for 1D primitives: velocity, transverse components of field
enum {IT=0, IVX=NCOMP, IVY=NCOMP+1, IVZ=NCOMP+2, IPR=NCOMP+3, IBY=NHYDRO, IBZ=NHYDRO+1};

// array indices for face-centered electric fields returned by Riemann solver
enum {X1E2=0, X1E3=1, X2E3=0, X2E1=1, X3E1=0, X3E2=1};

// array indices for metric in GR
enum {I00, I01, I02, I03, I11, I12, I13, I22, I23, I33, NMETRIC};

// array indices for radiation model
enum {ITAU=1, ISA=2, IPM=3};

// needed for arrays dimensioned over grid directions
enum CoordinateDirection {X1DIR=0, X2DIR=1, X3DIR=2};

// needed wherever MPI communications are used.  Must be < 32 and unique
enum Athena_MPI_Tag {TAG_HYDRO=0, TAG_FIELD=1, TAG_RAD=2, TAG_CHEM=3, TAG_HYDFLX=4,
  TAG_FLDFLX=5, TAG_RADFLX=6, TAG_CHMFLX=7, TAG_AMR=8, TAG_FLDFLX_POLE=9, TAG_WTLIM=10,
  TAG_PARTICLE=11, TAG_PARTICLE_NUM=12};

// neeed to identify the phase of a molecule
enum PhaseID {GAS=0, LIQUID=1, SOLID=2};

//----------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime

typedef void (*BValFunc_t)(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
  FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
typedef int (*AMRFlagFunc_t)(MeshBlock *pmb);
typedef Real (*MeshGenFunc_t)(Real x, RegionSize rs);
typedef void (*SrcTermFunc_t)(MeshBlock *pmb, const Real time, const Real dt, const int step,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
typedef Real (*TimeStepFunc_t)(MeshBlock *pmb);
typedef Real (*HistoryOutputFunc_t)(MeshBlock *pmb, int iout);
typedef bool (*ParticleUpdateFunc_t)(MeshBlock *pmb, Real const time, Real const dt,
  Particle& pt, AthenaArray<Real> const& prim, AthenaArray<Real> const& cons, 
  int kji[3], AthenaArray<Real>& cons_out);
typedef Real (*ReactionFunc_t)(Reaction const& rc, Real const prim[], Real);

#endif // ATHENA_HPP
