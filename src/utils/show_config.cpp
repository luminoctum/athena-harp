//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file show_config.cpp 

// C++ headers
#include <iostream>
#include <sstream>

// Athena headers
#include "../athena.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ShowConfig(void)
//  \brief prints diagnostic messages about the configuration of an Athena++ executable

void ShowConfig(void)
{
  std::cout<<"This Athena++ executable is configured with:" << std::endl;
  std::cout<<"  Problem generator:          " << PROBLEM_GENERATOR << std::endl;
  std::cout<<"  Coordinate system:          " << COORDINATE_SYSTEM << std::endl;
  std::cout<<"  Equation of state:          " << EQUATION_OF_STATE << std::endl;
  std::cout<<"  Riemann solver:             " << RIEMANN_SOLVER << std::endl;
  std::cout<<"  Reconstruction method:      " << RECONSTRUCTION_METHOD << std::endl;
  std::cout<<"  Hydro integrator:           " << HYDRO_TIME_INTEGRATOR << std::endl;
  if (MAGNETIC_FIELDS_ENABLED) {
    std::cout<<"  Magnetic fields:            ON" << std::endl;
  } else {
    std::cout<<"  Magnetic fields:            OFF" << std::endl;
  }
  if (RELATIVISTIC_DYNAMICS) {
    std::cout<<"  Relativistic dynamics:      ON" << std::endl;
  } else {
    std::cout<<"  Relativistic dynamics:      OFF" << std::endl;
  }
  if (GENERAL_RELATIVITY) {
    std::cout<<"  General relativity:         ON" << std::endl;
    if (FRAME_TRANSFORMATIONS) {
      std::cout<<"  Frame transformations:      ON" << std::endl;
    } else {
      std::cout<<"  Frame transformations:      OFF" << std::endl;
    }
  } else {
    std::cout<<"  General Relativity:         OFF" << std::endl;
  }
#ifdef MPI_PARALLEL
  std::cout<<"  MPI parallelism:            ON" << std::endl;
#else
  std::cout<<"  MPI parallelism:            OFF" << std::endl;
#endif
#ifdef OPENMP_PARALLEL
  std::cout<<"  OpenMP parallelism:         ON" << std::endl;
#else
  std::cout<<"  OpenMP parallelism:         OFF" << std::endl;
#endif
#ifdef HDF5OUTPUT
  std::cout<<"  HDF5 Output:                ON" << std::endl;
#else
  std::cout<<"  HDF5 Output:                OFF" << std::endl;
#endif
  std::cout<<"  Compiler:                   " << COMPILED_WITH << std::endl;
  std::cout<<"  Compilation command:        " << COMPILER_COMMAND << ' '
      << COMPILED_WITH_OPTIONS << std::endl;
  return;
}
