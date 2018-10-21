#ifndef GLOBALS_HPP
#define GLOBALS_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file globals.hpp
//  \brief namespace containing external global variables

namespace Globals
{
  extern int my_rank, nranks;

  extern double const Rgas;
  extern double const kBoltz;
  extern double const Lo;
}

#endif // GLOBALS_HPP
