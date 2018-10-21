//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file change_rundir.cpp 
//! \brief executes unix 'chdir' command to change dir in which Athena++ runs

// Athena headers
#include "../athena.hpp"

// C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>  // mkdir()
#include <unistd.h>    // chdir()

//----------------------------------------------------------------------------------------
//! \fn void ChangeRunDir(const char *pdir)
//  \brief change to input run directory; create if it does not exist yet

void ChangeRunDir(const char *pdir)
{
  std::stringstream msg;

  if (pdir == NULL || *pdir == '\0') return;

  mkdir(pdir, 0775);
  if(chdir(pdir)) {
    msg << "### FATAL ERROR in function [ChangeToRunDir]" << std::endl
        << "Cannot cd to directory '" << pdir << "'";
    throw std::runtime_error(msg.str().c_str());
  }

  return;
}
