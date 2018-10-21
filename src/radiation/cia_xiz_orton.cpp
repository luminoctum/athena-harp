// C/C++ headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "absorber.hpp"
#include "../math_funcs.hpp"
#include "../globals.hpp"
#include "../misc.hpp"

void CIAXiz::LoadCoefficient(std::string fname)
{
  std::stringstream msg;
  len_[1] = GetNumCols(fname) - 1;
  len_[0] = GetNumRows(fname) - 1;
  
  std::ifstream infile(fname.c_str(), std::ios::in);
  axis_.resize(len_[0] + len_[1]);
  kcoeff_.resize(len_[0]*len_[1]);
  Real junk;
  if (infile.is_open()) {
    infile >> junk;
    for (int j = 0; j < len_[1]; j++) {
      infile >> axis_[len_[0] + j];
    }
    for (int k = 0; k < len_[0]; k++) {
      infile >> axis_[k];
      for (int j = 0; j < len_[1]; j++)
        infile >> kcoeff_[k*len_[1] + j];
    }
    infile.close();
  } else {
    msg << "### FATAL ERROR in CIAXiz::LoadCoefficient: ";
    msg << "Cannot open file: " << fname << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
}

Real CIAXiz::AbsorptionCoefficient(Real wave, Real const prim[]) const
{
  // first axis is wavenumber, second is temperature
  Real val, coord[2] = {wave, prim[IT]};
  _interpn(&val, coord, kcoeff_.data(), axis_.data(), len_, 2);

  Real amagat = prim[IPR]/(Globals::kBoltz*prim[IT]*Globals::Lo);
  Real xt = 0.;
  for (int i = 1; i < NCOMP; ++i) xt += prim[i];
  Real x1 = id_[0] == 0 ? 1. - xt : prim[id_[0]];
  Real x2 = id_[1] == 0 ? 1. - xt : prim[id_[1]];

  return 100.*exp(-val)*x1*x2*amagat*amagat; // 1/cm -> 1/m
}
