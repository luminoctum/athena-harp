#ifndef RADIATION_HPP
#define RADIATION_HPP

// C/C++ headers
#include <string>

// Athena++ headers
#include "../athena.hpp"

class MeshBlock;
class ParameterInput;
class Absorber;
template<typename T> class AthenaArray;

class Radiation {
public:
  // data
  MeshBlock *pmy_block;
  std::string myname;
  Radiation *prev, *next;
  Real *wave, *level;
  int nlevel, nwave, nmom, nphi, ntau, numu;
  AthenaArray<Real> oppr, rad, uu;

  // functions
  Radiation(MeshBlock *pmb, ParameterInput *pin, std::string name);
  ~Radiation();
  Radiation* Next(int n);
  Radiation* AddRadiation(MeshBlock *pmb, ParameterInput *pin, std::string name);
  void InitAbsorbers(ParameterInput *pin, int type = 0);
  void SetOpticalProperties(AthenaArray<Real> const& prim,
      int k, int j, int is, int ie);
  void TotalFlux(AthenaArray<Real>& flux) const;
  void WriteTopFlux(std::string fname) const;
  void WriteTopRadiance(std::string fname) const;
  void WriteOpticalDepth(std::string fname) const;

protected:
  Absorber *pabs;
  Real *weight_;
};

void WriteHeatingRate(std::string fname, AthenaArray<Real> const& flux,
      AthenaArray<Real> const& hr, Real const* level);

#endif
