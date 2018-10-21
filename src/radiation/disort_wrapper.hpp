#ifndef DISORT_WRAPPER_HPP
#define DISORT_WRAPPER_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "../athena.hpp"


class ParameterInput;
template<typename T> class AthenaArray;

#ifdef USE_DISORT
// third party software DISORT headers
extern "C" {
  #include "cdisort.h"
}

enum AtmDir {TOP2BOT = true, BOT2TOP = false};

class DisortWrapper {
public:
  DisortWrapper(ParameterInput *pin, int nwave, AtmDir dir_);
  ~DisortWrapper();

  // data
  std::vector<disort_bc> bc;

  // functions
  void Run(AthenaArray<Real>& rad, AthenaArray<Real>& uu,
    AthenaArray<Real> const& oppr, Real *wave);

private:
  disort_state  ds;
  disort_output ds_out;
  AtmDir dir;
};
#endif // USE_DISORT

#endif
