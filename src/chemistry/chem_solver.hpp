#ifndef CHEM_SOLVER_HPP_
#define CHEM_SOLVER_HPP_

// C++ headers
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../eos/eos.hpp"
#include "reaction.hpp"

struct IsentropicSolver {
  EquationOfState *peos;
  ReactionGroup *prg;
  Real entropy;
  Real *prim;
  int max_iter;
  char c;

  IsentropicSolver() : max_iter(10), c('T') {}
  Real operator()(Real v) {
    if (c == 'T') prim[IT] = v;
    else prim[IPR] = v;
    int err = prg->EquilibrateTP(prim, max_iter);
    if (err) {
      std::stringstream msg;
      msg << "### FATAL ERROR: EquilibrateTP doesn't converge." << std::endl
          << "Try to increase max_iter." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    return peos->Entropy(prim) - entropy;
  }
};

struct IsoenergeticSolver {
  EquationOfState *peos;
  ReactionGroup *prg;
  Real energy;
  Real *prim;
  int max_iter;
  char c;

  IsoenergeticSolver() : max_iter(10), c('T') {}
  Real operator()(Real v) {
    if (c == 'T') prim[IT] = v;
    else prim[IPR] = v;
    prg->EquilibrateTP(prim, max_iter);
    return peos->Energy(prim) - energy;
  }
};

struct IsoenthalpicSolver {
  EquationOfState *peos;
  ReactionGroup *prg;
  Real enthalpy;
  Real *prim;
  int max_iter;
  char c;

  IsoenthalpicSolver() : max_iter(10), c('T') {}
  Real operator()(Real v) {
    if (c == 'T') prim[IT] = v;
    else prim[IPR] = v;
    prg->EquilibrateTP(prim, max_iter);
    return peos->Enthalpy(prim) - enthalpy;
  }
};

#endif
