//! \file sw.cpp
//  \brief global shallow water model

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->u(IDN,k,j,i) = 10.;
        if (pcoord->x2v(j) > 0.)
          phydro->u(IM1,k,j,i) = -10.;
        else
          phydro->u(IM1,k,j,i) = 10.;

        if (pcoord->x1v(i) > 0. && pcoord->x1v(i) < 5.)
          phydro->u(IDN,k,j,i) += 2.;
      }
}
