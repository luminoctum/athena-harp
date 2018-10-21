// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../globals.hpp"

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // pertubation wavelength
  Real kx = 20.*(PI)/(pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real ky = 20.*(PI)/(pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  long int iseed = -1;

  Real xsrf = pin->GetReal("problem", "xsrf");
  Real xtop = pin->GetReal("problem", "xtop");
  Real psrf = pin->GetReal("problem", "psrf");
  Real tsrf = pin->GetReal("problem", "tsrf");
  Real grav = phydro->psrc->GetG2();
  Real prim1[NCOMP];
  AthenaArray<Real>& prim = phydro->w;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        prim1[1] = xsrf + (xtop - xsrf)*(j-js)/(je-js);
        Real cp = peos->HeatCapacityP(prim1);
        Real mu = peos->Mass(prim1);
        Real dz = pcoord->dx2f(j);
        prim(1,k,j,i)   = prim1[1];
        prim(IT,k,j,i)  = j == js ? tsrf : prim(IT,k,j-1,i) + mu*grav/cp*dz;
        prim(IPR,k,j,i) = j == js ? psrf :
          prim(IPR,k,j-1,i)*pow(prim(IT,k,j,i)/prim(IT,k,j-1,i), cp/Globals::Rgas);
        prim(IVX,k,j,i) = prim(IVZ,k,j,i) = 0.;
        prim(IVY,k,j,i) = 0.05*ran2(&iseed)*(1.0+cos(kx*pcoord->x1v(i)))*(1.0+cos(ky*pcoord->x2v(j)))/4.0;
      }
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je,
  ks, ke);
}
