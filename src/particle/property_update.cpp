// Athena++ headers
#include "particle.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../math_funcs.hpp" // _interpn

void ParticleGroup::PropertyUpdate(Real time, Real dt, AthenaArray<Real> const& prim,
  AthenaArray<Real> const& cons, AthenaArray<Real>& cons_out)
{
  AthenaArray<Real> v1, v2, v3;
  Real loc[3];
  int  kji[3] = {0, 0, 0};

  MeshBlock *pmb = pmy_block;

  v1.InitWithShallowSlice(prim,4,IM1,1);
  v2.InitWithShallowSlice(prim,4,IM2,1);
  v3.InitWithShallowSlice(prim,4,IM3,1);

  Real x1min = pmb->block_size.x1min;
  Real x1max = pmb->block_size.x1max;
  Real x2min = pmb->block_size.x2min;
  Real x2max = pmb->block_size.x2max;
  Real x3min = pmb->block_size.x3min;
  Real x3max = pmb->block_size.x3max;

  size_t i = 0, j = q.size();
  int ox1 = 0, ox2 = 0, ox3 = 0, fi1 = 0, fi2 = 0;
  Particle tmp;

  while (i < j) {
    loc[0] = q[i].x3;
    loc[1] = q[i].x2;
    loc[2] = q[i].x1;

    _interpn(&q[i].v1, loc, v1.data(), coordinates_, lengths_, 3);
    if (lengths_[1] > 1)
      _interpn(&q[i].v2, loc, v2.data(), coordinates_, lengths_, 3);
    else
      q[i].v2 = 0.;
    if (lengths_[0] > 1)
      _interpn(&q[i].v3, loc, v3.data(), coordinates_, lengths_, 3);
    else
      q[i].v3 = 0.;

    if (lengths_[0] > 1)
      kji[0] = pmb->ks + _locate(coordinates_, loc[0], lengths_[0]);
    if (lengths_[1] > 1)
      kji[1] = pmb->js + _locate(coordinates_ + lengths_[0], loc[1], lengths_[1]);
    kji[2] = pmb->is + _locate(coordinates_ + lengths_[0] + lengths_[1], loc[2], lengths_[2]);

    bool alive;
    if (particle_fn_ != NULL)
      alive = particle_fn_(pmb, time, dt, q[i], prim, cons, kji, cons_out);
    else
      alive = true;

    // take care of reflective boundary condition
    if (pmb->block_bcs[0] == REFLECTING_BNDRY && q[i].x1 < x1min) {
      q[i].x1 = 2*x1min - q[i].x1;
      q[i].v1 = - q[i].v1;
    }
    if (pmb->block_bcs[1] == REFLECTING_BNDRY && q[i].x1 > x1max) {
      q[i].x1 = 2*x1max - q[i].x1;
      q[i].v1 = - q[i].v1;
    }
    ox1 = q[i].x1 < x1min ? -1 : (q[i].x1 > x1max ? 1 : 0);

    if (pmb->block_size.nx2 > 1) {
      if (pmb->block_bcs[2] == REFLECTING_BNDRY && q[i].x2 < x2min) {
        q[i].x2 = 2*x2min - q[i].x2;
        q[i].v2 = - q[i].v2;
      }
      if (pmb->block_bcs[3] == REFLECTING_BNDRY && q[i].x2 > x2max) {
        q[i].x2 = 2*x2max - q[i].x2;
        q[i].v2 = - q[i].v2;
      }
      ox2 = q[i].x2 < x2min ? -1 : (q[i].x2 > x2max ? 1 : 0);
    }

    if (pmb->block_size.nx3 > 1) {
      if (pmb->block_bcs[4] == REFLECTING_BNDRY && q[i].x3 < x3min) {
        q[i].x3 = 2*x3min - q[i].x3;
        q[i].v3 = - q[i].v3;
      }
      if (pmb->block_bcs[5] == REFLECTING_BNDRY && q[i].x3 > x3max) {
        q[i].x3 = 2*x3max - q[i].x3;
        q[i].v3 = - q[i].v3;
      }
      ox3 = q[i].x3 < x3min ? -1 : (q[i].x3 > x3max ? 1 : 0);
    }

    if (pmb->pmy_mesh->multilevel) {
      // reserved implementation for multilevel, fi1, fi2
    }

    int id = FindBufferID(ox1, ox2, ox3, fi1, fi2, pmb->pmy_mesh->maxneighbor_);

    if (alive && (id == -1)) { // particle is alive and inside domain
      i++;
    } else {  // particle deseased or moved out of the domain
      tmp = q[i];
      q[i] = q[j - 1];
      q[j - 1] = tmp;
      bufid.push_back(alive ? id : -1); // Note that bufid is reversed
      j--;
    }
  }
}
