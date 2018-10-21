//! \file test_particle.cpp
//  \breif test program for particle
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../utils/utils.hpp" // ran2
#include "../particle/particle.hpp"

bool ParticleTranslate(MeshBlock *pmb, Particle &pt, int cid[3], Real const time, Real const dt)
{
  pt.x1 += pt.v1 * dt;
  pt.x2 += pt.v2 * dt;
  return true;
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserParticleUpdateFunction(ParticleTranslate);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.;
        phydro->u(IM1,k,j,i) = 1.;
      }

  ppg = new ParticleGroup(this, "tracer1");
  ppg->AddParticleGroup(this, "tracer2");

  std::vector<Particle>& tr1 = ppg->GetParticle("tracer1");
  std::vector<Particle>& tr2 = ppg->GetParticle("tracer2");

  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;

  long int iseed = -1;

  int ntracer = pin->GetInteger("problem", "ntracer");

  for (int n = 0; n < ntracer; ++n) {
    Particle tracer;
    tracer.time = 0.;
    tracer.x1 = x1min + (x1max - x1min) * ran2(&iseed);
    tracer.x2 = 0.;
    tracer.x3 = 0.;
    tr1.push_back(tracer);
    tr2.push_back(tracer);
  }

  Real size = block_size.x1max - block_size.x1min;

  // test boundary condition
  for (int n = 0; n < ntracer; ++n) {
    Particle tracer;
    tracer.time = 0.;
    if (n < ntracer/2) {
      tracer.x1 = x1min - 0.5 + (x1max - x1min) * ran2(&iseed);
    } else {
      tracer.x1 = x1min + 0.5 + (x1max - x1min) * ran2(&iseed);
    }

    if (tracer.x1 - block_size.x1max > size)
      tracer.x1 -= size;
    if (block_size.x1min - tracer.x1 > size)
      tracer.x1 += size;

    tracer.x2 = 0.;
    tracer.x3 = 0.;
    tr1.push_back(tracer);
    tr2.push_back(tracer);
  }
}
