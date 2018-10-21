// C++ headers
#include <sstream>
#include <stdexcept>
#include <iostream>

// Athena++ headers
#include "particle.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

Particle::Particle() :
  time(0.), x1(0.), x2(0.), x3(0.),
  v1(0.), v2(0.), v3(0.)
{
#if NREAL_PARTICLE_DATA > 0
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i)
    rdata[i] = 0.;
#endif

#if NINT_PARTICLE_DATA > 0
  for (int i = 0; i < NINT_PARTICLE_DATA; ++i)
    idata[i] = 0;
#endif
}

Particle::~Particle() {}

Particle::Particle(Particle const& other)
{
  if (this == &other) return;
  *this = other;
}

Particle& Particle::operator=(Particle const& other)
{
  time = other.time;
  x1 = other.x1;
  x2 = other.x2;
  x3 = other.x3;
  v1 = other.v1;
  v2 = other.v2;
  v3 = other.v3;
#if NREAL_PARTICLE_DATA > 0
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i)
    rdata[i] = other.rdata[i];
#endif

#if NINT_PARTICLE_DATA > 0
  for (int i = 0; i < NINT_PARTICLE_DATA; ++i)
    idata[i] = other.idata[i];
#endif
  return *this;
}

std::ostream& operator<<(std::ostream &os, Particle const& pt)
{
  os << "time: "<< pt.time << std::endl
     << "x1: " << pt.x1 << " v1: " << pt.v1 << std::endl
     << "x2: " << pt.x2 << " v2: " << pt.v2 << std::endl
     << "x3: " << pt.x3 << " v3: " << pt.v3 << std::endl;
#if NREAL_PARTICLE_DATA > 0
  os << "rdata: ";
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i)
    os << pt.rdata[i] << " ";
  os << std::endl;
#endif

#if NINT_PARTICLE_DATA > 0
  os << "idata: ";
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i)
    os << pt.idata[i] << " ";
#endif
  return os;
}

// constructor, initializes data structure and parameters
ParticleGroup::ParticleGroup(MeshBlock *pmb, std::string name, ParticleUpdateFunc_t func):
  pmy_block(pmb), myname(name), particle_fn_(func)
{
  prev = NULL;
  next = NULL;

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  coordinates_ = new Real [ncells3 + ncells2 + ncells1];

  for (int k = 0; k < ncells3; ++k)
    coordinates_[k] = pmb->pcoord->x3v(k);
  for (int j = 0; j < ncells2; ++j)
    coordinates_[ncells3 + j] = pmb->pcoord->x2v(j);
  for (int i = 0; i < ncells1; ++i)
    coordinates_[ncells3 + ncells2 + i] = pmb->pcoord->x1v(i);

  lengths_[0] = ncells3;
  lengths_[1] = ncells2;
  lengths_[2] = ncells1;
}

// destructor
ParticleGroup::~ParticleGroup()
{
  delete[] coordinates_;
  if (prev != NULL) prev->next = next;
  if (next != NULL) next->prev = prev;
}

// functions
ParticleGroup* ParticleGroup::AddParticleGroup(std::string name, ParticleUpdateFunc_t
func)
{
  std::stringstream msg;
  ParticleGroup *p = this;
  if (p == NULL) {
    msg << "### FATAL ERROR in AddParticleGroup: ParticleGroup is empty, use new ParticleGroup instead" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  while (p->next != NULL) p = p->next;
  p->next = new ParticleGroup(pmy_block, name, func);
  p->next->prev = p;
  p->next->next = NULL;

  return p->next;
}

std::vector<Particle>& ParticleGroup::GetParticle(std::string name)
{
  std::stringstream msg;
  ParticleGroup *p = this;

  while ((p != NULL) && (p->myname != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in GetParticle : ParticleGroup " << name << " not found" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return p->q;
}

std::vector<Particle> const& ParticleGroup::GetParticle(std::string name) const
{
  std::stringstream msg;
  ParticleGroup const *p = this;

  while ((p != NULL) && (p->myname != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in GetParticle : ParticleGroup " << name << " not found" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return p->q;
}

int ParticleGroup::TotalNumber()
{
  int ntotal = 1;
  ParticleGroup *p = this;
  while (p->next != NULL) {
    p = p->next;
    ntotal++;
  }
  return ntotal;
}
