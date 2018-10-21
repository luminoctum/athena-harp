#ifndef PARTICLE_HPP
#define PARTICLE_HPP

// C++ headers
#include <vector>
#include <string>
#include <iosfwd>

// Athena++ classes headers
#include "../athena.hpp"

class MeshBlock;
template<typename T> class AthenaArray;
//class ParticleTableOutput;

struct Particle {
  Real time, x1, x2, x3;
  Real v1, v2, v3;

  #if NREAL_PARTICLE_DATA > 0
    Real rdata[NREAL_PARTICLE_DATA];
  #endif

  #if NINT_PARTICLE_DATA > 0
    int  idata[NINT_PARTICLE_DATA];
  #endif

  Particle();
  ~Particle();
  Particle(Particle const& other);
  Particle& operator=(Particle const& other);
};
std::ostream& operator<<(std::ostream &os, Particle const& pt);

class ParticleGroup {
  //friend class ParticleTableOutput;
public:
  ParticleGroup(MeshBlock *pmb, std::string name, ParticleUpdateFunc_t func = NULL);
  ~ParticleGroup();
  
  // data
  MeshBlock* pmy_block;
  std::string myname;
  ParticleGroup *prev, *next;
  std::vector<Particle> q;
  std::vector<int> bufid;

  // functions
  ParticleGroup* AddParticleGroup(std::string name, ParticleUpdateFunc_t func = NULL);
  std::vector<Particle>& GetParticle(std::string name);
  std::vector<Particle> const& GetParticle(std::string name) const;
  void PropertyUpdate(Real time, Real dt, AthenaArray<Real> const & prim,
    AthenaArray<Real> const& cons, AthenaArray<Real>& cons_out);
  int TotalNumber();

protected:
  ParticleUpdateFunc_t particle_fn_;
  Real *coordinates_;
  int lengths_[3];
};

#endif
