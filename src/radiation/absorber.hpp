#ifndef ABSORBER_HPP_
#define ABSORBER_HPP_

// C++ header
#include <string>
#include <vector>
#include <iosfwd>

// Athena++ headers
#include "../athena.hpp"

/**@file
 * @brief This file contains declaration of Absorber
*
* **Author** : Cheng Li, California Institute of Technology <br>
* **Contact** : cli@gps.caltech.edu <br>
* **Revision history** :
* - June 21 2016, start documenting this file
* - July 28 2016, merge scatter into absorber
* - June 24 2017, adapt to Athena++ framework
*/

class Molecule;
template<typename T> class AthenaArray;

class Absorber {
  friend std::ostream& operator<<(std::ostream &os, Absorber const& ab);
public:
  // data
  std::string myname;
  Absorber *prev, *next;
  
  // functions
  Absorber(std::string name = "");
  Absorber(std::string name, Molecule *pmol, std::string mols = "");
  virtual ~Absorber();
  template<typename A> Absorber* AddAbsorber(A const& a);
  virtual void SaveCoefficient(std::string fname) const {}
  virtual void LoadCoefficient(std::string fname) {}
  virtual Real AbsorptionCoefficient(Real wave, Real const prim[]) const { return 0.; }
  virtual Real SingleScateringAlbedo(Real wave, Real const prim[]) const { return 0.; }
  virtual void Momentum(Real wave, Real const prim[], Real *pp, int np) const {}

protected:
  int id_[NCOMP];  /**< id of dependent molecules */
};

template<typename A>
inline Absorber* Absorber::AddAbsorber(A const& a) {
  A* pa = new A(a);
  Absorber *p = this;
  while (p->next != NULL) p = p->next;
  p->next = pa;
  p->next->prev = p;
  p->next->next = NULL;
  return p->next;
}

class HitranAbsorber: public Absorber {
  friend std::ostream& operator<<(std::ostream &os, HitranAbsorber const& ab);
  //friend HitranAbsorber const& MakeCKAbsorber<>(HitranAbsorber const& albl,
  //  int const *ck_axis, Real const *ck_wave, int nbins);
public:
  HitranAbsorber() {}
  HitranAbsorber(std::string name, Molecule *pmol) : Absorber(name, pmol, name) {}
  virtual ~HitranAbsorber() {}
  void LoadCoefficient(std::string fname);
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;

protected:
  int len_[3];                  /**< length of interpolation axis */
  std::vector<Real> axis_;    /**< interpolation axis */
  std::vector<Real> kcoeff_;  /**< absorption coefficient */
  AthenaArray<Real>   refatm_;  /**< reference atmosphere */
  Real RefTemp_(Real pres) const;
};

class CIAXiz: public Absorber {
public:
  CIAXiz(std::string name, Molecule *pmol) : Absorber(name, pmol, name) {}
  virtual ~CIAXiz() {}
  void LoadCoefficient(std::string fname);
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;

protected:
  int len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};

class CIAOrton: public Absorber {
public:
  CIAOrton(std::string name, Molecule *pmol) : Absorber(name, pmol, name) {}
  virtual ~CIAOrton() {}
  void LoadCoefficient(std::string fname);
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;

protected:
  int len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};

class N2N2CIA: public Absorber {
public:
  N2N2CIA(Molecule *pmol) : Absorber("N2_N2", pmol, "N2") {}
  virtual ~N2N2CIA() {}
  void LoadCoefficient(std::string fname);
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;

protected:
  int rt_len_[2];
  int fd1_len_[2];
  int fd2_len_[2];
  std::vector<Real> rt_axis_;
  std::vector<Real> fd1_axis_;
  std::vector<Real> fd2_axis_;
  std::vector<Real> rt_;
  std::vector<Real> fd1_;
  std::vector<Real> fd2_;
};

class O2O2CIA: public Absorber {
public:
  O2O2CIA(Molecule *pmol) : Absorber("O2_O2", pmol, "O2") {}
  virtual ~O2O2CIA() {}
  void LoadCoefficient(std::string fname);
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;

protected:
  int fd_len_[2];
  int a1dg_x3sg00_len_[2];
  int a1dg_x3sg10_len_[2];
  int ab_len_[2];
  int other_len_[2];
  std::vector<Real> fd_axis_;
  std::vector<Real> a1dg_x3sg00_axis_;
  std::vector<Real> a1dg_x3sg10_axis_;
  std::vector<Real> ab_axis_;
  std::vector<Real> other_axis_;
  std::vector<Real> fd_;
  std::vector<Real> a1dg_x3sg00_;
  std::vector<Real> a1dg_x3sg10_;
  std::vector<Real> ab_;
  std::vector<Real> other_;
};

class MwrAbsorberCIA: public Absorber {
public:
  MwrAbsorberCIA(Molecule *pmol) : Absorber("CIA", pmol, "H2 He CH4") {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
};

class MwrAbsorberNH3Hanley: public Absorber {
public:
  MwrAbsorberNH3Hanley(Molecule *pmol) : Absorber("NH3", pmol, "H2 He NH3 H2O") {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
};

class MwrAbsorberNH3Bellotti: public Absorber {
public:
  MwrAbsorberNH3Bellotti(Molecule *pmol) : Absorber("NH3", pmol, "H2 He NH3 H2O") {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
};

class MwrAbsorberPH3Hoffman: public Absorber {
public:
  MwrAbsorberPH3Hoffman(Molecule *pmol) : Absorber("PH3", pmol, "H2 He PH3") {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
};


class MwrAbsorberH2OKarpowicz: public Absorber {
public:
  MwrAbsorberH2OKarpowicz(Molecule *pmol) : Absorber("H2O", pmol, "H2 He H2O") {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
};

class MwrAbsorberH2S: public Absorber {
public:
  MwrAbsorberH2S(Molecule *pmol) : Absorber("H2S", pmol, "H2 He H2S") {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
};

class MwrAbsorberFreeFree: public Absorber {
public:
  MwrAbsorberFreeFree() : Absorber("free-free") {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
};

#endif
