#include <iostream>
#include "../chemistry/molecule.hpp"
#include "../misc.hpp"
#include "absorber.hpp"

std::ostream& operator<<(std::ostream &os, Absorber const& ab)
{
  os << "name = " << ab.myname << std::endl;
  os << "dependence id = ";
  for (int n = 0; n < NCOMP; ++n) 
    if (ab.id_[n] != -1) os << ab.id_[n] << " ";
  os << std::endl;
  return os;
}

Absorber::Absorber(std::string name): 
  myname(name), prev(NULL), next(NULL)
{
  std::fill(id_, id_ + NCOMP, -1);
}

Absorber::Absorber(std::string name, Molecule *pmol, std::string mols):
  myname(name), prev(NULL), next(NULL)
{
  std::fill(id_, id_ + NCOMP, -1);
  std::vector<std::string> names;
  SplitString(mols, names);
  for (size_t n = 0; n < names.size(); ++n)
    id_[n] = pmol->GetMoleculeId(names[n]);
}

Absorber::~Absorber()
{
  if (prev != NULL) prev->next = next;
  if (next != NULL) next->prev = prev;
}
