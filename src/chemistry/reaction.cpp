// C++ headers
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cstdlib>

// Athena++ headers
#include "molecule.hpp"
#include "reaction.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../math_funcs.hpp"

Reaction::Reaction()
{
  for (int i = 0; i < NREACTOR; ++i) {
    reactor[i] = -1;
    measure[i] = 0.;
  }
}

Reaction::~Reaction() {}

Reaction::Reaction(Reaction const& other)
{
  if (this == &other) return;
  *this = other;
}

Reaction& Reaction::operator=(Reaction const& other)
{
  name = other.name;
  tag = other.tag;
  comment = other.comment;
  for (int i = 0; i < NREACTOR; ++i) {
    reactor[i] = other.reactor[i];
    measure[i] = other.measure[i];
  }
  coeff = other.coeff;
  return *this;
}

void Reaction::SetFromString(std::string str, Molecule *pmol, std::string tag_)
{
  std::stringstream msg;
  std::istringstream ss(str);
  std::string token, type;
  bool in_reagent = true;
  bool in_product = false;
  bool in_coeff = false;
  bool in_comment = false;
  int nreagent = 0;
  int ntotal = 0;
  std::string molecule[NREACTOR];

  // clear coefficients
  coeff.clear();

  while (ss.good()) {
    ss >> token;
    if (token == "--" || token == "->") {
      type = token;
      in_reagent = false;
      in_product = true;
      in_coeff = false;
      in_comment = false;
      ss >> token;
    } 
    if (token == "&") {
      in_reagent = false;
      in_product = false;
      in_coeff = true;
      in_comment = false;
      ss >> token;
    }
    if (token == "!") {
      in_reagent = false;
      in_product = false;
      in_coeff = false;
      in_comment = true;
      //ss >> token;
    }
    if (token == "+") continue;
    if (in_reagent) {
      size_t idigit = token.find_first_not_of("0123456789.");
      if (idigit == 0) {
        measure[ntotal] = -1.;
        reactor[ntotal] = pmol->GetMoleculeId(token);
        molecule[ntotal] = token;
      } else {
        std::string digit = token.substr(0, idigit);
        std::string symbol = token.substr(idigit);
        measure[ntotal] = -atof(digit.c_str());
        reactor[ntotal] = pmol->GetMoleculeId(symbol);
        molecule[ntotal] = symbol;
      }
      nreagent++;
      ntotal++;
      if (ntotal >= NREACTOR) {
        msg << "### FATAL ERROR in Reaction:SetFromString, number of reagents exceeds NREACTOR" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
    if (in_product) {
      size_t idigit = token.find_first_not_of("0123456789.");
      if (idigit == 0) {
        measure[ntotal] = 1.;
        reactor[ntotal] = pmol->GetMoleculeId(token);
        molecule[ntotal] = token;
      } else {
        std::string digit = token.substr(0, idigit);
        std::string symbol = token.substr(idigit);
        measure[ntotal] = atof(digit.c_str());
        reactor[ntotal] = pmol->GetMoleculeId(symbol);
        molecule[ntotal] = symbol;
      }
      ntotal++;
      if (ntotal >= NREACTOR) {
        msg << "### FATAL ERROR in Reaction:SetFromString, number of reagents exceeds NREACTOR" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
    if (in_coeff)
      coeff.push_back(atof(token.c_str()));
    if (in_comment) {
      char line[256];
      ss.getline(line, 256);
      comment = line;
      int first = comment.find_first_not_of(" \t"),
          last = comment.find_last_not_of(" \t");
      comment = comment.substr(first, (last - first + 1));
    }
  }
  name = molecule[0] + " ";
  for (int i = 1; i < nreagent; i++)
    name += "+ " + molecule[i] + " ";
  name += type + " " + molecule[nreagent] + " ";
  for (int i = nreagent + 1; i < ntotal; i++)
    name += "+ " + molecule[i] + " ";
  name.erase(name.size() - 1);
  tag = tag_;

  // reset quantities other than ntotal
  for (int i = ntotal; i < NREACTOR; ++i) {
    reactor[i] = -1;
    measure[i] = 0.;
  }
}


std::ostream& operator<<(std::ostream &os, Reaction const& rc)
{
  os << rc.name;
  for (size_t i = 0; i < rc.coeff.size(); ++i)
    os << std::setw(12) << rc.coeff[i];
  if (rc.comment != "")
    os << " ! " << rc.comment;
  if (rc.tag != "")
    os << ", " << rc.tag;
  return os;
}

ReactionGroup::ReactionGroup(MeshBlock *pmb, std::string _name):
  pmy_block(pmb), name(_name), mixr_floor_(1.E-20)
{
  prev = NULL;
  next = NULL;
}

ReactionGroup::~ReactionGroup()
{
  if (prev != NULL) prev->next = next;
  if (next != NULL) next->prev = prev;
}

// functions
ReactionGroup* ReactionGroup::AddReactionGroup(std::string name)
{
  std::stringstream msg;
  ReactionGroup *p = this;
  if (p == NULL) {
    msg << "### FATAL ERROR in AddReactionGroup: ReactionGroup is empty, use new ReactionGroup instead" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  while (p->next != NULL) p = p->next;
  p->next = new ReactionGroup(pmy_block, name);
  p->next->prev = p;
  p->next->next = NULL;

  return p->next;
}

ReactionGroup* ReactionGroup::AddReaction(ParameterInput *pin, std::string block,
  std::string tag, Molecule *pmol, ReactionFunc_t pfunc)
{
  std::string rname = pin->GetString(block, tag);
  Reaction rc;
  rc.SetFromString(rname, pmol, tag);
  fns_.push_back(pfunc);
  rts_.push_back(rc);
}

ReactionGroup* ReactionGroup::GetReactionGroup(std::string name)
{
  std::stringstream msg;
  ReactionGroup *p = this;
  while ((p != NULL) && (p->name != name)) p = p->next;
  return p;
}

Reaction& ReactionGroup::GetReaction(std::string tag)
{
  for (size_t r = 0; r < rts_.size(); ++r)
    if (rts_[r].tag == tag)
      return rts_[r];
  std::stringstream msg;
  msg << "### FATAL ERROR in GetReaction: " << tag << " not found" << std::endl;
  throw std::runtime_error(msg.str().c_str());
}

Reaction const& ReactionGroup::GetReaction(std::string tag) const
{
  for (size_t r = 0; r < rts_.size(); ++r)
    if (rts_[r].tag == tag)
      return rts_[r];
  std::stringstream msg;
  msg << "### FATAL ERROR in GetReaction: " << tag << " not found" << std::endl;
  throw std::runtime_error(msg.str().c_str());
}

void ReactionGroup::SetReactionFunction(std::string tag, ReactionFunc_t pfunc)
{
  for (size_t r = 0; r < rts_.size(); ++r)
    if (rts_[r].tag == tag) {
      fns_[r] = pfunc;
      return;
    }
  std::stringstream msg;
  msg << "### FATAL ERROR in SetReactionFunction: " << tag << " not found" << std::endl;
  throw std::runtime_error(msg.str().c_str());
}

void ReactionGroup::CalculateReactionRates(std::vector<Real>& rates, Real time,
  AthenaArray<Real> const& prim, int i, int j, int k) const
{
  Real prim1[NHYDRO];
  for (int n = 0; n < NHYDRO; ++n)
    prim1[n] = prim(n,k,j,i);
  rates.resize(rts_.size());
  for (size_t r = 0; r < fns_.size(); ++r)
    rates[r] = fns_[r](rts_[r], prim1, time);
}

Real ReactionGroup::EvolveOneTimeStep(AthenaArray<Real>& prim, Real& time, Real dtmax,
  int is, int ie, int js, int je, int ks, int ke)
{
  Real dt = dtmax;
  // allocate space for net reaction rate
  nrate_.resize(NCOMP*(ie-is+1)*(je-js+1)*(ke-ks+1));
  std::fill(nrate_.begin(), nrate_.end(), 0.);
  int count = 0.;
  Real prim1[NHYDRO];

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        // copy primitive variable
        Real x1 = 1.;
        for (int n = 0; n < NHYDRO; ++n) {
          prim1[n] = prim(n,k,j,i);
          if (n > 0 && n < NCOMP) x1 -= prim1[n];
        }
        // calculate net reaction rate for each component
        for (size_t r = 0; r < rts_.size(); ++r) {
          Real rate = fns_[r](rts_[r], prim1, time);
          for (int n = 0; rts_[r].reactor[n] != -1; ++n)
            nrate_[count + rts_[r].reactor[n]] += rate*rts_[r].measure[n];
        }
        // calculate time step;
        if (nrate_[count] < 0.)
          dt = _min(dt, x1/fabs(nrate_[count]));
        for (int n = 1; n < NCOMP; ++n)
          if (nrate_[count+n] < 0.)
            dt = _min(dt, prim1[n]/fabs(nrate_[count+n]));
        count += NCOMP;
      }

  // evolve chemical system for one time step
  count = 0;
  Real norm = 0.;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        ++count; // skip the first one, which is temperature
        for (int n = 1; n < NCOMP; ++n, ++count) {
          prim(n,k,j,i) += nrate_[count]*dt;
          norm += fabs(nrate_[count]*dt);
        }
      }
  time += dt;
  norm /= NCOMP*(ie-is+1)*(je-js+1)*(ke-ks+1);
  return norm;
}

int ReactionGroup::EquilibrateTP(Real prim[], int max_iter) const
{
  int iter = 1;
  Real max_rate;

  do {
    max_rate = 0.;
    for (size_t r = 0; r < rts_.size(); ++r) {
      Real rate = fns_[r](rts_[r], prim, 0.);
      if (fabs(rate) > TINY_NUMBER)
        for (int n = 0; rts_[r].reactor[n] != -1; ++n)
          prim[rts_[r].reactor[n]] += rate*rts_[r].measure[n];
      max_rate = std::max(max_rate, fabs(rate));
    }
    iter++;
  } while (max_rate > TINY_NUMBER && iter < max_iter);

  // apply mixing ratio floor for gases to prevent NAN when evaluating entropy
  for (int n = 1; n < NGAS; ++n)
    prim[n] = std::max(prim[n], mixr_floor_);

  if (iter == max_iter) return 1; // not converged
  else return 0;  // converged and reacted
}

Real NullReaction(Reaction const& rc, Real const prim[], Real time)
{
  return 0.;
}
