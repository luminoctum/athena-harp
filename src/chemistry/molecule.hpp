#ifndef MOLECULE_HPP_
#define MOLECULE_HPP_

// C++ headers
#include <string>
#include <iosfwd>

// Athena++ headers
#include "../athena.hpp"

/** @file
 * @brief This file contains class Molecule.
 * **Written By**   : Cheng Li, California Institute of Technology <br>
 * **Contact**      : cli@gps.caltech.edu <br>
 * **Revision**     :
 * - Sep 3 2015, first version.
 * - Dec 4 2015, add a macro to use simple thermodynamic formulation
 * - Dec 9 2015, 1) change Real vector to a single array in m_shomate
 *               2) remove m_antoine, instead using idealized vapor function
 * - Feb 11 2016, define a typedef MoleculeList
 */

/** @brief Molecule stored the necessary information for calculating the
 * thermodynamic properties of an air parcel.
 *
 * This class has the following private members storing the properties:
 * - name, the name of this molecule, solids and liquids should add suffix "(s)"
 *   and "(l)" to the corresponding molecule.
 * - mu, molecular weight [g/mol]
 * - cp, solid or liquid heat capacity [J/(mol K)]
 * - latent, latent heat at triple point [kJ/mol]
 * - entropy, standard entropy [J/(mol K)]
 * - enthalpy, standard enthalpy [kJ/mol]
 * - gibbs, standard gibbs free energy [kJ/mol]
 * - tr, triple point temperature [K]
 * - pr, triple point pressure [pa]
 * - tc, critical point temperature [K]
 * - pc, critical point pressure [bar]
 * - cliq, liquid heat capacity [J/(mol K)]
 * - enliq, vaporization heat of liquids at triple point [kJ/mol]
 * - csld, solid heat capacity [J/(mol K)]
 * - ensld, sublimation heat of solids at triple point [kJ/mol]
 * - beta, slope of idealized saturation vapor pressure curve,
 *   beta = (L(T) - (cp - ci) * T) / (R * Tr)
 * - gamma, correction factor of saturation vapor pressure curve,
 *   gamma = (ci - cp) / R
 *
 * The following members store coefficients in thermodynamic expressions:
 * - shomate, 7 shomate expression coefficients
 * - shomate_sp, temperature separations in shomate expressions
 */

#define NSHOMATE 7 // number of Shomate coefficients
#define MAXSHOMATE 3 // maximum number of shomate expressions

// Forward declarations
class ParameterInput;

class Molecule {
  friend std::ostream& operator<<(std::ostream &os, Molecule const& mol);
public:
  Molecule(std::string name);
  Molecule(ParameterInput *pin);
  virtual ~Molecule();

  // data
  std::string myname;
  Real mu, tr, pr, tc, pc;
  PhaseID phase;  // Gas, Liquid, Solid
  Molecule *prev, *next;

  // functions
  Molecule* AddMolecule(std::string name);
  virtual void LoadChemistryFile(std::string chemfile);
  int TotalNumber();
  virtual Real Cp(Real T) const;
  Real Cv(Real T) const;
  virtual Real Enthalpy(Real T) const;
  virtual Real Entropy(Real T) const;
  virtual Real Latent(Real T) const;

  Real SatVaporPres(Real T) const;
  Real SatVaporTemp(Real P,
    Real Tmin = 50., Real Tmax = 2000., Real precision = 1.E-8) const;
  void SetPhase(PhaseID phase);

  Molecule* GetMolecule(std::string name);
  Molecule const* GetMolecule(std::string name) const;
  int GetMoleculeId(std::string name) const;

protected:
  Real m_cp, m_latent, m_entropy, m_enthalpy, 
       m_gibbs, m_cliq, m_enliq, m_csld, m_ensld,
       m_beta, m_gamma;

  int  m_nshomate;
  Real m_shomate[MAXSHOMATE][NSHOMATE];
  Real m_shomate_sp[MAXSHOMATE + 1];
};

struct SatVaporTempSolver {
  Real pres;
  Molecule const* pmol;
  Real operator()(Real temp) {
    return pmol->SatVaporPres(temp) - pres;
  }
};

#endif
