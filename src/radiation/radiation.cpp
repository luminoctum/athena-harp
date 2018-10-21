// C/C++ headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>  // atof

// Athena++ header
#include "../parameter_input.hpp"
#include "../misc.hpp"
#include "absorber.hpp"
#include "radiation.hpp"

// External library headers
#ifdef NETCDFOUTPUT
  #include <netcdf.h>
#endif

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin, std::string name):
  pmy_block(pmb), myname(name)
{
  std::stringstream msg;

  prev = NULL;
  next = NULL;
  pabs = NULL;
  nlevel = pin->GetInteger("radiation", "nlevel");
  nmom = pin->GetInteger("radiation", "nmom");
  nwave = 1;
  bool kname_is_wave = false;

  std::string kname = pin->GetString("radiation", myname);

  if (kname.size() < 3)
    kname_is_wave = true;
  else if (kname.substr(kname.size() - 3, 3) == ".nc")
    kname_is_wave = false;
  else
    kname_is_wave = true;

  if (kname_is_wave) {
    std::vector<std::string> str;
    SplitString(kname, str);
    nwave = str.size();
    if (nwave < 1) {
      msg << "### FATAL ERROR in Radiation: number of wave less than 1" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    wave = new Real[nwave];
    for (int i = 0; i < nwave; ++i)
      wave[i] = atof(str[i].c_str());
  } else {
  #ifdef NETCDFOUTPUT
    int fileid, dimid, varid;
    nc_open(kname.c_str(), NC_NETCDF4, &fileid);
    nc_inq_dimid(fileid, "Wavenumber", &dimid);
    nc_inq_dimlen(fileid, dimid, (size_t*)&nwave);
  #endif
    if (nwave < 1) {
      msg << "### FATAL ERROR in Radiation: number of wave less than 1" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    wave = new Real[nwave];
  #ifdef NETCDFOUTPUT
    nc_inq_varid(fileid, "Wavenumber", &varid);
    nc_get_var_double(fileid, varid, wave);
  #endif
  }

  level = new Real[nlevel];
  weight_ = new Real[nwave];
  weight_[0] = 0.5*(wave[1] - wave[0]);
  weight_[nwave - 1] = 0.5*(wave[nwave - 1] - wave[nwave - 2]);
  for (int i = 1; i < nwave - 1; ++i)
    weight_[i] = 0.5*(wave[i + 1] - wave[i - 1]);
  
  // IT,ITAU,ISA,IPM(0,1,...,nmom)
  oppr.NewAthenaArray(3 + nmom + 1, nwave, nlevel);

  bool usrang = pin->GetOrAddBoolean("radiation", "usrang", 0);
  bool usrtau = pin->GetOrAddBoolean("radiation", "usrtau", 0);

  if (usrtau)
    ntau = pin->GetStringSize("radiation", "utau");
  else
    ntau = nlevel;
  rad.NewAthenaArray(nwave, ntau, NRADIANT);

  if (usrang) {
    numu = pin->GetStringSize("radiation", "umu");
    nphi = pin->GetStringSize("radiation", "phi");
  } else {
    numu = pin->GetInteger("radiation", "nstr");
    nphi = 1;
  }
  uu.NewAthenaArray(nwave, nphi, ntau, numu);
}

Radiation::~Radiation()
{
  if (prev != NULL) prev->next = next;
  if (next != NULL) next->prev = prev;
  while (pabs != NULL) {
    Absorber *p = pabs;
    pabs = pabs->next;
    delete p;
  }
  delete[] wave;
  delete[] level;
  delete[] weight_;
}

Radiation* Radiation::AddRadiation(MeshBlock *pmb, ParameterInput *pin, std::string name)
{
  Radiation* prad = new Radiation(pmb, pin, name);
  Radiation* p = this;
  while (p->next != NULL) p = p->next;
  p->next = prad;
  p->next->prev = p;
  p->next->next = NULL;
  return p->next;
}

Radiation* Radiation::Next(int n) {
  std::stringstream msg;
  Radiation* p = this;
  for (int i = 0; i < n; ++i) {
    p = p->next;
    if (p == NULL) {
      msg << "### FATAL ERROR in Radiation::Next(n) "
          << "past the end of the linked list." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
  return p;
}

// default way of setting optical properties
void __attribute__((weak)) Radiation::SetOpticalProperties(
    AthenaArray<Real> const& prim, int k, int j, int is, int ie)
{
  std::stringstream msg;
  // check dimension consistency
  if (nlevel != ie - is + 1) {
    msg << "### FATAL ERROR in Radiation::SetOpticalProperties "
        << "nlevel != ie - is + 1" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // set tau, ssalb, pmom, etc...
  std::fill(oppr.data(), oppr.data() + oppr.GetSize(), 0.);
  Absorber *q = pabs;
  Real prim1[NHYDRO];
  while (q != NULL) {
    for (int m = 0; m < nwave; ++m)
      for (int i = 0; i < nlevel; ++i) {
        for (int n = 0; n < NHYDRO; ++n)
          prim1[n] = prim(n,k,j,is + i);
        oppr(ITAU,m,i) += q->AbsorptionCoefficient(wave[m], prim1);
        oppr(IT,m,i) = prim1[IT];
      }
    q = q->next;
  }

  // absorption coefficients -> optical thickness
  for (int m = 0; m < nwave; ++m)
    for (int i = 0; i < nlevel; ++i) {
      if (i != nlevel - 1)
        oppr(ITAU,m,i) = 0.5*(oppr(ITAU,m,i) + oppr(ITAU,m,i + 1))
          *fabs(level[i + 1] - level[i]);
      else
        oppr(ITAU,m,i) = 0.;
    }
}

void Radiation::TotalFlux(AthenaArray<Real>& flux) const
{
  // count how many bands
  int b = 0;
  int max_ntau = 0;
  Radiation const *p = this;
  while (p != NULL) {
    max_ntau = std::max(max_ntau, p->ntau);
    p = p->next;
    b++;
  }

  // check dimension consistancy, reallocate memory if necessary
  if (flux.GetDim1() != b && flux.GetDim2() < max_ntau) {
    flux.DeleteAthenaArray();
    flux.NewAthenaArray(max_ntau, b);
  }

  b = 0;
  p = this;
  while (p != NULL) {
    for (int j = 0; j < p->ntau; ++j)
      flux(j,b) = 0.;
    for (int i = 0; i < p->nwave; ++i)
      for (int j = 0; j < p->ntau; ++j)
        // flup - rfldn - rfldir
        flux(j,b) += p->weight_[i]*(p->rad(i,j,2) - p->rad(i,j,1) - p->rad(i,j,0));
    p = p->next;
    b++;
  }
}

void Radiation::WriteTopFlux(std::string fname) const
{
  std::cout << "Top flux written into file: " << fname << std::endl;
  std::ofstream out(fname.c_str(), std::ios::out);
  out << std::left << std::setw(20) << std::setprecision(8) << "# Wavenumber(cm-1)"
      << std::left << std::setw(20) << "Top rfldir[]"
      << std::left << std::setw(20) << "Top rfldn[]"
      << std::left << std::setw(20) << "Top flup[]"
      << std::endl;

  Radiation const *p = this;
  while (p != NULL) {
    for (int i = 0; i < p->nwave; ++i)
      out << std::left << std::setw(20) << p->wave[i]
          << std::left << std::setw(20) << p->rad(i,0,0)
          << std::left << std::setw(20) << p->rad(i,0,1)
          << std::left << std::setw(20) << p->rad(i,0,2)
          << std::endl;
    p = p->next;
  }
}

void Radiation::WriteOpticalDepth(std::string fname) const 
{
  std::cout << "Optical depth written into file: " << fname << std::endl;
  std::ofstream out(fname.c_str(), std::ios::out);

  Radiation const *p = this;
  while (p != NULL) {
    out << std::left << "# Band: " << p->myname << std::endl;
    out << std::setw(40) << "# Number of wavenumbers : " << p->nwave << std::endl;
    out << std::setw(40) << "# Number of levels: " << p->nlevel << std::endl;
    out << std::right;
    for (int i = 0; i < p->nlevel; ++i)
      out << std::setw(16) << std::setprecision(8) << p->level[i]/1.E3;
    out << std::endl;

    for (int i = 0; i < p->nwave; ++i) {
      out << std::setw(16) << std::setprecision(8) << p->wave[i];
      for (int j = 0; j < p->nlevel - 1; ++j)
        out << std::setw(16) << std::setprecision(8) << p->oppr(ITAU,i,j);
      out << std::endl;
    }
    p = p->next;
  }
}

void Radiation::WriteTopRadiance(std::string fname) const
{
  std::cout << "Top radiance written into file: " << fname << std::endl;
  std::ofstream out(fname.c_str(), std::ios::out);
  out << std::left << std::setw(20) << std::setprecision(8)
      << "# Wavenumber/Frequency" << std::endl;

  Radiation const *p = this;
  while (p != NULL) {
    for (int i = 0; i < p->nwave; ++i) {
      out << std::left << std::setw(20) << p->wave[i];
      for (int j = 0; j < p->nphi; ++j)
        for (int k = 0; k < p->numu; ++k)
          out << std::left << std::setw(20) << p->uu(i,j,0,k);
      out << std::endl;
    }
    p = p->next;
  }
}

void WriteHeatingRate(std::string fname, AthenaArray<Real> const& flux,
      AthenaArray<Real> const& hrate, Real const* level)
{
  int nband = flux.GetDim1();
  int nlevel = flux.GetDim2();

  // print mean intensity, flux divergence and heating rate
  std::cout << "Heating rate written into file: " << fname << std::endl;
  std::ofstream out(fname.c_str(), std::ios::out);
  out  << std::left << std::setw(6) << "Level"
       << std::left << std::setw(12) << "Height [km]";

  for (int b = 0; b < nband; ++b) {
    char bname[80];
    sprintf(bname, "B%d flux [w/m^2]", b + 1);
    out << std::left << std::setw(16) << bname
        << std::left << std::setw(20) << "Heating rate [w/kg]";
  }

  if (nband > 1) {
    out << std::left << std::setw(16) << "Flux [w/m^2]"
        << std::left << std::setw(20) << "Heating rate [w/kg]"
        << std::endl;
  } else {
    out << std::endl;
  }

  for (int i = 0; i < nlevel; ++i) {
    out  << std::left << std::setw(6) << i + 1 
         << std::left << std::setw(12) << level[i]/1.E3;

    Real total_flux = 0., total_hrate = 0.;
    for (int b = 0; b < nband; ++b) {
      out << std::left << std::setw(16) << std::setprecision(8) << flux(i,b)
          << std::left << std::setw(20) << std::setprecision(8) << hrate(i,b);
      total_flux += flux(i,b);
      total_hrate += hrate(i,b);
    }

    if (nband > 1) {
      out << std::left << std::setw(16) << std::setprecision(8) << total_flux
          << std::left << std::setw(20) << std::setprecision(8) << total_hrate 
          << std::endl;
    } else {
      out << std::endl;
    }
  }
}
