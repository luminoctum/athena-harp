// C++ headers
#include <stdexcept>
#include <sstream>
#include <cstdlib>  // atof

// Athena++ headers
#include "../parameter_input.hpp"
#include "../athena_arrays.hpp"
#include "../misc.hpp"

#ifdef USE_DISORT
#include "disort_wrapper.hpp"

DisortWrapper::DisortWrapper(ParameterInput *pin, int nwave, AtmDir dir_):
  dir(dir_)
{
  ds.nlyr = pin->GetInteger("radiation", "nlevel") - 1;
  ds.nstr = pin->GetInteger("radiation", "nstr");
  ds.nmom = pin->GetInteger("radiation", "nmom");
  ds.nphase = pin->GetInteger("radiation", "nphase");
  ds.accur = pin->GetOrAddReal("radiation", "accur", 0.);

  ds.bc.btemp = pin->GetReal("radiation", "btemp");
  ds.bc.ttemp = pin->GetReal("radiation", "ttemp");
  ds.bc.fbeam = pin->GetReal("radiation", "fbeam");
  ds.bc.umu0  = pin->GetReal("radiation", "umu0");
  ds.bc.phi0  = pin->GetReal("radiation", "phi0");
  ds.bc.fisot = pin->GetReal("radiation", "fisot");
  ds.bc.albedo = pin->GetReal("radiation", "albedo");
  ds.bc.temis = pin->GetReal("radiation", "temis");

  ds.flag.ibcnd = pin->GetOrAddBoolean("radiation", "ibcnd", 0);
  ds.flag.usrtau = pin->GetOrAddBoolean("radiation", "usrtau", 0);
  ds.flag.usrang = pin->GetOrAddBoolean("radiation", "usrang", 0);
  ds.flag.lamber = pin->GetOrAddBoolean("radiation", "lamber", 1);
  ds.flag.planck = pin->GetOrAddBoolean("radiation", "planck", 0);
  ds.flag.spher = pin->GetOrAddBoolean("radiation", "spher", 0);
  ds.flag.onlyfl = pin->GetOrAddBoolean("radiation", "onlyfl", 1);
  ds.flag.quiet = pin->GetOrAddBoolean("radiation", "quiet", 1);
  ds.flag.intensity_correction = pin->GetOrAddBoolean("radiation", "intensity_correction", 1);
  ds.flag.old_intensity_correction = pin->GetOrAddBoolean("radiation", "old_intensity_correction", 0);
  ds.flag.general_source = pin->GetOrAddBoolean("radiation", "general_source", 0);
  ds.flag.output_uum = pin->GetOrAddBoolean("radiation", "output_uum", 0);

  for (int i = 0; i < 5; ++i) ds.flag.prnt[i] = 0;

  if (ds.flag.usrtau) {
    std::vector<std::string> utau;
    SplitString(pin->GetString("radiation", "utau"), utau);
    ds.ntau = utau.size();
    for (int i = 0; i < ds.ntau; ++i)
      ds.utau[i] = atof(utau[i].c_str());
  }

  std::vector<std::string> umu, phi;
  if (ds.flag.usrang) {
    SplitString(pin->GetString("radiation", "umu"), umu);
    SplitString(pin->GetString("radiation", "phi"), phi);
    ds.numu = umu.size();
    ds.nphi = phi.size();
  } else {
    ds.numu = 0;
    ds.nphi = 0;
  }

  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds, &ds_out);

  for (int i = 0; i < ds.numu; ++i)
    ds.umu[i] = atof(umu[i].c_str());
  for (int i = 0; i < ds.nphi; ++i)
    ds.phi[i] = atof(phi[i].c_str());

  bc.resize(nwave);
  for (int i = 0; i < nwave; ++i)
    bc[i] = ds.bc;
}

DisortWrapper::~DisortWrapper()
{
  c_disort_state_free(&ds);
  c_disort_out_free(&ds, &ds_out);
}

void DisortWrapper::Run(AthenaArray<Real>& rad, AthenaArray<Real>& uu,
  AthenaArray<Real> const& oppr, Real *wave)
{
  std::stringstream msg;
  if (ds.flag.ibcnd != 0) {
    msg << "### FATAL ERROR in DisortWrapper::Radiance: ";
    msg << "ibcnd != 0" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // run disort
  int nwave = oppr.GetDim2();
  int nlevel = oppr.GetDim1();
  for (int i = 0; i < nwave; ++i) {
    // boundary condition
    ds.bc = bc[i];

    // source function
    if (ds.flag.planck)
      for (int j = 0; j < nlevel; ++j)
        ds.temper[j] = dir ? oppr(IT,i,j) : oppr(IT,i,nlevel - 1 - j);
    ds.wvnmlo = wave[i];
    ds.wvnmhi = wave[i];

    // absorption
    for (int j = 0; j < nlevel - 1; ++j)
      ds.dtauc[j] = dir ? oppr(ITAU,i,j) : oppr(ITAU,i,nlevel - 2 - j);

    // single scatering albedo
    for (int j = 0; j < nlevel; ++j)
      ds.ssalb[j] = dir ? oppr(ISA,i,j) : oppr(ISA,i,nlevel - 1 - j);

    // Legendre coefficients
    //for (int j = 0; j < nlevel; ++j)
    //  for (int k = 0; k <= ds.nmom; ++k)
    //    ds.pmom[j*(ds.nmom + 1) + k] = oppr(IPM + k,i,nlevel - 1 - j);

    c_disort(&ds, &ds_out);

    for (int j = 0; j < ds.ntau; ++j) {
      int j1 = dir ? j : ds.ntau - 1 - j;
      rad(i,j1,0) = ds_out.rad[j].rfldir;
      rad(i,j1,1) = ds_out.rad[j].rfldn;
      rad(i,j1,2) = ds_out.rad[j].flup;
      rad(i,j1,3) = ds_out.rad[j].dfdt;
      rad(i,j1,4) = ds_out.rad[j].uavg;
      rad(i,j1,5) = ds_out.rad[j].uavgdn;
      rad(i,j1,6) = ds_out.rad[j].uavgup;
      rad(i,j1,7) = ds_out.rad[j].uavgso;
    }

    if (!ds.flag.onlyfl) {
      int count = 0;
      for (int j = 0; j < ds.nphi; ++j)
        for (int k = 0; k < ds.ntau; ++k)
          for (int l = 0; l < ds.numu; ++l, ++count)
            uu(i,j,k,l) = ds_out.uu[count];
    }
  }
}

#endif // USE_DISORT
