//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file netcdf.cpp
//  \brief writes output data in NETCDF format.
//  Data is written in RECTILINEAR_GRID geometry, in BINARY format, and in FLOAT type
//  Writes one file per MeshBlock.

// C/C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../coordinates/coordinates.hpp"
#include "outputs.hpp"

// Only proceed if NETCDF output enabled
#ifdef NETCDFOUTPUT

// External library headers
#include <netcdf.h> // nc_[create|close], nc_enddef, nc_strerror,
                    // nc_def_[dim|var], nc_put_var_float, nc_put_vara_float

//----------------------------------------------------------------------------------------
// NetcdfOutput constructor
// destructor - not needed for this derived class

NetcdfOutput::NetcdfOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

//----------------------------------------------------------------------------------------
//! \fn void NetcdfOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes OutputData in NETCDF format, one
//         MeshBlock per file

void NetcdfOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
{
  MeshBlock *pmb=pm->pblock;

  // Loop over MeshBlocks
  while (pmb != NULL) {
    // set start/end array indices depending on whether ghost zones are included
    out_is=pmb->is; out_ie=pmb->ie;
    out_js=pmb->js; out_je=pmb->je;
    out_ks=pmb->ks; out_ke=pmb->ke;
    if (output_params.include_ghost_zones) {
      out_is -= NGHOST; out_ie += NGHOST;
      if (out_js != out_je) {out_js -= NGHOST; out_je += NGHOST;}
      if (out_ks != out_ke) {out_ks -= NGHOST; out_ke += NGHOST;}
    }

    // set ptrs to data in OutputData linked list, then slice/sum as needed
    LoadOutputData(pmb);
    if (TransformOutputData(pmb) == false) {
      ClearOutputData();  // required when LoadOutputData() is used.
      pmb=pmb->next;
      continue;
    } // skip if slice was out of range

    // create filename: "file_basename"+ "."+"blockid"+"."+"file_id"+"."+XXXXX+".nc",
    // where XXXXX = 5-digit file_number
    std::string fname;
    char number[6];
    sprintf(number,"%05d",output_params.file_number);
    char blockid[12];
    sprintf(blockid,"block%d",pmb->gid);

    fname.assign(output_params.file_basename);
    fname.append(".");
    fname.append(blockid);
    fname.append(".");
    fname.append(output_params.file_id);
    fname.append(".");
    fname.append(number);
    fname.append(".nc");

    // 1. open file for output
    std::stringstream msg;
    int err, ifile;

    nc_create(fname.c_str(), NC_CLASSIC_MODEL, &ifile);

    // 2. coordinate structure
    int ncells1 = out_ie - out_is + 1;
    int ncells2 = out_je - out_js + 1;
    int ncells3 = out_ke - out_ks + 1;
    int ncoord1 = ncells1; if (ncells1 > 1) ncoord1++;
    int ncoord2 = ncells2; if (ncells2 > 1) ncoord2++;
    int ncoord3 = ncells3; if (ncells3 > 1) ncoord3++;

    // 2. define coordinate
    int idt, idx1, idx2, idx3, idx1b, idx2b, idx3b;
    // time
    nc_def_dim(ifile, "time", NC_UNLIMITED, &idt);

    nc_def_dim(ifile, "x1", ncells1, &idx1);
    if (ncells1 > 1)
      nc_def_dim(ifile, "x1b", ncoord1, &idx1b);

    nc_def_dim(ifile, "x2", ncells2, &idx2);
    if (ncells2 > 1)
      nc_def_dim(ifile, "x2b", ncoord2, &idx2b);

    nc_def_dim(ifile, "x3", ncells3, &idx3);
    if (ncells3 > 1)
      nc_def_dim(ifile, "x3b", ncoord3, &idx3b);

    // 3. define variables
    int ivt, ivx1, ivx2, ivx3, ivx1b, ivx2b, ivx3b;
    int loc[4] = {pmb->loc.lx1, pmb->loc.lx2, pmb->loc.lx3, pmb->loc.level};
    int pos[4];

    nc_def_var(ifile, "time", NC_FLOAT, 1, &idt, &ivt);

    nc_def_var(ifile, "x1", NC_FLOAT, 1, &idx1, &ivx1);
    pos[0] = 1; pos[1] = pmb->pmy_mesh->mesh_size.nx1;
    pos[2] = ncells1 * loc[0] + 1; pos[3] = ncells1 * (loc[0] + 1);
    nc_put_att_int(ifile, ivx1, "domain_decomposition", NC_INT, 4, pos);
    if (ncells1 > 1) {
      nc_def_var(ifile, "x1b", NC_FLOAT, 1, &idx1b, &ivx1b);
      pos[0]--; pos[2]--;
      nc_put_att_int(ifile, ivx1b, "domain_decomposition", NC_INT, 4, pos);
    }

    nc_def_var(ifile, "x2", NC_FLOAT, 1, &idx2, &ivx2);
    pos[0] = 1; pos[1] = pmb->pmy_mesh->mesh_size.nx2;
    pos[2] = ncells2 * loc[1] + 1; pos[3] = ncells2 * (loc[1] + 1);
    nc_put_att_int(ifile, ivx2, "domain_decomposition", NC_INT, 4, pos);
    if (ncells2 > 1) {
      nc_def_var(ifile, "x2b", NC_FLOAT, 1, &idx2b, &ivx2b);
      pos[0]--; pos[2]--;
      nc_put_att_int(ifile, ivx2b, "domain_decomposition", NC_INT, 4, pos);
    }

    nc_def_var(ifile, "x3", NC_FLOAT, 1, &idx3, &ivx3);
    pos[0] = 1; pos[1] = pmb->pmy_mesh->mesh_size.nx3;
    pos[2] = ncells3 * loc[2] + 1; pos[3] = ncells3 * (loc[2] + 1);
    nc_put_att_int(ifile, ivx3, "domain_decomposition", NC_INT, 4, pos);
    if (ncells3 > 1) {
      nc_def_var(ifile, "x3b", NC_FLOAT, 1, &idx3b, &ivx3b);
      pos[0]--; pos[2]--;
      nc_put_att_int(ifile, ivx3b, "domain_decomposition", NC_INT, 4, pos);
    }

    OutputData *pdata = pfirst_data_;

    // count total variables (vector variables are expanded into flat scalars)
    int total_vars = 0;
    while (pdata != NULL) {
      if (pdata->type == "SCALARS")
        ++total_vars;
      else  // VECTORS
        total_vars += pdata->data.GetDim4();
      
      pdata = pdata->pnext;
    }

    int iaxis[4] = {idt, idx3, idx2, idx1};
    int *var_ids = new int [total_vars];
    int *ivar = var_ids;

    pdata = pfirst_data_;
    while (pdata != NULL) {
      if (pdata->type == "SCALARS") {
        nc_def_var(ifile, pdata->name.c_str(), NC_FLOAT, 4, iaxis, ivar++);
      } else { // VECTORS
        for (int n = 0; n < pdata->data.GetDim4(); ++n) {
          char c[16]; sprintf(c, "%d", n + 1);
          std::string name = pdata->name + c;
          nc_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis, ivar++);
        }
      }

      pdata = pdata->pnext;
    }

    nc_enddef(ifile);

    // 4. write variables
    float *data = new float[ncells1 * ncells2 * ncells3];
    size_t start[4] = {0, 0, 0, 0};
    size_t count[4] = {1, ncells3, ncells2, ncells1};

    float time = (float)pm->time;
    nc_put_vara_float(ifile, ivt, start, count, &time);

    for (int i = out_is; i <= out_ie; ++i)
      data[i-out_is] = (float)(pmb->pcoord->x1v(i));
    nc_put_var_float(ifile, ivx1, data);

    if (ncells1 > 1) {
      for (int i = out_is; i <= out_ie + 1; ++i)
        data[i-out_is] = (float)(pmb->pcoord->x1f(i));
      nc_put_var_float(ifile, ivx1b, data);
    }

    for (int j = out_js; j <= out_je; ++j)
      data[j-out_js] = (float)(pmb->pcoord->x2v(j));
    nc_put_var_float(ifile, ivx2, data);

    if (ncells2 > 1) {
      for (int j = out_js; j <= out_je + 1; ++j)
        data[j-out_js] = (float)(pmb->pcoord->x2f(j));
      nc_put_var_float(ifile, ivx2b, data);
    }

    for (int k = out_ks; k <= out_ke; ++k)
      data[k-out_ks] = (float)(pmb->pcoord->x3v(k));
    nc_put_var_float(ifile, ivx3, data);

    if (ncells3 > 1) {
      for (int k = out_ks; k <= out_ke + 1; ++k)
        data[k-out_ks] = (float)(pmb->pcoord->x1f(k));
      nc_put_var_float(ifile, ivx3b, data);
    }

    ivar = var_ids;
    pdata = pfirst_data_;
    while (pdata != NULL) {
      int nvar = pdata->data.GetDim4();
      for (int n = 0; n < nvar; n++) {
        float *it = data;
        for (int k = out_ks; k <= out_ke; ++k)
          for (int j = out_js; j <= out_je; ++j)
            for (int i = out_is; i <= out_ie; ++i, ++it)
              *it = (float)pdata->data(n, k, j, i);
        nc_put_vara_float(ifile, *(ivar++), start, count, data);
      }

      pdata = pdata->pnext;
    }

    // 5. close nc file
    nc_close(ifile);

    ClearOutputData();  // required when LoadOutputData() is used.
    delete [] data;
    delete [] var_ids;

    pmb=pmb->next;
  }  // end loop over MeshBlocks

  // increment counters
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);

  return;
}

#endif // NETCDFOUTPUT
