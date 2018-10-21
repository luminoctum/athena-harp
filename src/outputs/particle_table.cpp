//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file particle_table.cpp
//  \brief writes output data as a particle table. Writes one file per ParticleGroup per Meshblock.

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
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "outputs.hpp"
#include "../particle/particle.hpp"

//----------------------------------------------------------------------------------------
// ParticleTableOutput constructor
// destructor not required for this derived class

ParticleTableOutput::ParticleTableOutput(OutputParameters oparams)
  : FormattedTableOutput(oparams)
{
}

//----------------------------------------------------------------------------------------
//! \fn void ParticleTableOutput:::WriteOutputFile(Mesh *pm)
//  \brief writes OutputData to file in tabular format using C style fprintf
//         Writes one file per ParticleGroup per MeshBlock

void ParticleTableOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
{
  MeshBlock *pmb=pm->pblock;

  // Loop over MeshBlocks
  while (pmb != NULL) {
    ParticleGroup *ppg = pmb->ppg;
    // Loop over ParticleGroup
    while (ppg != NULL) {
      // create filename: "file_basename"+"."+"name"+"."+"blockid"+"."+?????+".pat",
      // where ????? = 5-digit file_number
      std::string fname;
      char number[6];
      sprintf(number,"%05d",output_params.file_number);
      char blockid[12];
      sprintf(blockid,"block%d",pmb->gid);
  
      fname.assign(output_params.file_basename);
      fname.append(".");
      fname.append(ppg->myname);
      fname.append(".");
      fname.append(blockid);
      fname.append(".");
      fname.append(number);
      fname.append(".pat");

      // open file for output
      FILE *pfile;
      std::stringstream msg;
      if ((pfile = fopen(fname.c_str(),"w")) == NULL){
        msg << "### FATAL ERROR in function [ParticleTableOutput::WriteOutputFile]"
            <<std::endl<< "Output file '" <<fname<< "' could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }

      // print file header
      fprintf(pfile,"# Athena++ data at time=%e",pm->time);
      fprintf(pfile,"  cycle=%d",pmb->pmy_mesh->ncycle);
      fprintf(pfile,"  particle=%s",ppg->myname.c_str());
      fprintf(pfile,"  number of particles=%ld \n",ppg->q.size());

      // write x1, x2, x3 column headers
      fprintf(pfile,"#");
      fprintf(pfile,"%15s","t");
      fprintf(pfile,"%16s","x1");
      fprintf(pfile,"%16s","x2");
      fprintf(pfile,"%16s","x3");
      for (int j = 0; j < NREAL_PARTICLE_DATA; ++j)
        fprintf(pfile, "%14s%02d", "RDATA", j);
      for (int j = 0; j < NINT_PARTICLE_DATA; ++j)
        fprintf(pfile, "%14s%02d", "IDATA", j);
      fprintf(pfile,"\n"); // terminate line

      // loop over all particles
      for (size_t i = 0; i < ppg->q.size(); ++i) {
        fprintf(pfile, output_params.data_format.c_str(), ppg->q[i].time);
        fprintf(pfile, output_params.data_format.c_str(), ppg->q[i].x1);
        fprintf(pfile, output_params.data_format.c_str(), ppg->q[i].x2);
        fprintf(pfile, output_params.data_format.c_str(), ppg->q[i].x3);
        #if NREAL_PARTICLE_DATA > 0
          for (size_t j = 0; j < NREAL_PARTICLE_DATA; ++j)
            fprintf(pfile, output_params.data_format.c_str(), ppg->q[i].rdata[j]);
        #endif

        #if NINT_PARTICLE_DATA > 0
          for (size_t j = 0; j < NINT_PARTICLE_DATA; ++j)
            fprintf(pfile, output_params.data_format.c_str(), ppg->q[i].idata[j]);
        #endif
        fprintf(pfile,"\n"); // terminate line
      }

      // close file, and next variable
      fclose(pfile);
      ppg = ppg->next;
    } // end loop over ParticleGroup
    pmb = pmb->next;
  }  // end loop over MeshBlocks

  // increment counters
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);

  return;
}
