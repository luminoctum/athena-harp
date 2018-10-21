//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh.cpp
//  \brief implementation of functions in Mesh class

// C/C++ headers
#include <cfloat>     // FLT_MAX
#include <cmath>      // std::abs(), pow()
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <algorithm>  // sort
#include <iomanip>
#include <stdlib.h>
#include <string.h>  // memcpy

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp" 
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../parameter_input.hpp"
#include "../outputs/io_wrapper.hpp"
#include "../utils/buffer_utils.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../particle/particle.hpp"
#include "mesh_refinement.hpp"
#include "meshblock_tree.hpp"
#include "mesh.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
MPI_Datatype MPI_PARTICLE;
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
// Mesh constructor, builds mesh at start of calculation using parameters in input file

Mesh::Mesh(ParameterInput *pin, int mesh_test)
{
// MPI PARTICLE
#ifdef MPI_PARALLEL
  int counts[2] = {7 + NREAL_PARTICLE_DATA, NINT_PARTICLE_DATA};
  MPI_Datatype types[2] = {MPI_ATHENA_REAL, MPI_INT};
  MPI_Aint disps[2];

  Particle p;

  MPI_Address(&p.time, disps);
  #if NINT_PARTICLE_DATA > 0
    MPI_Address(&p.idata, disps + 1);
  #endif

  disps[1] -= disps[0];
  disps[0] = 0;

  if (NINT_PARTICLE_DATA > 0)
    MPI_Type_struct(2, counts, disps, types, &MPI_PARTICLE);
  else
    MPI_Type_contiguous(counts[0], types[0], &MPI_PARTICLE);
  MPI_Type_commit(&MPI_PARTICLE);
#endif

  std::stringstream msg;
  RegionSize block_size;
  MeshBlockTree *neibt;
  MeshBlock *pfirst;
  enum BoundaryFlag block_bcs[6];
  int nbmax, dim;

  // mesh test
  if(mesh_test>0) Globals::nranks=mesh_test;

  // read time and cycle limits from input file
  start_time = pin->GetOrAddReal("time","start_time",0.0);
  tlim       = pin->GetReal("time","tlim");
  cfl_number = pin->GetReal("time","cfl_number");
  time = start_time;
  dt   = (FLT_MAX*0.4);
  nbnew=0; nbdel=0;

  nlim = pin->GetOrAddInteger("time","nlim",-1);
  ncycle = 0;
  nint_user_mesh_data_=0;
  nreal_user_mesh_data_=0;
  nuser_history_output_=0;

  // read number of OpenMP threads for mesh
  num_mesh_threads_ = pin->GetOrAddInteger("mesh","num_threads",1);
  if (num_mesh_threads_ < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Number of OpenMP threads must be >= 1, but num_threads=" 
        << num_mesh_threads_ << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // read number of grid cells in root level of mesh from input file.  
  mesh_size.nx1 = pin->GetInteger("mesh","nx1");
  if (mesh_size.nx1 < 4) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx1 must be >= 4, but nx1=" 
        << mesh_size.nx1 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  mesh_size.nx2 = pin->GetInteger("mesh","nx2");
  if (mesh_size.nx2 < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx2 must be >= 1, but nx2=" 
        << mesh_size.nx2 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  mesh_size.nx3 = pin->GetInteger("mesh","nx3");
  if (mesh_size.nx3 < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file nx3 must be >= 1, but nx3=" 
        << mesh_size.nx3 << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (mesh_size.nx2 == 1 && mesh_size.nx3 > 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "In mesh block in input file: nx2=1, nx3=" << mesh_size.nx3 
        << ", 2D problems in x1-x3 plane not supported" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  dim=1;
  if(mesh_size.nx2>1) dim=2;
  if(mesh_size.nx3>1) dim=3;

  // check cfl_number
  if(cfl_number > 1.0 && dim==1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The CFL number must be smaller than 1.0 in 1D simulation" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(cfl_number > 0.5 && dim==2) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The CFL number must be smaller than 0.5 in 2D simulation" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(cfl_number > 1.0/3.0 && dim==3) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The CFL number must be smaller than 1/3 in 3D simulation" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // read physical size of mesh (root level) from input file.  
  mesh_size.x1min = pin->GetReal("mesh","x1min");
  mesh_size.x2min = pin->GetReal("mesh","x2min");
  mesh_size.x3min = pin->GetReal("mesh","x3min");

  mesh_size.x1max = pin->GetReal("mesh","x1max");
  mesh_size.x2max = pin->GetReal("mesh","x2max");
  mesh_size.x3max = pin->GetReal("mesh","x3max");

  if (mesh_size.x1max <= mesh_size.x1min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x1max must be larger than x1min: x1min=" << mesh_size.x1min 
        << " x1max=" << mesh_size.x1max << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (mesh_size.x2max <= mesh_size.x2min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x2max must be larger than x2min: x2min=" << mesh_size.x2min 
        << " x2max=" << mesh_size.x2max << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (mesh_size.x3max <= mesh_size.x3min) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Input x3max must be larger than x3min: x3min=" << mesh_size.x3min 
        << " x3max=" << mesh_size.x3max << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // read ratios of grid cell size in each direction
  block_size.x1rat = mesh_size.x1rat = pin->GetOrAddReal("mesh","x1rat",1.0);
  block_size.x2rat = mesh_size.x2rat = pin->GetOrAddReal("mesh","x2rat",1.0);
  block_size.x3rat = mesh_size.x3rat = pin->GetOrAddReal("mesh","x3rat",1.0);

  if (std::abs(mesh_size.x1rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x1rat <= 1.1, x1rat=" 
        << mesh_size.x1rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (std::abs(mesh_size.x2rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x2rat <= 1.1, x2rat=" 
        << mesh_size.x2rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (std::abs(mesh_size.x3rat - 1.0) > 0.1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Ratio of cell sizes must be 0.9 <= x3rat <= 1.1, x3rat=" 
        << mesh_size.x3rat << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // read BC flags for each of the 6 boundaries in turn.
  mesh_bcs[INNER_X1] = GetBoundaryFlag(pin->GetOrAddString("mesh","ix1_bc","none"));
  mesh_bcs[OUTER_X1] = GetBoundaryFlag(pin->GetOrAddString("mesh","ox1_bc","none"));
  mesh_bcs[INNER_X2] = GetBoundaryFlag(pin->GetOrAddString("mesh","ix2_bc","none"));
  mesh_bcs[OUTER_X2] = GetBoundaryFlag(pin->GetOrAddString("mesh","ox2_bc","none"));
  mesh_bcs[INNER_X3] = GetBoundaryFlag(pin->GetOrAddString("mesh","ix3_bc","none"));
  mesh_bcs[OUTER_X3] = GetBoundaryFlag(pin->GetOrAddString("mesh","ox3_bc","none"));

  // read MeshBlock parameters
  block_size.nx1 = pin->GetOrAddInteger("meshblock","nx1",mesh_size.nx1);
  if(dim>=2)
    block_size.nx2 = pin->GetOrAddInteger("meshblock","nx2",mesh_size.nx2);
  else
    block_size.nx2=mesh_size.nx2;
  if(dim==3)
    block_size.nx3 = pin->GetOrAddInteger("meshblock","nx3",mesh_size.nx3);
  else
    block_size.nx3=mesh_size.nx3;

  // check consistency of the block and mesh
  if(mesh_size.nx1%block_size.nx1 != 0
  || mesh_size.nx2%block_size.nx2 != 0
  || mesh_size.nx3%block_size.nx3 != 0) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "the mesh must be evenly divisible by the meshblock" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(block_size.nx1 <4 || (block_size.nx2<4 && dim>=2)
     || (block_size.nx3<4 && dim==3)) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "block_size must be larger than or equal to 4 meshes." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // calculate the number of the blocks
  nrbx1=mesh_size.nx1/block_size.nx1;
  nrbx2=mesh_size.nx2/block_size.nx2;
  nrbx3=mesh_size.nx3/block_size.nx3;
  nbmax=(nrbx1>nrbx2)?nrbx1:nrbx2;
  nbmax=(nbmax>nrbx3)?nbmax:nrbx3;


  //initialize user-enrollable functions
  if(mesh_size.x1rat!=1.0)
    use_meshgen_fn_[X1DIR]=true;
  else
    use_meshgen_fn_[X1DIR]=false;
  if(mesh_size.x2rat!=1.0)
    use_meshgen_fn_[X2DIR]=true;
  else
    use_meshgen_fn_[X2DIR]=false;
  if(mesh_size.x3rat!=1.0)
    use_meshgen_fn_[X3DIR]=true;
  else
    use_meshgen_fn_[X3DIR]=false;
  MeshGenerator_[X1DIR]=DefaultMeshGeneratorX1;
  MeshGenerator_[X2DIR]=DefaultMeshGeneratorX2;
  MeshGenerator_[X3DIR]=DefaultMeshGeneratorX3;
  for(int dir=0; dir<6; dir++)
    BoundaryFunction_[dir]=NULL;
  AMRFlag_=NULL;
  UserSourceTerm_=NULL;
  UserTimeStep_=NULL;
  //particle_fn_ = NULL;

  // calculate the logical root level and maximum level
  for (root_level=0; (1<<root_level)<nbmax; root_level++);
  current_level=root_level;

  // create the root grid
  tree.CreateRootGrid(nrbx1,nrbx2,nrbx3,root_level);

  // SMR / AMR: create finer grids here
  multilevel=false;
  adaptive=false;
  if (pin->GetOrAddString("mesh","refinement","none")=="adaptive")
    adaptive=true, multilevel=true;
  else if (pin->GetOrAddString("mesh","refinement","none")=="static")
    multilevel=true;
  if (adaptive==true) {
    max_level = pin->GetOrAddInteger("mesh","numlevel",1)+root_level-1;
    if(max_level > 63) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The number of the refinement level must be smaller than "
          << 63-root_level+1 << "." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  } else {
    max_level = 63;
  }

  InitUserMeshData(pin);

  if(multilevel==true) {
    if(block_size.nx1%2==1 || (block_size.nx2%2==1 && block_size.nx2>1)
                           || (block_size.nx3%2==1 && block_size.nx3>1)) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
      << "The size of MeshBlock must be divisible by 2 in order to use SMR or AMR."
      << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    InputBlock *pib = pin->pfirst_block;
    while (pib != NULL) {
      if (pib->block_name.compare(0,10,"refinement") == 0) {
        RegionSize ref_size;
        ref_size.x1min=pin->GetReal(pib->block_name,"x1min");
        ref_size.x1max=pin->GetReal(pib->block_name,"x1max");
        if(dim>=2) {
          ref_size.x2min=pin->GetReal(pib->block_name,"x2min");
          ref_size.x2max=pin->GetReal(pib->block_name,"x2max");
        }
        else {
          ref_size.x2min=mesh_size.x2min;
          ref_size.x2max=mesh_size.x2max;
        }
        if(dim>=3) {
          ref_size.x3min=pin->GetReal(pib->block_name,"x3min");
          ref_size.x3max=pin->GetReal(pib->block_name,"x3max");
        }
        else {
          ref_size.x3min=mesh_size.x3min;
          ref_size.x3max=mesh_size.x3max;
        }
        int ref_lev=pin->GetReal(pib->block_name,"level");
        int lrlev=ref_lev+root_level;
        if(lrlev>current_level) current_level=lrlev;
        // range check
        if(ref_lev<1) {
          msg << "### FATAL ERROR in Mesh constructor" << std::endl
              << "Refinement level must be larger than 0 (root level = 0)" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }
        if(lrlev > max_level) {
          msg << "### FATAL ERROR in Mesh constructor" << std::endl
              << "Refinement level exceeds the maximum level (specify maxlevel in <mesh> if adaptive)."
              << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }
        if(ref_size.x1min > ref_size.x1max || ref_size.x2min > ref_size.x2max
        || ref_size.x3min > ref_size.x3max)  {
          msg << "### FATAL ERROR in Mesh constructor" << std::endl
              << "Invalid refinement region is specified."<<  std::endl;
          throw std::runtime_error(msg.str().c_str());
        }
        if(ref_size.x1min < mesh_size.x1min || ref_size.x1max > mesh_size.x1max
        || ref_size.x2min < mesh_size.x2min || ref_size.x2max > mesh_size.x2max
        || ref_size.x3min < mesh_size.x3min || ref_size.x3max > mesh_size.x3max) {
          msg << "### FATAL ERROR in Mesh constructor" << std::endl
              << "Refinement region must be smaller than the whole mesh." << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }
        // find the logical range in the ref_level
        // note: if this is too slow, this should be replaced with bi-section search.
        long int lx1min=0, lx1max=0, lx2min=0, lx2max=0, lx3min=0, lx3max=0;
        long int lxmax=nrbx1*(1L<<ref_lev);
        for(lx1min=0;lx1min<lxmax;lx1min++) {
          if(MeshGenerator_[X1DIR]((Real)(lx1min+1)/lxmax,mesh_size)>ref_size.x1min)
            break;
        }
        for(lx1max=lx1min;lx1max<lxmax;lx1max++) {
          if(MeshGenerator_[X1DIR]((Real)(lx1max+1)/lxmax,mesh_size)>=ref_size.x1max)
            break;
        }
        if(lx1min%2==1) lx1min--;
        if(lx1max%2==0) lx1max++;
        if(dim>=2) { // 2D or 3D
          lxmax=nrbx2*(1L<<ref_lev);
          for(lx2min=0;lx2min<lxmax;lx2min++) {
            if(MeshGenerator_[X2DIR]((Real)(lx2min+1)/lxmax,mesh_size)>ref_size.x2min)
              break;
          }
          for(lx2max=lx2min;lx2max<lxmax;lx2max++) {
            if(MeshGenerator_[X2DIR]((Real)(lx2max+1)/lxmax,mesh_size)>=ref_size.x2max)
              break;
          }
          if(lx2min%2==1) lx2min--;
          if(lx2max%2==0) lx2max++;
        }
        if(dim==3) { // 3D
          lxmax=nrbx3*(1L<<ref_lev);
          for(lx3min=0;lx3min<lxmax;lx3min++) {
            if(MeshGenerator_[X3DIR]((Real)(lx3min+1)/lxmax,mesh_size)>ref_size.x3min)
              break;
          }
          for(lx3max=lx3min;lx3max<lxmax;lx3max++) {
            if(MeshGenerator_[X3DIR]((Real)(lx3max+1)/lxmax,mesh_size)>=ref_size.x3max)
              break;
          }
          if(lx3min%2==1) lx3min--;
          if(lx3max%2==0) lx3max++;
        }
        // create the finest level
        if(dim==1) {
          for(long int i=lx1min; i<lx1max; i+=2) {
            LogicalLocation nloc;
            nloc.level=lrlev, nloc.lx1=i, nloc.lx2=0, nloc.lx3=0;
            int nnew;
            tree.AddMeshBlock(tree,nloc,dim,mesh_bcs,nrbx1,nrbx2,nrbx3,root_level,nnew);
          }
        }
        if(dim==2) {
          for(long int j=lx2min; j<lx2max; j+=2) {
            for(long int i=lx1min; i<lx1max; i+=2) {
              LogicalLocation nloc;
              nloc.level=lrlev, nloc.lx1=i, nloc.lx2=j, nloc.lx3=0;
              int nnew;
              tree.AddMeshBlock(tree,nloc,dim,mesh_bcs,nrbx1,nrbx2,nrbx3,root_level,nnew);
            }
          }
        }
        if(dim==3) {
          for(long int k=lx3min; k<lx3max; k+=2) {
            for(long int j=lx2min; j<lx2max; j+=2) {
              for(long int i=lx1min; i<lx1max; i+=2) {
                LogicalLocation nloc;
                nloc.level=lrlev, nloc.lx1=i, nloc.lx2=j, nloc.lx3=k;
                int nnew;
                tree.AddMeshBlock(tree,nloc,dim,mesh_bcs,nrbx1,nrbx2,nrbx3,root_level,nnew);
              }
            }
          }
        }
      }
      pib=pib->pnext;
    }
  }

  maxneighbor_=BufferID(dim, multilevel);

  // initial mesh hierarchy construction is completed here

  tree.CountMeshBlock(nbtotal);
  loclist=new LogicalLocation[nbtotal];
  tree.GetMeshBlockList(loclist,NULL,nbtotal);

#ifdef MPI_PARALLEL
  // check if there are sufficient blocks
  if(nbtotal < Globals::nranks) {
    if(mesh_test==0) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "Too few mesh blocks: nbtotal ("<< nbtotal <<") < nranks ("
          << Globals::nranks << ")" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    } else { // test
      std::cout << "### Warning in Mesh constructor" << std::endl
          << "Too few mesh blocks: nbtotal ("<< nbtotal <<") < nranks ("
          << Globals::nranks << ")" << std::endl;
    }
  }
#endif

  ranklist=new int[nbtotal];
  nslist=new int[Globals::nranks];
  nblist=new int[Globals::nranks];
  costlist=new Real[nbtotal];
  if(adaptive==true) { // allocate arrays for AMR
    nref = new int [Globals::nranks];
    nderef = new int [Globals::nranks];
    rdisp = new int [Globals::nranks];
    ddisp = new int [Globals::nranks];
    bnref = new int [Globals::nranks];
    bnderef = new int [Globals::nranks];
    brdisp = new int [Globals::nranks];
    bddisp = new int [Globals::nranks];
  }

  // initialize cost array with the simplest estimate; all the blocks are equal
  for(int i=0;i<nbtotal;i++) costlist[i]=1.0;

  LoadBalance(costlist, ranklist, nslist, nblist, nbtotal);

  // Output some diagnostic information to terminal

  // Output MeshBlock list and quit (mesh test only); do not create meshes
  if (mesh_test>0) {
    if (Globals::my_rank==0) OutputMeshStructure(dim);
    return;
  }

  // create MeshBlock list for this process
  int nbs=nslist[Globals::my_rank];
  int nbe=nbs+nblist[Globals::my_rank]-1;
  // create MeshBlock list for this process
  for(int i=nbs;i<=nbe;i++) {
    SetBlockSizeAndBoundaries(loclist[i], block_size, block_bcs);
    // create a block and add into the link list
    if(i==nbs) {
      pblock = new MeshBlock(i, i-nbs, loclist[i], block_size, block_bcs, this, pin);
      pfirst = pblock;
    }
    else {
      pblock->next = new MeshBlock(i, i-nbs, loclist[i], block_size, block_bcs, this, pin);
      pblock->next->prev = pblock;
      pblock = pblock->next;
    }

    pblock->SearchAndSetNeighbors(tree, ranklist, nslist);
  }
  pblock=pfirst;

}

//----------------------------------------------------------------------------------------
// Mesh constructor for restarts. Load the restart file

Mesh::Mesh(ParameterInput *pin, IOWrapper& resfile, int mesh_test)
{
// MPI PARTICLE
#ifdef MPI_PARALLEL
  int counts[2] = {7 + NREAL_PARTICLE_DATA, NINT_PARTICLE_DATA};
  MPI_Datatype types[2] = {MPI_ATHENA_REAL, MPI_INT};
  MPI_Aint disps[2];

  Particle p;

  MPI_Address(&p.time, disps);
  #if NINT_PARTICLE_DATA > 0
    MPI_Address(&p.idata, disps + 1);
  #endif

  disps[1] -= disps[0];
  disps[0] = 0;

  if (NINT_PARTICLE_DATA > 0)
    MPI_Type_struct(2, counts, disps, types, &MPI_PARTICLE);
  else
    MPI_Type_contiguous(counts[0], types[0], &MPI_PARTICLE);
  MPI_Type_commit(&MPI_PARTICLE);
#endif

  std::stringstream msg;
  RegionSize block_size;
  enum BoundaryFlag block_bcs[6];
  MeshBlock *pfirst;
  int i, j, dim;
  IOWrapperSize_t *offset, datasize, listsize, headeroffset;

  // mesh test
  if(mesh_test>0) Globals::nranks=mesh_test;

  // read time and cycle limits from input file
  start_time = pin->GetOrAddReal("time","start_time",0.0);
  tlim       = pin->GetReal("time","tlim");
  cfl_number = pin->GetReal("time","cfl_number");
  nlim = pin->GetOrAddInteger("time","nlim",-1);
  nint_user_mesh_data_=0;
  nreal_user_mesh_data_=0;
  nuser_history_output_=0;

  nbnew=0; nbdel=0;

  // read number of OpenMP threads for mesh
  num_mesh_threads_ = pin->GetOrAddInteger("mesh","num_threads",1);
  if (num_mesh_threads_ < 1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Number of OpenMP threads must be >= 1, but num_threads=" 
        << num_mesh_threads_ << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // read BC flags for each of the 6 boundaries
  mesh_bcs[INNER_X1] = GetBoundaryFlag(pin->GetOrAddString("mesh","ix1_bc","none"));
  mesh_bcs[OUTER_X1] = GetBoundaryFlag(pin->GetOrAddString("mesh","ox1_bc","none"));
  mesh_bcs[INNER_X2] = GetBoundaryFlag(pin->GetOrAddString("mesh","ix2_bc","none"));
  mesh_bcs[OUTER_X2] = GetBoundaryFlag(pin->GetOrAddString("mesh","ox2_bc","none"));
  mesh_bcs[INNER_X3] = GetBoundaryFlag(pin->GetOrAddString("mesh","ix3_bc","none"));
  mesh_bcs[OUTER_X3] = GetBoundaryFlag(pin->GetOrAddString("mesh","ox3_bc","none"));

  // get the end of the header
  headeroffset=resfile.GetPosition();
  // read the restart file
  // the file is already open and the pointer is set to after <par_end>
  IOWrapperSize_t headersize = sizeof(int)*3+sizeof(Real)*2
                             + sizeof(RegionSize)+sizeof(IOWrapperSize_t);
  char *headerdata = new char[headersize];
  if(Globals::my_rank==0) { // the master process reads the header data
    if(resfile.Read(headerdata,1,headersize)!=headersize) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The restart file is broken." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
#ifdef MPI_PARALLEL
  // then broadcast the header data
  MPI_Bcast(headerdata, headersize, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
  IOWrapperSize_t hdos = 0;
  memcpy(&nbtotal, &(headerdata[hdos]), sizeof(int));
  hdos+=sizeof(int);
  memcpy(&root_level, &(headerdata[hdos]), sizeof(int));
  hdos+=sizeof(int);
  current_level=root_level;
  memcpy(&mesh_size, &(headerdata[hdos]), sizeof(RegionSize));
  hdos+=sizeof(RegionSize);
  memcpy(&time, &(headerdata[hdos]), sizeof(Real));
  hdos+=sizeof(Real);
  memcpy(&dt, &(headerdata[hdos]), sizeof(Real));
  hdos+=sizeof(Real);
  memcpy(&ncycle, &(headerdata[hdos]), sizeof(int));
  hdos+=sizeof(int);
  memcpy(&datasize, &(headerdata[hdos]), sizeof(IOWrapperSize_t));
  hdos+=sizeof(IOWrapperSize_t);

  dim=1;
  if(mesh_size.nx2>1) dim=2;
  if(mesh_size.nx3>1) dim=3;

  // check cfl_number
  if(cfl_number > 1.0 && dim==1) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The CFL number must be smaller than 1.0 in 1D simulation" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(cfl_number > 0.5 && dim==2) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The CFL number must be smaller than 0.5 in 2D simulation" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(cfl_number > 1.0/3.0 && dim==3) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The CFL number must be smaller than 1/3 in 3D simulation" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  //initialize
  loclist=new LogicalLocation[nbtotal];
  offset=new IOWrapperSize_t[nbtotal];
  costlist=new Real[nbtotal];
  ranklist=new int[nbtotal];
  nslist=new int[Globals::nranks];
  nblist=new int[Globals::nranks];

  block_size.nx1 = pin->GetOrAddReal("meshblock","nx1",mesh_size.nx1);
  block_size.nx2 = pin->GetOrAddReal("meshblock","nx2",mesh_size.nx2);
  block_size.nx3 = pin->GetOrAddReal("meshblock","nx3",mesh_size.nx3);

  // calculate the number of the blocks
  nrbx1=mesh_size.nx1/block_size.nx1;
  nrbx2=mesh_size.nx2/block_size.nx2;
  nrbx3=mesh_size.nx3/block_size.nx3;

  //initialize user-enrollable functions
  if(mesh_size.x1rat!=1.0)
    use_meshgen_fn_[X1DIR]=true;
  else
    use_meshgen_fn_[X1DIR]=false;
  if(mesh_size.x2rat!=1.0)
    use_meshgen_fn_[X2DIR]=true;
  else
    use_meshgen_fn_[X2DIR]=false;
  if(mesh_size.x3rat!=1.0)
    use_meshgen_fn_[X3DIR]=true;
  else
    use_meshgen_fn_[X3DIR]=false;
  MeshGenerator_[X1DIR]=DefaultMeshGeneratorX1;
  MeshGenerator_[X2DIR]=DefaultMeshGeneratorX2;
  MeshGenerator_[X3DIR]=DefaultMeshGeneratorX3;
  for(int dir=0; dir<6; dir++)
    BoundaryFunction_[dir]=NULL;
  AMRFlag_=NULL;
  UserSourceTerm_=NULL;
  UserTimeStep_=NULL;

  multilevel=false;
  adaptive=false;
  if (pin->GetOrAddString("mesh","refinement","none")=="adaptive")
    adaptive=true, multilevel=true;
  else if (pin->GetOrAddString("mesh","refinement","none")=="static")
    multilevel=true;
  if (adaptive==true) {
    max_level = pin->GetOrAddInteger("mesh","numlevel",1)+root_level-1;
    if(max_level > 63) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The number of the refinement level must be smaller than "
          << 63-root_level+1 << "." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  } else {
    max_level = 63;
  }

  InitUserMeshData(pin);

  // read user Mesh data
  IOWrapperSize_t udsize = 0;
  for(int n=0; n<nint_user_mesh_data_; n++)
    udsize+=iuser_mesh_data[n].GetSizeInBytes();
  for(int n=0; n<nreal_user_mesh_data_; n++)
    udsize+=ruser_mesh_data[n].GetSizeInBytes();
  if(udsize!=0) {
    char *userdata = new char[udsize];
    if(Globals::my_rank==0) { // only the master process reads the ID list
      if(resfile.Read(userdata,1,udsize)!=udsize) {
        msg << "### FATAL ERROR in Mesh constructor" << std::endl
            << "The restart file is broken." << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
#ifdef MPI_PARALLEL
    // then broadcast the ID list
    MPI_Bcast(userdata, udsize, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

    IOWrapperSize_t udoffset=0;
    for(int n=0; n<nint_user_mesh_data_; n++) {
      memcpy(iuser_mesh_data[n].data(), &(userdata[udoffset]),
             iuser_mesh_data[n].GetSizeInBytes());
      udoffset+=iuser_mesh_data[n].GetSizeInBytes();
    }
    for(int n=0; n<nreal_user_mesh_data_; n++) {
      memcpy(ruser_mesh_data[n].data(), &(userdata[udoffset]),
             ruser_mesh_data[n].GetSizeInBytes());
      udoffset+=ruser_mesh_data[n].GetSizeInBytes();
    }
    delete [] userdata;
  }

  // read the ID list
  listsize=sizeof(LogicalLocation)+sizeof(Real);
  //allocate the idlist buffer
  char *idlist = new char [listsize*nbtotal];
  if(Globals::my_rank==0) { // only the master process reads the ID list
    if(resfile.Read(idlist,listsize,nbtotal)!=nbtotal) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "The restart file is broken." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
#ifdef MPI_PARALLEL
  // then broadcast the ID list
  MPI_Bcast(idlist, listsize*nbtotal, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  int os=0;
  for(int i=0;i<nbtotal;i++) {
    memcpy(&(loclist[i]), &(idlist[os]), sizeof(LogicalLocation));
    os+=sizeof(LogicalLocation);
    memcpy(&(costlist[i]), &(idlist[os]), sizeof(Real));
    os+=sizeof(Real);
    if(loclist[i].level>current_level) current_level=loclist[i].level;
  }
  delete [] idlist;

  // calculate the header offset and seek
  headeroffset+=headersize+udsize+listsize*nbtotal;
  if(Globals::my_rank!=0)
    resfile.Seek(headeroffset);

  maxneighbor_=BufferID(dim, multilevel);

  // rebuild the Block Tree
  for(int i=0;i<nbtotal;i++)
    tree.AddMeshBlockWithoutRefine(loclist[i],nrbx1,nrbx2,nrbx3,root_level);
  int nnb;
  // check the tree structure, and assign GID
  tree.GetMeshBlockList(loclist, NULL, nnb);
  if(nnb!=nbtotal) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "Tree reconstruction failed. The total numbers of the blocks do not match. ("
        << nbtotal << " != " << nnb << ")" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

#ifdef MPI_PARALLEL
  if(nbtotal < Globals::nranks) {
    if(mesh_test==0) {
      msg << "### FATAL ERROR in Mesh constructor" << std::endl
          << "Too few mesh blocks: nbtotal ("<< nbtotal <<") < nranks ("
          << Globals::nranks << ")" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    else { // test
      std::cout << "### Warning in Mesh constructor" << std::endl
          << "Too few mesh blocks: nbtotal ("<< nbtotal <<") < nranks ("
          << Globals::nranks << ")" << std::endl;
      return;
    }
  }
#endif

  if(adaptive==true) { // allocate arrays for AMR
    nref = new int [Globals::nranks];
    nderef = new int [Globals::nranks];
    rdisp = new int [Globals::nranks];
    ddisp = new int [Globals::nranks];
    bnref = new int [Globals::nranks];
    bnderef = new int [Globals::nranks];
    brdisp = new int [Globals::nranks];
    bddisp = new int [Globals::nranks];
  }

  LoadBalance(costlist, ranklist, nslist, nblist, nbtotal);

  // Output MeshBlock list and quit (mesh test only); do not create meshes
  if(mesh_test>0) {
    if(Globals::my_rank==0) OutputMeshStructure(dim);
    delete [] offset;
    return;
  }

  // allocate data buffer
  int nb=nblist[Globals::my_rank];
  int nbs=nslist[Globals::my_rank];
  int nbe=nbs+nb-1;
  char *mbdata = new char [datasize*nb];
  // load MeshBlocks (parallel)
  if(resfile.Read_at_all(mbdata, datasize, nb, headeroffset+nbs*datasize)!=nb) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The restart file is broken or input parameters are inconsistent."
        << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  for(i=nbs;i<=nbe;i++) {
    int buff_os = datasize * (i-nbs);
    SetBlockSizeAndBoundaries(loclist[i], block_size, block_bcs);
    // create a block and add into the link list
    if(i==nbs) {
      pblock = new MeshBlock(i, i-nbs, this, pin, loclist[i], block_size,
                             block_bcs, costlist[i], mbdata+buff_os);
      pfirst = pblock;
    }
    else {
      pblock->next = new MeshBlock(i, i-nbs, this, pin, loclist[i], block_size,
                                   block_bcs, costlist[i], mbdata+buff_os);
      pblock->next->prev = pblock;
      pblock = pblock->next;
    }
    pblock->SearchAndSetNeighbors(tree, ranklist, nslist);
  }
  pblock=pfirst;
  delete [] mbdata;
  // check consistency
  if(datasize!=pblock->GetBlockSizeInBytes()) {
    msg << "### FATAL ERROR in Mesh constructor" << std::endl
        << "The restart file is broken or input parameters are inconsistent."
        << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // clean up
  delete [] offset;
}

//----------------------------------------------------------------------------------------
// destructor

Mesh::~Mesh()
{
  while(pblock->prev != NULL) // should not be true
    delete pblock->prev;
  while(pblock->next != NULL)
    delete pblock->next;
  delete pblock;
  delete [] nslist;
  delete [] nblist;
  delete [] ranklist;
  delete [] costlist;
  delete [] loclist;
  if(adaptive==true) { // deallocate arrays for AMR
    delete [] nref;
    delete [] nderef;
    delete [] rdisp;
    delete [] ddisp;
    delete [] bnref;
    delete [] bnderef;
    delete [] brdisp;
    delete [] bddisp;
  }
  // delete user Mesh data
  for(int n=0; n<nreal_user_mesh_data_; n++)
    ruser_mesh_data[n].DeleteAthenaArray();
  if(nreal_user_mesh_data_>0) delete [] ruser_mesh_data;
  for(int n=0; n<nint_user_mesh_data_; n++)
    iuser_mesh_data[n].DeleteAthenaArray();
  if(nint_user_mesh_data_>0) delete [] iuser_mesh_data;

#ifdef MPI_PARALLEL
  MPI_Type_free(&MPI_PARTICLE);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::OutputMeshStructure(int dim)
//  \brief print the mesh structure information

void Mesh::OutputMeshStructure(int dim)
{
  RegionSize block_size;
  enum BoundaryFlag block_bcs[6];
  FILE *fp;

  // open 'mesh_structure.dat' file
  if(dim>=2) {
    if ((fp = fopen("mesh_structure.dat","wb")) == NULL) {
      std::cout << "### ERROR in function Mesh::OutputMeshStructure" << std::endl
                << "Cannot open mesh_structure.dat" << std::endl;
      return;
    }
  }

  // Write overall Mesh structure to stdout and file
  std::cout << std::endl;
  std::cout << "Root grid = " << nrbx1 << " x " << nrbx2 << " x " << nrbx3
            << " MeshBlocks" << std::endl;
  std::cout << "Total number of MeshBlocks = " << nbtotal << std::endl;
  std::cout << "Number of physical refinement levels = " 
            << (current_level - root_level) << std::endl;
  std::cout << "Number of logical  refinement levels = " << current_level << std::endl;

  // compute/output number of blocks per level, and cost per level
  int nb_per_plevel[max_level];
  int cost_per_plevel[max_level];
  for (int i=0; i<=max_level; ++i) {
    nb_per_plevel[i]=0;
    cost_per_plevel[i]=0;
  }
  for (int i=0; i<nbtotal; i++) {
    nb_per_plevel[(loclist[i].level - root_level)]++;
    cost_per_plevel[(loclist[i].level - root_level)] += costlist[i];
  }
  for(int i=root_level;i<=max_level;i++) {
    if(nb_per_plevel[i-root_level]!=0) {
      std::cout << "  Physical level = " << i-root_level << " (logical level = " << i 
                << "): " << nb_per_plevel[i-root_level] << " MeshBlocks, cost = " 
                << cost_per_plevel[i-root_level] <<  std::endl;
    }
  }

  // compute/output number of blocks per rank, and cost per rank
  std::cout << "Number of parallel ranks = " << Globals::nranks << std::endl;
  int nb_per_rank[Globals::nranks];
  int cost_per_rank[Globals::nranks];
  for (int i=0; i<Globals::nranks; ++i) {
    nb_per_rank[i]=0;
    cost_per_rank[i]=0;
  }
  for (int i=0; i<nbtotal; i++) {
    nb_per_rank[ranklist[i]]++;
    cost_per_rank[ranklist[i]] += costlist[i];
  }
  for (int i=0; i<Globals::nranks; ++i) {
    std::cout << "  Rank = " << i << ": " << nb_per_rank[i] <<" MeshBlocks, cost = "
              << cost_per_rank[i] << std::endl;
  }

  // output relative size/locations of meshblock to file, for plotting
  Real mincost=FLT_MAX, maxcost=0.0, totalcost=0.0;
  for (int i=root_level; i<=max_level; i++) {
    Real dx=1.0/(Real)(1L<<i);
    for (int j=0; j<nbtotal; j++) {
      if(loclist[j].level==i) {
        SetBlockSizeAndBoundaries(loclist[j], block_size, block_bcs);
        long int &lx1=loclist[j].lx1;
        long int &lx2=loclist[j].lx2;
        long int &lx3=loclist[j].lx3;
        int &ll=loclist[j].level;
        mincost=std::min(mincost,costlist[i]);
        maxcost=std::max(maxcost,costlist[i]);
        totalcost+=costlist[i];
        fprintf(fp,"#MeshBlock %d on rank=%d with cost=%g\n",j,ranklist[j],costlist[j]);
        fprintf(fp,"#  Logical level %d, location = (%ld %ld %ld)\n",ll,lx1,lx2,lx3);
        if(dim==2) {
          fprintf(fp, "%g %g\n", block_size.x1min, block_size.x2min);
          fprintf(fp, "%g %g\n", block_size.x1max, block_size.x2min);
          fprintf(fp, "%g %g\n", block_size.x1max, block_size.x2max);
          fprintf(fp, "%g %g\n", block_size.x1min, block_size.x2max);
          fprintf(fp, "%g %g\n", block_size.x1min, block_size.x2min);
          fprintf(fp, "\n\n");
        }
        if(dim==3) {
          fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2min, block_size.x3min);
          fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2min, block_size.x3min);
          fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2max, block_size.x3min);
          fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2max, block_size.x3min);
          fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2min, block_size.x3min);
          fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2min, block_size.x3max);
          fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2min, block_size.x3max);
          fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2min, block_size.x3min);
          fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2min, block_size.x3max);
          fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2max, block_size.x3max);
          fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2max, block_size.x3min);
          fprintf(fp, "%g %g %g\n", block_size.x1max, block_size.x2max, block_size.x3max);
          fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2max, block_size.x3max);
          fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2max, block_size.x3min);
          fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2max, block_size.x3max);
          fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2min, block_size.x3max);
          fprintf(fp, "%g %g %g\n", block_size.x1min, block_size.x2min, block_size.x3min);
          fprintf(fp, "\n\n");
        }
      }
    }
  }

  // close file, final outputs
  if(dim>=2) fclose(fp);
  std::cout << "Load Balancing:" << std::endl;
  std::cout << "  Minimum cost = " << mincost << ", Maximum cost = " << maxcost
            << ", Average cost = " << totalcost/nbtotal << std::endl << std::endl;
  std::cout << "See the 'mesh_structure.dat' file for a complete list"
            << " of MeshBlocks." << std::endl;
  std::cout << "Use 'python ../vis/python/plot_mesh.py' or gnuplot"
            << " to visualize mesh structure." << std::endl << std::endl;

  return;
}

//----------------------------------------------------------------------------------------
// \!fn void Mesh::NewTimeStep(void)
// \brief function that loops over all MeshBlocks and find new timestep
//        this assumes that phydro->NewBlockTimeStep is already called

void Mesh::NewTimeStep(void)
{
  MeshBlock *pmb = pblock;
  Real min_dt=pmb->new_block_dt;
  pmb=pmb->next;
  while (pmb != NULL)  {
    min_dt=std::min(min_dt,pmb->new_block_dt);
    pmb=pmb->next;
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,&min_dt,1,MPI_ATHENA_REAL,MPI_MIN,MPI_COMM_WORLD);
#endif
  // set it
  dt=std::min(min_dt,2.0*dt);
  if (time < tlim && tlim-time < dt)  // timestep would take us past desired endpoint
    dt = tlim-time;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserBoundaryFunction(enum BoundaryFace dir, BValHydro_t my_bc)
//  \brief Enroll a user-defined boundary function

void Mesh::EnrollUserBoundaryFunction(enum BoundaryFace dir, BValFunc_t my_bc)
{
  std::stringstream msg;
  if(dir<0 || dir>5) {
    msg << "### FATAL ERROR in EnrollBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(mesh_bcs[dir]!=USER_BNDRY) {
    msg << "### FATAL ERROR in EnrollUserBoundaryFunction" << std::endl
        << "The boundary condition flag must be set to the string 'user' in the "
        << " <mesh> block in the input file to use user-enrolled BCs" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  BoundaryFunction_[dir]=my_bc;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserRefinementCondition(AMRFlagFunc_t amrflag)
//  \brief Enroll a user-defined function for checking refinement criteria

void Mesh::EnrollUserRefinementCondition(AMRFlagFunc_t amrflag)
{
  if(adaptive==true)
    AMRFlag_=amrflag;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserMeshGenerator(enum CoordinateDirection,MeshGenFunc_t my_mg)
//  \brief Enroll a user-defined function for Mesh generation

void Mesh::EnrollUserMeshGenerator(enum CoordinateDirection dir, MeshGenFunc_t my_mg)
{
  std::stringstream msg;
  if(dir<0 || dir>3) {
    msg << "### FATAL ERROR in EnrollUserMeshGenerator function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  use_meshgen_fn_[dir]=true;
  MeshGenerator_[dir]=my_mg;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserExplicitSourceFunction(SrcTermFunc_t my_func)
//  \brief Enroll a user-defined source function

void Mesh::EnrollUserExplicitSourceFunction(SrcTermFunc_t my_func)
{
  UserSourceTerm_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserTimeStepFunction(TimeStepFunc_t my_func)
//  \brief Enroll a user-defined time step function

void Mesh::EnrollUserTimeStepFunction(TimeStepFunc_t my_func)
{
  UserTimeStep_ = my_func;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserParticleUpdateFunction(ParticleUpdateFunc_t my_func)
//  \brief Enroll a user-defined particle update function

//void Mesh::EnrollUserParticleUpdateFunction(ParticleUpdateFunc_t my_func)
//{
//  particle_fn_ = my_func;
//  return;
//}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::AllocateUserHistoryOutput(int n)
//  \brief set the number of user-defined history outputs

void Mesh::AllocateUserHistoryOutput(int n)
{
  nuser_history_output_ = n;
  user_history_output_names_ = new std::string[n];
  user_history_func_ = new HistoryOutputFunc_t[n];
  for(int i=0; i<n; i++) user_history_func_[i] = NULL;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::EnrollUserHistoryOutput(int i, HistoryOutputFunc_t my_func,
//                                         const char *name)
//  \brief Enroll a user-defined history output function and set its name

void Mesh::EnrollUserHistoryOutput(int i, HistoryOutputFunc_t my_func, const char *name)
{
  std::stringstream msg;
  if(i>=nuser_history_output_) {
    msg << "### FATAL ERROR in EnrollUserHistoryOutput function" << std::endl
        << "The number of the user-defined history output (" << i << ") "
        << "exceeds the declared number (" << nuser_history_output_ << ")." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  user_history_output_names_[i] = name;
  user_history_func_[i] = my_func;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::AllocateRealUserMeshDataField(int n)
//  \brief Allocate Real AthenaArrays for user-defned data in Mesh

void Mesh::AllocateRealUserMeshDataField(int n)
{
  if(nreal_user_mesh_data_!=0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::AllocateRealUserMeshDataField"
        << std::endl << "User Mesh data arrays are already allocated" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  nreal_user_mesh_data_=n;
  ruser_mesh_data = new AthenaArray<Real>[n];
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::AllocateIntUserMeshDataField(int n)
//  \brief Allocate integer AthenaArrays for user-defned data in Mesh

void Mesh::AllocateIntUserMeshDataField(int n)
{
  if(nint_user_mesh_data_!=0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::AllocateIntUserMeshDataField"
        << std::endl << "User Mesh data arrays are already allocated" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  nint_user_mesh_data_=n;
  iuser_mesh_data = new AthenaArray<int>[n];
  return;
}

//----------------------------------------------------------------------------------------
// \!fn void Mesh::Initialize(int res_flag, ParameterInput *pin)
// \brief  initialization before the main loop

void Mesh::Initialize(int res_flag, ParameterInput *pin)
{
  MeshBlock *pmb;
  Hydro *phydro;
  Field *pfield;
  BoundaryValues *pbval;
  std::stringstream msg;
  int inb=nbtotal;

  bool iflag=true;
  do {
    if(res_flag==0) {
      pmb = pblock;
      while (pmb != NULL)  {
        pmb->ProblemGenerator(pin);
        pmb->pbval->CheckBoundary();
        pmb=pmb->next;
      }
    }

    // prepare to receive conserved variables
    pmb = pblock;
    while (pmb != NULL)  {
      pmb->pbval->Initialize();
      pmb->pbval->StartReceivingForInit(true);
      pmb=pmb->next;
    }

    // send conserved variables
    pmb = pblock;
    while (pmb != NULL)  {
      phydro=pmb->phydro;
      pfield=pmb->pfield;
      pmb->pbval->SendHydroBoundaryBuffers(phydro->u, true);
      if (MAGNETIC_FIELDS_ENABLED)
        pmb->pbval->SendFieldBoundaryBuffers(pfield->b);
      pmb=pmb->next;
    }

    // wait to receive conserved variables
    pmb = pblock;
    while (pmb != NULL)  {
      phydro=pmb->phydro;
      pfield=pmb->pfield;
      pbval=pmb->pbval;
      pbval->ReceiveHydroBoundaryBuffersWithWait(phydro->u, true);
      if (MAGNETIC_FIELDS_ENABLED)
        pbval->ReceiveFieldBoundaryBuffersWithWait(pfield->b);
      pmb->pbval->ClearBoundaryForInit(true);
      pmb=pmb->next;
    }

    // With AMR/SMR GR send primitives to enable cons->prim before prolongation
    if (GENERAL_RELATIVITY && multilevel) {

      // prepare to receive primitives
      pmb = pblock;
      while (pmb != NULL) {
        pmb->pbval->StartReceivingForInit(false);
        pmb=pmb->next;
      }

      // send primitives
      pmb = pblock;
      while (pmb != NULL) {
        phydro=pmb->phydro;
        pmb->pbval->SendHydroBoundaryBuffers(phydro->w, false);
        pmb=pmb->next;
      }

      // wait to receive AMR/SMR GR primitives
      pmb = pblock;
      while (pmb != NULL) {
        phydro=pmb->phydro;
        pfield=pmb->pfield;
        pbval=pmb->pbval;
        pbval->ReceiveHydroBoundaryBuffersWithWait(phydro->w, false);
        pmb->pbval->ClearBoundaryForInit(false);
        pmb=pmb->next;
      }
    }

    // Now do prolongation, compute primitives, apply BCs
    pmb = pblock;
    while (pmb != NULL)  {
      phydro=pmb->phydro;
      pfield=pmb->pfield;
      pbval=pmb->pbval;
      if(multilevel==true)
        pbval->ProlongateBoundaries(phydro->w, phydro->u, pfield->b, pfield->bcc,
                                    time, 0.0);

      int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
      if(pmb->nblevel[1][1][0]!=-1) is-=NGHOST;
      if(pmb->nblevel[1][1][2]!=-1) ie+=NGHOST;
      if(pmb->block_size.nx2 > 1) {
        if(pmb->nblevel[1][0][1]!=-1) js-=NGHOST;
        if(pmb->nblevel[1][2][1]!=-1) je+=NGHOST;
      }
      if(pmb->block_size.nx3 > 1) {
        if(pmb->nblevel[0][1][1]!=-1) ks-=NGHOST;
        if(pmb->nblevel[2][1][1]!=-1) ke+=NGHOST;
      }
      pmb->peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b, 
                                      phydro->w, pfield->bcc, pmb->pcoord,
                                      is, ie, js, je, ks, ke);
      pbval->ApplyPhysicalBoundaries(phydro->w, phydro->u, pfield->b, pfield->bcc,
                                     time, 0.0);
      pmb=pmb->next;
    }

    if((res_flag==0) && (adaptive==true)) {
      iflag=false;
      int onb=nbtotal;
      pmb = pblock;
      while (pmb != NULL)  {
        pmb->pmr->CheckRefinementCondition();
        pmb=pmb->next;
      }
      AdaptiveMeshRefinement(pin);
      if(nbtotal==onb) iflag=true;
      else if(nbtotal < onb && Globals::my_rank==0) {
         std::cout << "### Warning in Mesh::Initialize" << std::endl
         << "The number of MeshBlocks decreased during AMR grid initialization." << std::endl
         << "Possibly the refinement criteria have a problem." << std::endl;
      }
      if(nbtotal > 2*inb && Globals::my_rank==0) {
         std::cout << "### Warning in Mesh::Initialize" << std::endl
         << "The number of MeshBlocks increased more than twice during initialization."<< std::endl
         << "More computing power than you expected may be required." << std::endl;
      }
    }
  } while(iflag==false);

  // calculate the first time step
  pmb = pblock;
  while (pmb != NULL)  {
    pmb->phydro->NewBlockTimeStep();
    pmb=pmb->next;
  }
  NewTimeStep();
  return;
}

//----------------------------------------------------------------------------------------
//! \fn MeshBlock* Mesh::FindMeshBlock(int tgid)
//  \brief return the MeshBlock whose gid is tgid

MeshBlock* Mesh::FindMeshBlock(int tgid)
{
  MeshBlock *pbl=pblock;
  while(pbl!=NULL)
  {
    if(pbl->gid==tgid)
      break;
    pbl=pbl->next;
  }
  return pbl;
}

//----------------------------------------------------------------------------------------
// \!fn void Mesh::LoadBalance(Real *clist, int *rlist, int *slist, int *nlist, int nb)
// \brief Calculate distribution of MeshBlocks based on the cost list

void Mesh::LoadBalance(Real *clist, int *rlist, int *slist, int *nlist, int nb)
{
  std::stringstream msg;
  Real totalcost=0, maxcost=0.0, mincost=(FLT_MAX);

  for(int i=0; i<nb; i++) {
    totalcost+=clist[i];
    mincost=std::min(mincost,clist[i]);
    maxcost=std::max(maxcost,clist[i]);
  }
  int j=(Globals::nranks)-1;
  Real targetcost=totalcost/Globals::nranks;
  Real mycost=0.0;
  // create rank list from the end: the master node should have less load
  for(int i=nb-1;i>=0;i--) {
    if(targetcost==0.0) {
      msg << "### FATAL ERROR in LoadBalance" << std::endl
          << "There is at least one process which has no MeshBlock" << std::endl
          << "Decrease the number of processes or use smaller MeshBlocks." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    mycost+=clist[i];
    rlist[i]=j;
    if(mycost >= targetcost && j>0) {
      j--;
      totalcost-=mycost;
      mycost=0.0;
      targetcost=totalcost/(j+1);
    }
  }
  slist[0]=0;
  j=0;
  for(int i=1;i<nb;i++) { // make the list of nbstart and nblocks
    if(rlist[i]!=rlist[i-1]) {
      nlist[j]=i-nslist[j];
      slist[++j]=i;
    }
  }
  nlist[j]=nb-slist[j];

#ifdef MPI_PARALLEL
  if(nb % Globals::nranks != 0 && adaptive == false
  && maxcost == mincost && Globals::my_rank==0) {
    std::cout << "### Warning in LoadBalance" << std::endl
              << "The number of MeshBlocks cannot be divided evenly. "
              << "This will cause a poor load balance." << std::endl;
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
// \!fn void Mesh::SetBlockSizeAndBoundaries(LogicalLocation loc,
//                 RegionSize &block_size, enum BundaryFlag *block_bcs)
// \brief Set the physical part of a block_size structure and block boundary conditions

void Mesh::SetBlockSizeAndBoundaries(LogicalLocation loc, RegionSize &block_size,
                                     enum BoundaryFlag *block_bcs)
{
  long int &lx1=loc.lx1;
  long int &lx2=loc.lx2;
  long int &lx3=loc.lx3;
  int &ll=loc.level;
  // calculate physical block size, x1
  if(lx1==0) {
    block_size.x1min=mesh_size.x1min;
    block_bcs[INNER_X1]=mesh_bcs[INNER_X1];
  }
  else {
    Real rx=(Real)lx1/(Real)(nrbx1<<(ll-root_level));
    block_size.x1min=MeshGenerator_[X1DIR](rx,mesh_size);
    block_bcs[INNER_X1]=BLOCK_BNDRY;
  }
  if(lx1==(nrbx1<<(ll-root_level))-1) {
    block_size.x1max=mesh_size.x1max;
    block_bcs[OUTER_X1]=mesh_bcs[OUTER_X1];
  }
  else {
    Real rx=(Real)(lx1+1)/(Real)(nrbx1<<(ll-root_level));
    block_size.x1max=MeshGenerator_[X1DIR](rx,mesh_size);
    block_bcs[OUTER_X1]=BLOCK_BNDRY;
  }

  // calculate physical block size, x2
  if(mesh_size.nx2 == 1) {
    block_size.x2min=mesh_size.x2min;
    block_size.x2max=mesh_size.x2max;
    block_bcs[INNER_X2]=mesh_bcs[INNER_X2];
    block_bcs[OUTER_X2]=mesh_bcs[OUTER_X2];
  }
  else {
    if(lx2==0) {
      block_size.x2min=mesh_size.x2min;
      block_bcs[INNER_X2]=mesh_bcs[INNER_X2];
    }
    else {
      Real rx=(Real)lx2/(Real)(nrbx2<<(ll-root_level));
      block_size.x2min=MeshGenerator_[X2DIR](rx,mesh_size);
      block_bcs[INNER_X2]=BLOCK_BNDRY;
    }
    if(lx2==(nrbx2<<(ll-root_level))-1) {
      block_size.x2max=mesh_size.x2max;
      block_bcs[OUTER_X2]=mesh_bcs[OUTER_X2];
    }
    else {
      Real rx=(Real)(lx2+1)/(Real)(nrbx2<<(ll-root_level));
      block_size.x2max=MeshGenerator_[X2DIR](rx,mesh_size);
      block_bcs[OUTER_X2]=BLOCK_BNDRY;
    }
  }

  // calculate physical block size, x3
  if(mesh_size.nx3 == 1) {
    block_size.x3min=mesh_size.x3min;
    block_size.x3max=mesh_size.x3max;
    block_bcs[INNER_X3]=mesh_bcs[INNER_X3];
    block_bcs[OUTER_X3]=mesh_bcs[OUTER_X3];
  }
  else {
    if(lx3==0) {
      block_size.x3min=mesh_size.x3min;
      block_bcs[INNER_X3]=mesh_bcs[INNER_X3];
    }
    else {
      Real rx=(Real)lx3/(Real)(nrbx3<<(ll-root_level));
      block_size.x3min=MeshGenerator_[X3DIR](rx,mesh_size);
      block_bcs[INNER_X3]=BLOCK_BNDRY;
    }
    if(lx3==(nrbx3<<(ll-root_level))-1) {
      block_size.x3max=mesh_size.x3max;
      block_bcs[OUTER_X3]=mesh_bcs[OUTER_X3];
    }
    else {
      Real rx=(Real)(lx3+1)/(Real)(nrbx3<<(ll-root_level));
      block_size.x3max=MeshGenerator_[X3DIR](rx,mesh_size);
      block_bcs[OUTER_X3]=BLOCK_BNDRY;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// \!fn void Mesh::AdaptiveMeshRefinement(ParameterInput *pin)
// \brief Main function for adaptive mesh refinement

void Mesh::AdaptiveMeshRefinement(ParameterInput *pin)
{
  MeshBlock *pmb;

  // collect refinement flags from all the meshblocks
  // count the number of the blocks to be (de)refined
  nref[Globals::my_rank]=0;
  nderef[Globals::my_rank]=0;
  pmb=pblock;
  while(pmb!=NULL) {
    if(pmb->pmr->refine_flag_== 1) nref[Globals::my_rank]++;
    if(pmb->pmr->refine_flag_==-1) nderef[Globals::my_rank]++;
    pmb=pmb->next;
  }
#ifdef MPI_PARALLEL
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, nref,   1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, nderef, 1, MPI_INT, MPI_COMM_WORLD);
#endif

  // count the number of the blocks to be (de)refined and displacement
  int tnref=0, tnderef=0;
  for(int n=0; n<Globals::nranks; n++) {
    tnref  += nref[n];
    tnderef+= nderef[n];
  }
  if(tnref==0 && tnderef==0) // nothing to do
    return;

  int rd=0, dd=0;
  for(int n=0; n<Globals::nranks; n++) {
    bnref[n]   = nref[n]*sizeof(LogicalLocation);
    bnderef[n] = nderef[n]*sizeof(LogicalLocation);
    rdisp[n] = rd;
    ddisp[n] = dd;
    brdisp[n] = rd*sizeof(LogicalLocation);
    bddisp[n] = dd*sizeof(LogicalLocation);
    rd+=nref[n];
    dd+=nderef[n];
  }

  // allocate memory for the location arrays
  int nlbl=2, dim=1;
  if(mesh_size.nx2 > 1) nlbl=4, dim=2;
  if(mesh_size.nx3 > 1) nlbl=8, dim=3;
  LogicalLocation *lref, *lderef, *clderef;
  if(tnref!=0)
    lref = new LogicalLocation[tnref];
  if(tnderef>nlbl) {
    lderef = new LogicalLocation[tnderef];
    clderef = new LogicalLocation[tnderef/nlbl];
  }

  // collect the locations and costs
  int iref = rdisp[Globals::my_rank], ideref = ddisp[Globals::my_rank];
  pmb=pblock;
  while(pmb!=NULL) {
    if(pmb->pmr->refine_flag_== 1)
      lref[iref++]=pmb->loc;
    if(pmb->pmr->refine_flag_==-1 && tnderef>nlbl)
      lderef[ideref++]=pmb->loc;
    pmb=pmb->next;
  }
#ifdef MPI_PARALLEL
  if(tnref>0 && tnderef>nlbl) {
    MPI_Allgatherv(MPI_IN_PLACE, bnref[Globals::my_rank],   MPI_BYTE,
                   lref,   bnref,   brdisp, MPI_BYTE, MPI_COMM_WORLD);
    MPI_Allgatherv(MPI_IN_PLACE, bnderef[Globals::my_rank], MPI_BYTE,
                   lderef, bnderef, bddisp, MPI_BYTE, MPI_COMM_WORLD);
  }
  else if(tnref>0) {
    MPI_Allgatherv(MPI_IN_PLACE, bnref[Globals::my_rank],   MPI_BYTE,
                   lref,   bnref,   brdisp, MPI_BYTE, MPI_COMM_WORLD);
  }
  else if(tnderef>nlbl) {
    MPI_Allgatherv(MPI_IN_PLACE, bnderef[Globals::my_rank], MPI_BYTE,
                   lderef, bnderef, bddisp, MPI_BYTE, MPI_COMM_WORLD);
  }
#endif

  // calculate the list of the newly derefined blocks
  int ctnd=0;
  if(tnderef>nlbl) {
    int lk=0, lj=0;
    if(mesh_size.nx2 > 1) lj=1;
    if(mesh_size.nx3 > 1) lk=1;
    for(int n=0; n<tnderef; n++) {
      if((lderef[n].lx1&1L)==0 && (lderef[n].lx2&1L)==0 && (lderef[n].lx3&1L)==0) {
        int r=n, rr=0;
        for(long int k=0;k<=lk;k++) {
          for(long int j=0;j<=lj;j++) {
            for(long int i=0;i<=1;i++) {
              if((lderef[n].lx1+i)==lderef[r].lx1
              && (lderef[n].lx2+j)==lderef[r].lx2
              && (lderef[n].lx3+k)==lderef[r].lx3
              &&  lderef[n].level ==lderef[r].level)
                rr++;
              r++;
            }
          }
        }
        if(rr==nlbl) {
          clderef[ctnd].lx1  =(lderef[n].lx1>>1);
          clderef[ctnd].lx2  =(lderef[n].lx2>>1);
          clderef[ctnd].lx3  =(lderef[n].lx3>>1);
          clderef[ctnd].level=lderef[n].level-1;
          ctnd++;
        }
      }
    }
  }
  // sort the lists by level
  if(ctnd>1)
    std::sort(clderef, &(clderef[ctnd-1]), LogicalLocation::Greater);

  if(tnderef>nlbl)
    delete [] lderef;

  // Now the lists of the blocks to be refined and derefined are completed
  // Start tree manipulation
  // Step 1. perform refinement
  int nnew=0, ndel=0, ntot=0;
  for(int n=0; n<tnref; n++) {
    MeshBlockTree *bt=tree.FindMeshBlock(lref[n]);
    bt->Refine(tree, dim, mesh_bcs, nrbx1, nrbx2, nrbx3, root_level, nnew);
  }
  if(tnref!=0)
    delete [] lref;
  // Step 2. perform derefinement
  for(int n=0; n<ctnd; n++) {
    MeshBlockTree *bt=tree.FindMeshBlock(clderef[n]);
    bt->Derefine(tree, dim, mesh_bcs, nrbx1, nrbx2, nrbx3, root_level, ndel);
  }
  if(tnderef>nlbl)
    delete [] clderef;
  ntot=nbtotal+nnew-ndel;
  if(nnew==0 && ndel==0)
    return; // nothing to do
  // Tree manipulation completed
  nbnew+=nnew; nbdel+=ndel;

  // Block exchange
  // Step 1. construct new lists
  LogicalLocation *newloc = new LogicalLocation[ntot];
  int *newrank = new int[ntot];
  Real *newcost = new Real[ntot];
  int *newtoold = new int[ntot];
  int *oldtonew = new int[nbtotal];
  int nbtold=nbtotal;
  tree.GetMeshBlockList(newloc,newtoold,nbtotal);
  // create a list mapping the previous gid to the current one

  oldtonew[0]=0;
  int k=1;
  for(int n=1; n<ntot; n++) {
    if(newtoold[n]==newtoold[n-1]+1) { // normal
      oldtonew[k++]=n;
    }
    else if(newtoold[n]==newtoold[n-1]+nlbl) { // derefined
      for(int j=0; j<nlbl-1; j++)
        oldtonew[k++]=n-1;
      oldtonew[k++]=n;
    }
  }
  // fill the last block 
  for(;k<nbtold; k++)
    oldtonew[k]=ntot-1;

#ifdef MPI_PARALLEL
  // share the cost list
  MPI_Allgatherv(MPI_IN_PLACE, nblist[Globals::my_rank], MPI_INT,
                 costlist, nblist, nslist, MPI_INT, MPI_COMM_WORLD);
#endif

  current_level=0;
  for(int n=0; n<ntot; n++) {
    int on=newtoold[n];
    if(newloc[n].level>current_level) // set the current max level
      current_level=newloc[n].level;
    if(newloc[n].level>=loclist[on].level) // same or refined
      newcost[n]=costlist[on];
    else {
      Real acost=0.0;
      for(int l=0; l<nlbl; l++)
        acost+=costlist[on+l];
      newcost[n]=acost/nlbl;
    }
  }

  // store old nbstart and nbend
  int onbs=nslist[Globals::my_rank];
  int onbe=onbs+nblist[Globals::my_rank]-1;

  // Step 2. Calculate new load balance
  LoadBalance(newcost, newrank, nslist, nblist, ntot);

  int nbs=nslist[Globals::my_rank];
  int nbe=nbs+nblist[Globals::my_rank]-1;

  int f2, f3;
  int &bnx1=pblock->block_size.nx1;
  int &bnx2=pblock->block_size.nx2;
  int &bnx3=pblock->block_size.nx3;
  if(mesh_size.nx2>1) f2=1;
  else f2=0;
  if(mesh_size.nx3>1) f3=1;
  else f3=0;

#ifdef MPI_PARALLEL
  // Step 3. count the number of the blocks to be sent / received
  int nsend=0, nrecv=0;
  for(int n=nbs; n<=nbe; n++) { 
    int on=newtoold[n];
    if(loclist[on].level > newloc[n].level) { // f2c
      for(int k=0; k<nlbl; k++) {
        if(ranklist[on+k]!=Globals::my_rank)
          nrecv++;
      }
    }
    else {
      if(ranklist[on]!=Globals::my_rank)
        nrecv++;
    }
  }
  for(int n=onbs; n<=onbe; n++) { 
    int nn=oldtonew[n];
    if(loclist[n].level < newloc[nn].level) { // c2f
      for(int k=0; k<nlbl; k++) {
        if(newrank[nn+k]!=Globals::my_rank)
          nsend++;
      }
    }
    else {
      if(newrank[nn]!=Globals::my_rank)
        nsend++;
    }
  }

  // Step 4. calculate buffer sizes
  Real **sendbuf, **recvbuf;
  int bssame=bnx1*bnx2*bnx3*NHYDRO;
  int bsf2c=(bnx1/2)*((bnx2+1)/2)*((bnx3+1)/2)*NHYDRO;
  int bsc2f=(bnx1/2+2)*((bnx2+1)/2+2*f2)*((bnx3+1)/2+2*f3)*NHYDRO;
  if(MAGNETIC_FIELDS_ENABLED) {
    bssame+=(bnx1+1)*bnx2*bnx3+bnx1*(bnx2+f2)*bnx3+bnx1*bnx2*(bnx3+f3);
    bsf2c+=((bnx1/2)+1)*((bnx2+1)/2)*((bnx3+1)/2)
          +(bnx1/2)*(((bnx2+1)/2)+f2)*((bnx3+1)/2)
          +(bnx1/2)*((bnx2+1)/2)*(((bnx3+1)/2)+f3);
    bsc2f+=((bnx1/2)+1+2)*((bnx2+1)/2+2*f2)*((bnx3+1)/2+2*f3)
          +(bnx1/2+2)*(((bnx2+1)/2)+f2+2*f2)*((bnx3+1)/2+2*f3)
          +(bnx1/2+2)*((bnx2+1)/2+2*f2)*(((bnx3+1)/2)+f3+2*f3);
  }

  MPI_Request *req_send, *req_recv;
  // Step 5. allocate and start receiving buffers
  if(nrecv!=0) {
    recvbuf = new Real*[nrecv];
    req_recv = new MPI_Request[nrecv];
    int k=0;
    for(int n=nbs; n<=nbe; n++) { 
      int on=newtoold[n];
      LogicalLocation &oloc=loclist[on];
      LogicalLocation &nloc=newloc[n];
      if(oloc.level>nloc.level) { // f2c
        for(int l=0; l<nlbl; l++) {
          if(ranklist[on+l]==Globals::my_rank) continue;
          LogicalLocation &lloc=loclist[on+l];
          int ox1=lloc.lx1&1L, ox2=lloc.lx2&1L, ox3=lloc.lx3&1L;
          recvbuf[k] = new Real[bsf2c];
          int tag=CreateAMRMPITag(n-nbs, ox1, ox2, ox3);
          MPI_Irecv(recvbuf[k], bsf2c, MPI_ATHENA_REAL, ranklist[on+l],
                    tag, MPI_COMM_WORLD, &(req_recv[k]));
          k++;
        }
      }
      else { // same or c2f
        if(ranklist[on]==Globals::my_rank) continue;
        int size;
        if(oloc.level == nloc.level) size=bssame;
        else size=bsc2f;
        recvbuf[k] = new Real[size];
        int tag=CreateAMRMPITag(n-nbs, 0, 0, 0);
        MPI_Irecv(recvbuf[k], size, MPI_ATHENA_REAL, ranklist[on],
                  tag, MPI_COMM_WORLD, &(req_recv[k]));
        k++;
      }
    }
  }
  // Step 6. allocate, pack and start sending buffers
  if(nsend!=0) {
    sendbuf = new Real*[nsend];
    req_send = new MPI_Request[nsend];
    int k=0;
    for(int n=onbs; n<=onbe; n++) { 
      int nn=oldtonew[n];
      LogicalLocation &oloc=loclist[n];
      LogicalLocation &nloc=newloc[nn];
      MeshBlock* pb=FindMeshBlock(n);
      if(nloc.level==oloc.level) { // same
        if(newrank[nn]==Globals::my_rank) continue;
        sendbuf[k] = new Real[bssame];
        // pack
        int p=0;
        BufferUtility::Pack4DData(pb->phydro->u, sendbuf[k], 0, NHYDRO-1,
                       pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke, p);
        if(MAGNETIC_FIELDS_ENABLED) {
          BufferUtility::Pack3DData(pb->pfield->b.x1f, sendbuf[k],
                         pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
          BufferUtility::Pack3DData(pb->pfield->b.x2f, sendbuf[k],
                         pb->is, pb->ie, pb->js, pb->je+f2, pb->ks, pb->ke, p);
          BufferUtility::Pack3DData(pb->pfield->b.x3f, sendbuf[k],
                         pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke+f3, p);
        }
        int tag=CreateAMRMPITag(nn-nslist[newrank[nn]], 0, 0, 0);
        MPI_Isend(sendbuf[k], bssame, MPI_ATHENA_REAL, newrank[nn],
                  tag, MPI_COMM_WORLD, &(req_send[k]));
        k++;
      }
      else if(nloc.level>oloc.level) { // c2f
        for(int l=0; l<nlbl; l++) {
          if(newrank[nn+l]==Globals::my_rank) continue;
          LogicalLocation &lloc=newloc[nn+l];
          int ox1=lloc.lx1&1L, ox2=lloc.lx2&1L, ox3=lloc.lx3&1L;
          sendbuf[k] = new Real[bsc2f];
          // pack
          int is, ie, js, je, ks, ke;
          if(ox1==0) is=pb->is-1,                       ie=pb->is+pb->block_size.nx1/2;
          else       is=pb->is+pb->block_size.nx1/2-1,  ie=pb->ie+1;
          if(ox2==0) js=pb->js-f2,                      je=pb->js+pb->block_size.nx2/2;
          else       js=pb->js+pb->block_size.nx2/2-f2, je=pb->je+f2;
          if(ox3==0) ks=pb->ks-f3,                      ke=pb->ks+pb->block_size.nx3/2;
          else       ks=pb->ks+pb->block_size.nx3/2-f3, ke=pb->ke+f3;
          int p=0;
          BufferUtility::Pack4DData(pb->phydro->u, sendbuf[k], 0, NHYDRO-1,
                                    is, ie, js, je, ks, ke, p);
          if(MAGNETIC_FIELDS_ENABLED) {
            BufferUtility::Pack3DData(pb->pfield->b.x1f, sendbuf[k],
                                      is, ie+1, js, je, ks, ke, p);
            BufferUtility::Pack3DData(pb->pfield->b.x2f, sendbuf[k],
                                      is, ie, js, je+f2, ks, ke, p);
            BufferUtility::Pack3DData(pb->pfield->b.x3f, sendbuf[k],
                                      is, ie, js, je, ks, ke+f3, p);
          }
          int tag=CreateAMRMPITag(nn+l-nslist[newrank[nn+l]], 0, 0, 0);
          MPI_Isend(sendbuf[k], bsc2f, MPI_ATHENA_REAL, newrank[nn+l],
                    tag, MPI_COMM_WORLD, &(req_send[k]));
          k++;
        }
      }
      else { // f2c
        if(newrank[nn]==Globals::my_rank) continue;
        int ox1=oloc.lx1&1L, ox2=oloc.lx2&1L, ox3=oloc.lx3&1L;
        sendbuf[k] = new Real[bsf2c];
        // restrict and pack
        MeshRefinement *pmr=pb->pmr;
        pmr->RestrictCellCenteredValues(pb->phydro->u, pmr->coarse_cons_,
             0, NHYDRO-1, pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke);
        int p=0;
        BufferUtility::Pack4DData(pmr->coarse_cons_, sendbuf[k], 0, NHYDRO-1,
                       pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke, p);
        if(MAGNETIC_FIELDS_ENABLED) {
          pmr->RestrictFieldX1(pb->pfield->b.x1f, pmr->coarse_b_.x1f,
                               pb->cis, pb->cie+1, pb->cjs, pb->cje, pb->cks, pb->cke);
          BufferUtility::Pack3DData(pmr->coarse_b_.x1f, sendbuf[k],
                         pb->cis, pb->cie+1, pb->cjs, pb->cje, pb->cks, pb->cke, p);
          pmr->RestrictFieldX2(pb->pfield->b.x2f, pmr->coarse_b_.x2f,
                               pb->cis, pb->cie, pb->cjs, pb->cje+f2, pb->cks, pb->cke);
          BufferUtility::Pack3DData(pmr->coarse_b_.x2f, sendbuf[k],
                         pb->cis, pb->cie, pb->cjs, pb->cje+f2, pb->cks, pb->cke, p);
          pmr->RestrictFieldX3(pb->pfield->b.x3f, pmr->coarse_b_.x3f,
                               pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke+f3);
          BufferUtility::Pack3DData(pmr->coarse_b_.x3f, sendbuf[k],
                         pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke+f3, p);
        }
        int tag=CreateAMRMPITag(nn-nslist[newrank[nn]], ox1, ox2, ox3);
        MPI_Isend(sendbuf[k], bsf2c, MPI_ATHENA_REAL, newrank[nn],
                  tag, MPI_COMM_WORLD, &(req_send[k]));
        k++;
      }
    }
  }
#endif

  // Step 7. construct a new MeshBlock list
  // move the data within the node
  MeshBlock *newlist=NULL;
  RegionSize block_size=pblock->block_size;
  enum BoundaryFlag block_bcs[6];

  for(int n=nbs; n<=nbe; n++) {
    int on=newtoold[n];
    if((ranklist[on]==Globals::my_rank) && (loclist[on].level == newloc[n].level)) {
      // on the same node and same level -> just move it
      MeshBlock* pob=FindMeshBlock(on);
      if(pob->prev==NULL) pblock=pob->next;
      else pob->prev->next=pob->next;
      if(pob->next!=NULL) pob->next->prev=pob->prev;
      pob->next=NULL;
      if(n==nbs) { // first
        pob->prev=NULL;
        newlist=pmb=pob;
      }
      else {
        pmb->next=pob;
        pob->prev=pmb;
        pmb=pmb->next;
      }
      pmb->gid=n; pmb->lid=n-nbs;
    }
    else {
      // on a different level or node - create a new block
      SetBlockSizeAndBoundaries(newloc[n], block_size, block_bcs);
      if(n==nbs) { // first
        newlist = new MeshBlock(n, n-nbs, newloc[n], block_size, block_bcs, this, pin, true);
        pmb=newlist;
      }
      else {
        pmb->next = new MeshBlock(n, n-nbs, newloc[n], block_size, block_bcs, this, pin, true);
        pmb->next->prev=pmb;
        pmb=pmb->next;
      }
      // fill the conservative variables
      if((loclist[on].level>newloc[n].level)) { // fine to coarse
        for(int ll=0; ll<nlbl; ll++) {
          if(ranklist[on+ll]!=Globals::my_rank) continue;
          // on the same node - restriction
          MeshBlock* pob=FindMeshBlock(on+ll);
          MeshRefinement *pmr=pob->pmr;
          pmr->RestrictCellCenteredValues(pob->phydro->u, pmr->coarse_cons_,
               0, NHYDRO-1, pob->cis, pob->cie, pob->cjs, pob->cje, pob->cks, pob->cke);
          int is=pmb->is+(loclist[on+ll].lx1&1L)*pmb->block_size.nx1/2;
          int js=pmb->js+(loclist[on+ll].lx2&1L)*pmb->block_size.nx2/2;
          int ks=pmb->ks+(loclist[on+ll].lx3&1L)*pmb->block_size.nx3/2;
          AthenaArray<Real> &src=pmr->coarse_cons_;
          AthenaArray<Real> &dst=pmb->phydro->u;
          for(int nv=0; nv<NHYDRO; nv++) {
            for(int k=ks, fk=pob->cks; fk<=pob->cke; k++, fk++) {
              for(int j=js, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
                for(int i=is, fi=pob->cis; fi<=pob->cie; i++, fi++)
                  dst(nv, k, j, i)=src(nv, fk, fj, fi);
          }}}
          if(MAGNETIC_FIELDS_ENABLED) {
            pmr->RestrictFieldX1(pob->pfield->b.x1f, pmr->coarse_b_.x1f,
                         pob->cis, pob->cie+1, pob->cjs, pob->cje, pob->cks, pob->cke);
            pmr->RestrictFieldX2(pob->pfield->b.x2f, pmr->coarse_b_.x2f,
                         pob->cis, pob->cie, pob->cjs, pob->cje+f2, pob->cks, pob->cke);
            pmr->RestrictFieldX3(pob->pfield->b.x3f, pmr->coarse_b_.x3f,
                         pob->cis, pob->cie, pob->cjs, pob->cje, pob->cks, pob->cke+f3);
            FaceField &src=pmr->coarse_b_;
            FaceField &dst=pmb->pfield->b;
            for(int k=ks, fk=pob->cks; fk<=pob->cke; k++, fk++) {
              for(int j=js, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
                for(int i=is, fi=pob->cis; fi<=pob->cie+1; i++, fi++)
                  dst.x1f(k, j, i)=src.x1f(fk, fj, fi);
            }}
            for(int k=ks, fk=pob->cks; fk<=pob->cke; k++, fk++) {
              for(int j=js, fj=pob->cjs; fj<=pob->cje+f2; j++, fj++) {
                for(int i=is, fi=pob->cis; fi<=pob->cie; i++, fi++)
                  dst.x2f(k, j, i)=src.x2f(fk, fj, fi);
            }}
            if(pmb->block_size.nx2==1) {
              int ie=is+block_size.nx1/2-1;
              for(int i=is; i<=ie; i++)
                dst.x2f(pmb->ks, pmb->js+1, i)=dst.x2f(pmb->ks, pmb->js, i);
            }
            for(int k=ks, fk=pob->cks; fk<=pob->cke+f3; k++, fk++) {
              for(int j=js, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
                for(int i=is, fi=pob->cis; fi<=pob->cie; i++, fi++)
                  dst.x3f(k, j, i)=src.x3f(fk, fj, fi);
            }}
            if(pmb->block_size.nx3==1) {
              int ie=is+block_size.nx1/2-1, je=js+block_size.nx2/2-1;
              for(int j=js; j<=je; j++) {
                for(int i=is; i<=ie; i++)
                  dst.x3f(pmb->ks+1, j, i)=dst.x3f(pmb->ks, j, i);
              }
            }
          }
        }
      }
      else if((loclist[on].level < newloc[n].level) && (ranklist[on]==Globals::my_rank)) {
        // coarse to fine on the same node - prolongation
        if(ranklist[on]!=Globals::my_rank) continue;
        MeshBlock* pob=FindMeshBlock(on);
        MeshRefinement *pmr=pmb->pmr;
        int is=pob->cis-1, ie=pob->cie+1, js=pob->cjs-f2,
            je=pob->cje+f2, ks=pob->cks-f3, ke=pob->cke+f3;
        int cis=(newloc[n].lx1&1L)*pob->block_size.nx1/2+pob->is-1;
        int cjs=(newloc[n].lx2&1L)*pob->block_size.nx2/2+pob->js-f2;
        int cks=(newloc[n].lx3&1L)*pob->block_size.nx3/2+pob->ks-f3;
        AthenaArray<Real> &src=pob->phydro->u;
        AthenaArray<Real> &dst=pmr->coarse_cons_;
        // fill the coarse buffer
        for(int nv=0; nv<NHYDRO; nv++) {
          for(int k=ks, ck=cks; k<=ke; k++, ck++) {
            for(int j=js, cj=cjs; j<=je; j++, cj++) {
              for(int i=is, ci=cis; i<=ie; i++, ci++)
                dst(nv, k, j, i)=src(nv, ck, cj, ci);
        }}}
        pmr->ProlongateCellCenteredValues(dst, pmb->phydro->u, 0, NHYDRO-1,
                                          is, ie, js, je, ks, ke);
        if(MAGNETIC_FIELDS_ENABLED) {
          FaceField &src=pob->pfield->b;
          FaceField &dst=pmr->coarse_b_;
          for(int k=ks, ck=cks; k<=ke; k++, ck++) {
            for(int j=js, cj=cjs; j<=je; j++, cj++) {
              for(int i=is, ci=cis; i<=ie+1; i++, ci++)
                dst.x1f(k, j, i)=src.x1f(ck, cj, ci);
          }}
          for(int k=ks, ck=cks; k<=ke; k++, ck++) {
            for(int j=js, cj=cjs; j<=je+f2; j++, cj++) {
              for(int i=is, ci=cis; i<=ie; i++, ci++)
                dst.x2f(k, j, i)=src.x2f(ck, cj, ci);
          }}
          for(int k=ks, ck=cks; k<=ke+f3; k++, ck++) {
            for(int j=js, cj=cjs; j<=je; j++, cj++) {
              for(int i=is, ci=cis; i<=ie; i++, ci++)
                dst.x3f(k, j, i)=src.x3f(ck, cj, ci);
          }}
          pmr->ProlongateSharedFieldX1(dst.x1f, pmb->pfield->b.x1f,
                                       pob->is, ie+1, js, je, ks, ke);
          pmr->ProlongateSharedFieldX2(dst.x2f, pmb->pfield->b.x2f,
                                       is, ie, js, je+f2, ks, ke);
          pmr->ProlongateSharedFieldX3(dst.x3f, pmb->pfield->b.x3f,
                                       is, ie, js, je, ks, ke+f3);
          pmr->ProlongateInternalField(pmb->pfield->b, is, ie, js, je, ks, ke);
        }
      }
    }
  }

  // discard remaining MeshBlocks
  // they could be reused, but for the moment, just throw them away for simplicity
  if(pblock!=NULL) {
    while(pblock->next != NULL)
      delete pblock->next;
    delete pblock;
  }

  // Replace the MeshBlock list
  pblock=newlist;

  // Step 8. Receive the data and load into MeshBlocks
  // This is a test: try MPI_Waitall later.
#ifdef MPI_PARALLEL
  if(nrecv!=0) {
    int k=0;
    for(int n=nbs; n<=nbe; n++) { 
      int on=newtoold[n];
      LogicalLocation &oloc=loclist[on];
      LogicalLocation &nloc=newloc[n];
      MeshBlock *pb=FindMeshBlock(n);
      if(oloc.level==nloc.level) { // same
        if(ranklist[on]==Globals::my_rank) continue;
        MPI_Wait(&(req_recv[k]), MPI_STATUS_IGNORE);
        int p=0;
        BufferUtility::Unpack4DData(recvbuf[k], pb->phydro->u, 0, NHYDRO-1,
                       pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke, p);
        if(MAGNETIC_FIELDS_ENABLED) {
          FaceField &dst=pb->pfield->b;
          BufferUtility::Unpack3DData(recvbuf[k], dst.x1f,
                         pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
          BufferUtility::Unpack3DData(recvbuf[k], dst.x2f,
                         pb->is, pb->ie, pb->js, pb->je+f2, pb->ks, pb->ke, p);
          BufferUtility::Unpack3DData(recvbuf[k], dst.x3f,
                         pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke+f3, p);
          if(pb->block_size.nx2==1) {
            for(int i=pb->is; i<=pb->ie; i++)
              dst.x2f(pb->ks, pb->js+1, i)=dst.x2f(pb->ks, pb->js, i);
          }
          if(pb->block_size.nx3==1) {
            for(int j=pb->js; j<=pb->je; j++) {
              for(int i=pb->is; i<=pb->ie; i++)
                dst.x3f(pb->ks+1, j, i)=dst.x3f(pb->ks, j, i);
            }
          }
        }
        k++;
      }
      else if(oloc.level>nloc.level) { // f2c
        for(int l=0; l<nlbl; l++) {
          if(ranklist[on+l]==Globals::my_rank) continue;
          LogicalLocation &lloc=loclist[on+l];
          int ox1=lloc.lx1&1L, ox2=lloc.lx2&1L, ox3=lloc.lx3&1L;
          int p=0, is, ie, js, je, ks, ke;
          if(ox1==0) is=pb->is,                      ie=pb->is+pb->block_size.nx1/2-1;
          else       is=pb->is+pb->block_size.nx1/2, ie=pb->ie;
          if(ox2==0) js=pb->js,                      je=pb->js+pb->block_size.nx2/2-f2;
          else       js=pb->js+pb->block_size.nx2/2, je=pb->je;
          if(ox3==0) ks=pb->ks,                      ke=pb->ks+pb->block_size.nx3/2-f3;
          else       ks=pb->ks+pb->block_size.nx3/2, ke=pb->ke;
          MPI_Wait(&(req_recv[k]), MPI_STATUS_IGNORE);
          BufferUtility::Unpack4DData(recvbuf[k], pb->phydro->u, 0, NHYDRO-1,
                         is, ie, js, je, ks, ke, p);
          if(MAGNETIC_FIELDS_ENABLED) {
            FaceField &dst=pb->pfield->b;
            BufferUtility::Unpack3DData(recvbuf[k], dst.x1f,
                           is, ie+1, js, je, ks, ke, p);
            BufferUtility::Unpack3DData(recvbuf[k], dst.x2f,
                           is, ie, js, je+f2, ks, ke, p);
            BufferUtility::Unpack3DData(recvbuf[k], dst.x3f,
                           is, ie, js, je, ks, ke+f3, p);
            if(pb->block_size.nx2==1) {
              for(int i=is; i<=ie; i++)
                dst.x2f(pb->ks, pb->js+1, i)=dst.x2f(pb->ks, pb->js, i);
            }
            if(pb->block_size.nx3==1) {
              for(int j=js; j<=je; j++) {
                for(int i=is; i<=ie; i++)
                  dst.x3f(pb->ks+1, j, i)=dst.x3f(pb->ks, j, i);
              }
            }
          }
          k++;
        }
      }
      else { // c2f
        if(ranklist[on]==Globals::my_rank) continue;
        MeshRefinement *pmr=pb->pmr;
        int p=0;
        int is=pb->cis-1, ie=pb->cie+1, js=pb->cjs-f2,
            je=pb->cje+f2, ks=pb->cks-f3, ke=pb->cke+f3;
        MPI_Wait(&(req_recv[k]), MPI_STATUS_IGNORE);
        BufferUtility::Unpack4DData(recvbuf[k], pmr->coarse_cons_,
                                    0, NHYDRO-1, is, ie, js, je, ks, ke, p);
        pmr->ProlongateCellCenteredValues(pmr->coarse_cons_, pb->phydro->u, 0, NHYDRO-1,
                                          is, ie, js, je, ks, ke);
        if(MAGNETIC_FIELDS_ENABLED) {
          BufferUtility::Unpack3DData(recvbuf[k], pmr->coarse_b_.x1f,
                                      is, ie+1, js, je, ks, ke, p);
          BufferUtility::Unpack3DData(recvbuf[k], pmr->coarse_b_.x2f,
                                      is, ie, js, je+f2, ks, ke, p);
          BufferUtility::Unpack3DData(recvbuf[k], pmr->coarse_b_.x3f,
                                      is, ie, js, je, ks, ke+f3, p);
          pmr->ProlongateSharedFieldX1(pmr->coarse_b_.x1f, pb->pfield->b.x1f,
                                       is, ie+1, js, je, ks, ke);
          pmr->ProlongateSharedFieldX2(pmr->coarse_b_.x2f, pb->pfield->b.x2f,
                                       is, ie, js, je+f2, ks, ke);
          pmr->ProlongateSharedFieldX3(pmr->coarse_b_.x3f, pb->pfield->b.x3f,
                                       is, ie, js, je, ks, ke+f3);
          pmr->ProlongateInternalField(pb->pfield->b, is, ie, js, je, ks, ke);
        }
        k++;
      }
    }
  }
#endif

  // deallocate arrays
  delete [] loclist;
  delete [] ranklist;
  delete [] costlist;
  delete [] newtoold;
  delete [] oldtonew;
#ifdef MPI_PARALLEL
  if(nsend!=0) {
    MPI_Waitall(nsend, req_send, MPI_STATUSES_IGNORE);
    for(int n=0;n<nsend;n++)
      delete [] sendbuf[n];
    delete [] sendbuf;
    delete [] req_send;
  }
  if(nrecv!=0) {
    for(int n=0;n<nrecv;n++)
      delete [] recvbuf[n];
    delete [] recvbuf;
    delete [] req_recv;
  }
#endif

  // update the lists
  loclist = newloc;
  ranklist = newrank;
  costlist = newcost;

  // re-initialize the MeshBlocks
  pmb=pblock;
  while(pmb!=NULL) {
    pmb->SearchAndSetNeighbors(tree, ranklist, nslist);
    pmb=pmb->next;
  }
  Initialize(2, pin);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn unsigned int CreateAMRMPITag(int lid, int ox1, int ox2, int ox3)
//  \brief calculate an MPI tag for AMR block transfer
// tag = local id of destination (23) + ox1(1) + ox2(1) + ox3(1) + physics(5)

unsigned int Mesh::CreateAMRMPITag(int lid, int ox1, int ox2, int ox3)
{
  return (lid<<8) | (ox1<<7)| (ox2<<6) | (ox3<<5) | TAG_AMR;
}
