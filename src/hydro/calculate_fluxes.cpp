//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_fluxes.cpp
//  \brief Calculate hydro/MHD fluxes

// C/C++ headers
#include <algorithm>   // min,max

// Athena++ headers
#include "hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../bvals/bvals.hpp"
#include "../reconstruct/reconstruction.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes
//  \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void Hydro::CalculateFluxes(AthenaArray<Real> &w, FaceField &b,
                            AthenaArray<Real> &bcc, int step)
{
  MeshBlock *pmb=pmy_block;
  AthenaArray<Real> &x1flux=flux[X1DIR];
  AthenaArray<Real> &x2flux=flux[X2DIR];
  AthenaArray<Real> &x3flux=flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  AthenaArray<Real> b1,b2,b3,ei_x1f,ei_x2f,ei_x3f,w_x1f,w_x2f,w_x3f;
  if (MAGNETIC_FIELDS_ENABLED) {
    b1.InitWithShallowCopy(b.x1f);
    b2.InitWithShallowCopy(b.x2f);
    b3.InitWithShallowCopy(b.x3f);
    ei_x1f.InitWithShallowCopy(pmb->pfield->ei.x1f);
    ei_x2f.InitWithShallowCopy(pmb->pfield->ei.x2f);
    ei_x3f.InitWithShallowCopy(pmb->pfield->ei.x3f);
    w_x1f.InitWithShallowCopy(pmb->pfield->wght.x1f);
    w_x2f.InitWithShallowCopy(pmb->pfield->wght.x2f);
    w_x3f.InitWithShallowCopy(pmb->pfield->wght.x3f);
  }

  int tid=0;
  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif

  AthenaArray<Real> wl, wr, flx, dxw;
  wl.InitWithShallowSlice(wl_,3,tid,1);
  wr.InitWithShallowSlice(wr_,3,tid,1);
  flx.InitWithShallowSlice(flx_,3,tid,1);
  dxw.InitWithShallowSlice(dxw_,2,tid,1);

//----------------------------------------------------------------------------------------
// i-direction
  // set the loop limits
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb->block_size.nx2 > 1) {
      if(pmb->block_size.nx3 == 1) // 2D
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }
  for (int k=kl; k<=ku; ++k){ 
#pragma omp for schedule(static)
    for (int j=jl; j<=ju; ++j){

      // reconstruct L/R states
      if ((step == 1) && HYDRO_TIME_INTEGRATOR == "vl2") {
        pmb->precon->DonorCellX1(k,j,is,ie+1,w,bcc,wl,wr);
      } else {
        pmb->precon->HighResFuncX1(k,j,is,ie+1,w,bcc,wl,wr);
      }

      // compute fluxes
      RiemannSolver(k,j,is,ie+1,IVX,b1,wl,wr,flx);

      // store fluxes
      if(k>=ks && k<=ke && j>=js && j<=je) {
        for(int n=0; n<NHYDRO; n++) {
#pragma simd
          for(int i=is; i<=ie+1; i++)
            x1flux(n,k,j,i)=flx(n,i);
        }
      }

      // store electric fields, compute weights for GS07 CT algorithm
      // no correction to the EMFs is required; they are corrected later
      if (MAGNETIC_FIELDS_ENABLED) {
        pmb->pcoord->CenterWidth1(k,j,is,ie+1,dxw);
#pragma simd
        for (int i=is; i<=ie+1; ++i){
          ei_x1f(X1E3,k,j,i) = -flx(IBY,i); // flux(IBY) = (v1*b2 - v2*b1) = -EMFZ
          ei_x1f(X1E2,k,j,i) =  flx(IBZ,i); // flux(IBZ) = (v1*b3 - v3*b1) =  EMFY
          Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*flx(IDN,i)
                        / (dxw(i)*(wl(IDN,i) + wr(IDN,i)));
          Real tmp_min = std::min(0.5,v_over_c);
          w_x1f(k,j,i) = 0.5 + std::max(-0.5,tmp_min);
        }
      }
    }
  }

//----------------------------------------------------------------------------------------
// j-direction

  if (pmb->block_size.nx2 > 1) {
    // set the loop limits
    il=is, iu=ie, kl=ks, ku=ke;
    if (MAGNETIC_FIELDS_ENABLED) {
      if(pmb->block_size.nx3 == 1) // 2D
        il=is-1, iu=ie+1, kl=ks, ku=ke;
      else // 3D
        il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
    }
    for (int k=kl; k<=ku; ++k){
#pragma omp for schedule(static)
      for (int j=js; j<=je+1; ++j){

        // reconstruct L/R states at j
        if ((step == 1) && HYDRO_TIME_INTEGRATOR == "vl2") {
          pmb->precon->DonorCellX2(k,j,il,iu,w,bcc,wl,wr);
        } else {
          pmb->precon->HighResFuncX2(k,j,il,iu,w,bcc,wl,wr);
        }

        // compute fluxes at j
        RiemannSolver(k,j,il,iu,IVY,b2,wl,wr,flx); 

        // store fluxes
        if(k>=ks && k<=ke) {
          for(int n=0; n<NHYDRO; n++) {
#pragma simd
            for(int i=is; i<=ie; i++)
              x2flux(n,k,j,i)=flx(n,i);
          }
        }

        // store electric fields, compute weights for GS07 CT algorithm
        // no correction to the EMFs is required; they are corrected later
        if (MAGNETIC_FIELDS_ENABLED) {
          pmb->pcoord->CenterWidth2(k,j,il,iu,dxw);
#pragma simd
          for (int i=il; i<=iu; ++i){
            ei_x2f(X2E1,k,j,i) = -flx(IBY,i); // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
            ei_x2f(X2E3,k,j,i) =  flx(IBZ,i); // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
            Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*flx(IDN,i)
                          / (dxw(i)*(wl(IDN,i) + wr(IDN,i)));
            Real tmp_min = std::min(0.5,v_over_c);
            w_x2f(k,j,i) = 0.5 + std::max(-0.5,tmp_min);
          }
        }
      }
    }
  }

//----------------------------------------------------------------------------------------
// k-direction 

  if (pmb->block_size.nx3 > 1) {
    // set the loop limits
    il=is, iu=ie, jl=js, ju=je;
    if (MAGNETIC_FIELDS_ENABLED)
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;
#pragma omp for schedule(static)
    for (int k=ks; k<=ke+1; ++k){
      for (int j=jl; j<=ju; ++j){

        // reconstruct L/R states at k
        if ((step == 1) && HYDRO_TIME_INTEGRATOR == "vl2") {
          pmb->precon->DonorCellX3(k,j,il,iu,w,bcc,wl,wr);
        } else {
          pmb->precon->HighResFuncX3(k,j,il,iu,w,bcc,wl,wr);
        }

        // compute fluxes at k
        RiemannSolver(k,j,il,iu,IVZ,b3,wl,wr,flx);

        if(j>=js && j<=je) {
          for(int n=0; n<NHYDRO; n++) {
#pragma simd
            for(int i=is; i<=ie; i++)
              x3flux(n,k,j,i)=flx(n,i);
          }
        }

        // store electric fields, compute weights for GS07 CT algorithm
        // no correction to the EMFs is required; they are corrected later
        if (MAGNETIC_FIELDS_ENABLED) {
          pmb->pcoord->CenterWidth3(k,j,il,iu,dxw);
#pragma simd
          for (int i=il; i<=iu; ++i){
            ei_x3f(X3E2,k,j,i) = -flx(IBY,i); // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
            ei_x3f(X3E1,k,j,i) =  flx(IBZ,i); // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
            Real v_over_c = (1024.0)*(pmb->pmy_mesh->dt)*flx(IDN,i)
                          / (dxw(i)*(wl(IDN,i) + wr(IDN,i)));
            Real tmp_min = std::min(0.5,v_over_c);
            w_x3f(k,j,i) = 0.5 + std::max(-0.5,tmp_min);
          }
        }
      }
    }
  }

} // end of omp parallel region

  return;
}
