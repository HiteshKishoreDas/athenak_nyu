//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file turb.cpp
//  \brief Problem generator for turbulence
#include <cmath>
#include <iostream> // cout

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"


// User-defined history functions
void TurbulentHistory(HistoryData *pdata, Mesh *pm);

namespace{

struct ProblemData {

  Real gm1 = 2.0/3.0;
  bool Tfix_enabled = false;

  Real rho0 = 1.0;
  Real prs0 = 1.0;
  Real T0 = prs0 / rho0;

  Real xmin = -0.5;
  Real xmax = 0.5;

  // MHD variables
  Real beta = -1.0;
  int ifield = -1;

};
ProblemData prob_data;

void FatalProbInput(const std::string &message) {
  std::cout << "### FATAL ERROR in TRML_frame_tracking input" << std::endl
            << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void Tfix_source(Mesh *pm, const Real bdt) {

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;

  auto &size = pm->pmb_pack->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  const ProblemData data = prob_data;

  Real sum_Tvol = 0.0;
  Real sum_vol = 0.0;

  Kokkos::parallel_reduce(
      "T_collect",
      Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            pmbp->nmb_thispack * indcs.nx3 * indcs.nx2 * indcs.nx1),
      KOKKOS_LAMBDA(const int &idx, Real &block_Tvol, Real &block_vol) {
    const int nkji = indcs.nx3*indcs.nx2*indcs.nx1;
    const int nji = indcs.nx2*indcs.nx1;
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji + ks;
    int j = (idx - m*nkji - (k - ks)*nji)/indcs.nx1 + js;
    int i = (idx - m*nkji - (k - ks)*nji - (j - js)*indcs.nx1) + is;

    const Real density = u0(m, IDN, k, j, i);

    const Real ek = 0.5*(SQR(u0(m, IM1, k, j, i)) + SQR(u0(m, IM2, k, j, i)) +
             SQR(u0(m, IM3, k, j, i)))/density;

    const Real eint = u0(m, IEN, k, j, i) - ek;

    const Real cell_vol = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
    block_Tvol += eint / density * cell_vol;
    block_vol += cell_vol;

  }, Kokkos::Sum<Real>(sum_Tvol), Kokkos::Sum<Real>(sum_vol));

#if MPI_PARALLEL_ENABLED
  Real sums[2] = {sum_Tvol, sum_vol};
  MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  sum_Tvol = sums[0];
  sum_vol = sums[1];
#endif

  const Real T_avg = data.gm1 * sum_Tvol / sum_vol;
  const Real T_factor = data.T0 / T_avg;
  if (global_variable::my_rank == 0) {
    std::cout << "Volume-averaged temperature =" << T_avg << std::endl;
  }

  par_for("Tfix_cooling", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real density = u0(m, IDN, k, j, i);

    const Real ek = 0.5*(SQR(u0(m, IM1, k, j, i)) + SQR(u0(m, IM2, k, j, i)) +
             SQR(u0(m, IM3, k, j, i)))/density;
    const Real eint = u0(m, IEN, k, j, i) - ek;

    u0(m, IEN, k, j, i) = eint * T_factor + ek;
  });


} // Tfix_source


void ReadProbParameters(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;

  // TODO: Implement MHD version for turb_cond in the future
  // Just need to get the Tfix_source() implemented for MHD
  if ((pmbp->pmhd != nullptr)) {
    FatalProbInput("turb_cond hasn't have a MHD version implemented.");
  }

  if ((pmbp->phydro  == nullptr) && (pmbp->pmhd == nullptr)) {
    FatalProbInput("turb_cond requires a <hydro> or <mhd> block.");
  }
  if (!pmbp->phydro->peos->eos_data.is_ideal) {
    FatalProbInput("turb_cond requires an ideal-gas hydro EOS.");
  }

  ProblemData data;
  data.gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
  data.rho0 = pin->GetOrAddReal("problem", "rho_0", 1.0);
  data.prs0 = pin->GetOrAddReal("problem", "pgas_0", 1.0);
  data.xmin = pm->mesh_size.x1min;
  data.xmax = pm->mesh_size.x1max;

  data.Tfix_enabled = pin->GetOrAddBoolean("problem", "Tfix_enabled", true);

  if (data.rho0 <= 0.0 || data.prs0 <= 0.0) {
    FatalProbInput("Require rho_0 > 0 and pgas_0 > 0.");
  }
  data.T0 = data.prs0/data.rho0;

  if (pmbp->pmhd != nullptr){
    data.beta = pin->GetOrAddReal("problem","beta",1.0);
    data.ifield = pin->GetOrAddInteger("problem","ifield",2);
  }
  prob_data = data;

  if (global_variable::my_rank == 0) {
    std::cout << "turb_cond: rho0=" << data.rho0
              << ", T0=" << data.T0
              << ", Tfix_enabled=" << (data.Tfix_enabled? "on" : "off")
              << std::endl;
  }
} // ReadProbParameters



} // namespace

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::Turb_()
//  \brief Problem Generator for turbulence

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {

  ReadProbParameters(pin, pmy_mesh_);
  const ProblemData data = prob_data;
  if (data.Tfix_enabled) {
    user_srcs = true;
    user_srcs_func = Tfix_source;
  }

  user_hist_func = TurbulentHistory;

  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
       << "Turbulence problem generator can only be run with Hydro and/or MHD, but no "
       << "<hydro> or <mhd> block in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  // enroll user history function
  user_hist_func = TurbulentHistory;

  // capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;


  // Initialize Hydro variables -------------------------------
  if (pmbp->phydro != nullptr) {
    auto &u0 = pmbp->phydro->u0;
    EOS_Data &eos = pmbp->phydro->peos->eos_data;

    // Set initial conditions
    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = data.rho0;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;
      if (eos.is_ideal) {
        u0(m,IEN,k,j,i) = data.prs0/data.gm1 +
           0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
           SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
      }
    });
  }

  // Initialize MHD variables ---------------------------------
  if (pmbp->pmhd != nullptr) {
    if (data.ifield != 1 && data.ifield != 2) {
      std::cout << "### FATAL ERROR in " << __FILE__
                << " at line " << __LINE__ << std::endl
                << "Invalid <problem>/ifield = " << data.ifield
                << ", allowed values are 1 (zero-net-flux Bz) or 2 (uniform Bz)."
                << std::endl;
      exit(EXIT_FAILURE);
    }

    Real x1size = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
    Real kx = 2.0*(M_PI/x1size);
    auto &u0 = pmbp->pmhd->u0;
    auto &b0 = pmbp->pmhd->b0;
    auto &size = pmbp->pmb->mb_size;
    // EOS_Data &eos = pmbp->pmhd->peos->eos_data;

    Real B0 = std::sqrt(2.0*data.prs0/data.beta);

    // Set initial conditions
    par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      u0(m,IDN,k,j,i) = data.rho0;
      u0(m,IM1,k,j,i) = 0.0;
      u0(m,IM2,k,j,i) = 0.0;
      u0(m,IM3,k,j,i) = 0.0;

      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

      if (data.ifield == 1) {
        // zero-net-flux Bz
        b0.x1f(m,k,j,i) = 0.0;
        b0.x2f(m,k,j,i) = 0.0;
        b0.x3f(m,k,j,i) = B0*std::sin(kx*x1v);
        if (i==ie) {b0.x1f(m,k,j,i+1) = 0.0;}
        if (j==je) {b0.x2f(m,k,j+1,i) = 0.0;}
        if (k==ke) {b0.x3f(m,k+1,j,i) = B0*std::sin(kx*x1v);}
      } else if (data.ifield == 2) {
        // constant Bz
        b0.x1f(m,k,j,i) = 0.0;
        b0.x2f(m,k,j,i) = 0.0;
        b0.x3f(m,k,j,i) = B0;
        if (i==ie) {b0.x1f(m,k,j,i+1) = 0.0;}
        if (j==je) {b0.x2f(m,k,j+1,i) = 0.0;}
        if (k==ke) {b0.x3f(m,k+1,j,i) = B0;}
      }

      Real bz_cc = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
      u0(m,IEN,k,j,i) = data.prs0/data.gm1 + 0.5*bz_cc*bz_cc +
          0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
          SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
    });
  }

  return;
}


//----------------------------------------------------------------------------------------
// Function for computing history variables
// 0 = < T >
void TurbulentHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 1;
  pdata->label[0] = "Tsumvol";

  const ProblemData data = prob_data;

  // capture class variabels for kernel
  auto &u0_ = pm->pmb_pack->phydro->u0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  int &nhist_ = pdata->nhist;

  // loop over all MeshBlocks in this pack
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("HistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    // compute n,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;

    array_sum::GlobalSum hvars;

    // calculate mean B

    const Real density = u0_(m, IDN, k, j, i);
    const Real ek = 0.5*(SQR(u0_(m, IM1, k, j, i)) + SQR(u0_(m, IM2, k, j, i)) +
             SQR(u0_(m, IM3, k, j, i)))/density;
    const Real eint = u0_(m, IEN, k, j, i) - ek;

    hvars.the_array[0] = eint / density * data.gm1 * vol;

    // fill rest of the_array with zeros, if nhist < NHISTORY_VARIABLES
    for (int n=nhist_; n<NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }

    // sum into parallel reduce
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();

  // store data into hdata array
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }
  return;
}
