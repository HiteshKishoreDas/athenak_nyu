//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mass_removal_test.cpp
//  \brief Problem generator for testing mass removal
#include <iostream> // cout
#include <stdexcept>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"

#include "units/units.hpp"

//----------------------------------------------------------------------------------------
//! \struct pgen_trml
//! \brief Data structure holding all problem-specific parameters for the TRML simulation.
//!
//! This structure stores configuration parameters read from the input file and computed
//! derived quantities needed throughout the simulation. It is allocated once at
//! initialization and persists for the entire simulation run.

struct pgen_trml {

  // ====================================================================================
  // BOX DIMENSIONS
  // ====================================================================================

  Real ztop;                    //!< Top z-coordinate of the simulation box
  Real zbot;                    //!< Bottom z-coordinate of the simulation box

  // ====================================================================================
  // EQUATION OF STATE AND FLOOR VALUES
  // ====================================================================================
  Real gamma_adi;              //!< Adiabatic index (ratio of specific heats, gamma)
  Real dfloor;                 //!< Density floor (minimum allowed density)
  Real pfloor;                 //!< Pressure floor (minimum allowed pressure)
  Real tfloor;                 //!< Pressure floor (minimum allowed pressure)

  // ====================================================================================
  // INITIAL STATE PARAMETERS
  // ====================================================================================
  Real rho_0;                  //!< Reference density (hot phase density in code units)
  Real chi;               //!< Temperature contrast ratio (T_hot/T_cold = rho_cold/rho_hot)

  // ====================================================================================
  // INITIAL CLOUD PARAMETERS
  // ====================================================================================
  Real LboxR;                  //!< L_box/R_cloud ratio
  Real transition_R;           //!< R_cloud/transition layer thickness 

  Real grain_porosity = 1.0;

  // ====================================================================================
  // DUST MODEL PARAMETERS
  // ====================================================================================

  bool dust_model;             //! Flag to turn the dust model on/off 
  Real Z_gas;                  //! Initial Gas metallicity in Z_solar
  Real D_Z_init;               //! Initial dust-to-gas ratio
  Real dust_a;                 //! Min dust size
  Real dust_b;                 //! Max dust size
  Real dust_exp;               //! Dust size distribution exponent
  int Ndust_bins;              //! Number of dust size bins
  DvceArray1D<Real> dbins; //! Dust size bin centers (device view)
  DvceArray1D<Real> dedges; //! Dust size bin edges (device view)
  DvceArray1D<Real> dedges_pad; //! Dust size bin edges with ghosts (device view)

  Real rho_gr;                  //! Dust grain density in g/cc

  Real Z_solar;                //! Solar metallicity

  // ====================================================================================
  // TEMPERATURE THRESHOLDS
  // ====================================================================================
  Real T_cold;                 //!< Cold phase temperature (= pgas_0/rho_0/contrast)
  Real T_hot;                  //!< Hot phase temperature (= pgas_0/rho_0)
  Real T_floor;                //!< Temperature floor for adjusting cold gas

  // ====================================================================================
  // ADAPTIVE MESH REFINEMENT THRESHOLDS
  // ====================================================================================
  Real ddens_threshold; 

};

// Global pointer to the TRML parameters structure.
// Allocated in ProblemGenerator::UserProblem() and persists for the simulation lifetime.
// 
pgen_trml* ptrml = new pgen_trml();


//! \brief Add user-defined source terms (cooling/heating) to conserved variables.
//! Called every RK substep via user_srcs_func.
void AddUserSrcs(Mesh *pm, const Real bdt);

//! \brief Apply dust model to each cell
//! Includes thermal sputtering, accretion, shattering and coagulation

template <typename T>
KOKKOS_INLINE_FUNCTION
int searchsorted(const Real* a, int n, const T& x);

void AddDustSource(Mesh *pm, const Real bdt);
KOKKOS_INLINE_FUNCTION
Real sputtering(Real n_H, Real T, Real Z, Real grain_porosity);
KOKKOS_INLINE_FUNCTION
Real accretion(Real n_H, Real T, Real Z, Real grain_porosity);
KOKKOS_INLINE_FUNCTION
void rebin(Real* d_arr, const Real* edges, const int n, Real shift, Real rho_gr);


void TurbulentHistory(HistoryData *pdata, Mesh *pm);

void RefinementCondition(MeshBlockPack* pmbp);

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  // Enroll user functions 
  user_srcs_func = AddUserSrcs;
  user_hist_func = TurbulentHistory;
  user_ref_func = RefinementCondition;

  bool is_hydro = (pmbp->phydro != nullptr) ? true : false;
  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
  EOS_Data &eos = (is_mhd) ?
                pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  int &nscalars = (is_mhd) ?
                pmbp->pmhd->nscalars : pmbp->phydro->nscalars;
  int &nfluid = (is_mhd) ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;

  // Get code unit variables for temperature conversions
  Real mu = 1.0;
  Real KELVIN = 1.0;
  if (pmbp->punit != nullptr) {
    KELVIN = pmbp->punit->temperature_cgs();

    printf("==========================================\n");
    printf("Code unit conversion: mu = %e, KELVIN = %e\n", mu, KELVIN);
    printf("==========================================\n");

  } else if (global_variable::my_rank == 0) {
    std::cout << "WARNING: <units> block missing; assuming mu=1 and KELVIN=1 "
              << "for temperature conversions." << std::endl;
    std::exit(1);
  }
  
  // Read general parameters from input file
  ptrml->gamma_adi         = eos.gamma;
  ptrml->dfloor            = eos.dfloor;
  ptrml->pfloor            = eos.pfloor;
  ptrml->tfloor            = eos.tfloor;

  Real rho_0               = pin->GetReal("problem", "rho_0");
  ptrml->rho_0             = rho_0;

  ptrml->T_floor           = ptrml->tfloor * KELVIN;
  ptrml->T_hot             = pin->GetOrAddReal("problem", "T_hot", 1e6);
  ptrml->T_cold            = pin->GetOrAddReal("problem", "T_cold", 1e4);
  Real T_hot               = ptrml->T_hot;

  ptrml->chi               = ptrml->T_hot/ptrml->T_cold;

  // get ztop and bottom
  ptrml->ztop              = pin->GetReal("mesh", "x3max");
  ptrml->zbot              = pin->GetReal("mesh", "x3min");

  // For dust model
  ptrml->dust_model        = pin->GetBoolean("problem", "dust_model");
  ptrml->Z_gas             = pin->GetReal("problem", "Z_gas");  
  ptrml-> D_Z_init         = pin->GetReal("problem", "D_Z_init");  
  ptrml->dust_a            = pin->GetReal("problem", "dust_a"); // Min dust size
  ptrml->dust_b            = pin->GetReal("problem", "dust_b"); // Max dust size
  ptrml->dust_exp          = pin->GetReal("problem", "dust_exp"); // Dust size distribution exponent
  ptrml->rho_gr          = pin->GetReal("problem", "rho_gr"); // Dust grain density


  ptrml-> Z_solar          = 0.0134;  
  Real Z_gas               = ptrml->Z_gas;
  Real D_Z_init            = ptrml->D_Z_init;
  Real Z_solar             = ptrml->Z_solar;

  ptrml->LboxR              = pin->GetReal("problem", "LboxR");
  ptrml->transition_R       = pin->GetReal("problem", "transition_R");

  // Number of dust bins: total scalars minus metal and cloud tracer entries.
  // Must be initialized even on restarts because source terms/history use it.
  const int Ndust_bins = nscalars - 2;
  if (Ndust_bins <= 0) {
    throw std::runtime_error("Ndust_bins <= 0. Check <hydro>/nscalars in the input (need metal + tracer + >=1 dust bin)");
  }
  ptrml->Ndust_bins = Ndust_bins;

  // Create dust-bin midpoint array used by source terms.
  {
    const Real dust_a = ptrml->dust_a;
    const Real dust_b = ptrml->dust_b;

    DvceArray1D<Real> dbins("dust_bins", Ndust_bins);
    DvceArray1D<Real> dedges("dust_edges", Ndust_bins+1);
    DvceArray1D<Real> dedges_pad("dust_edges_pad", Ndust_bins+3);
    auto dbins_h = Kokkos::create_mirror_view(dbins);
    auto dedges_h = Kokkos::create_mirror_view(dedges);
    auto dedges_ph = Kokkos::create_mirror_view(dedges_pad);

    const Real r = std::pow(dust_b / dust_a, 1.0 / Ndust_bins);

    Real a_i = dust_a;
    for (int id = 0; id < Ndust_bins+1; ++id) {
      dedges_h(id) = a_i;
      dedges_ph(id+1) = a_i;
      a_i *= r;
    }
    for (int id = 0; id < Ndust_bins+1; ++id) {
      dbins_h(id) = 0.5*(dedges_h(id)+dedges_h(id+1));
    }
    
    Kokkos::deep_copy(dbins, dbins_h);
    Kokkos::deep_copy(dedges, dedges_h);
    Kokkos::deep_copy(dedges_pad, dedges_ph);

    ptrml->dbins = dbins;
    ptrml->dedges= dedges;
    ptrml->dedges_pad= dedges_pad;
  }


  // Read the density gradient threshold for refinement
  ptrml->ddens_threshold = pin->GetReal("problem", "ddens_max");

  if (restart) return;


  // Capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  // Initialize Hydro/MHD variables -------------------------------
  auto &w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto &u0 = (is_mhd) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  Real gm1 = eos.gamma - 1.0;
  auto &size = pmbp->pmb->mb_size;

  Real box = abs(ptrml->ztop - ptrml->zbot);
  Real radius = box/ptrml->LboxR;
  Real smoothing_thickness = radius/ptrml->transition_R;
  Real chi = ptrml->chi;

  int nmb1 = pmbp->nmb_thispack - 1;

  const int nmetal = nfluid;                 // index for metal
  const int ntracer = nfluid+1;              // index for tracer for cloud
  const int ndusti = nfluid+2;               // first dust bin num scalar

  printf("======================\n");
  printf("nmetal = %d, ntracer = %d, ndusti = %d\n", nmetal, ntracer, ndusti);
  printf("Ndust_bins = %d\n", Ndust_bins);
  printf("======================\n");

  // Set initial conditions
  if (global_variable::my_rank == 0) {
    std::cout << "Now initializing hydro/MHD variables." << "\n";
  }

  par_for("pgen_turb", DevExeSpace(),0,nmb1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;

    int nx1 = indcs.nx1;
    int nx2 = indcs.nx2;
    int nx3 = indcs.nx3;

    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    Real shape = 0.5 * (1.0+std::tanh((radius-std::sqrt(x1v*x1v + x2v*x2v + x3v*x3v))/smoothing_thickness));

    w0(m,IDN,k,j,i) = rho_0*(1.0 + (chi-1.0)*shape);
    w0(m,IVX,k,j,i) = 0.0;
    w0(m,IVY,k,j,i) = 0.0;
    w0(m,IVZ,k,j,i) = 0.0;
    if (eos.is_ideal) {
      w0(m,IEN,k,j,i) = T_hot * rho_0/ KELVIN  / gm1;
    }

    // add passive scalars
    Real tmp_dust_shape = 0.0;
    if(nscalars>0){

      // First scalar has to be metallicity
      w0(m,nmetal,k,j,i) = Z_gas * Z_solar * (0.1 + 0.9 * shape);

      // Cloud material tracer
      w0(m,ntracer,k,j,i) = 1.0 * shape;
      
      tmp_dust_shape = D_Z_init * Z_gas * Z_solar * shape / float(Ndust_bins);
      for (int id=0; id<Ndust_bins; id++){
        // Dust
        w0(m,ndusti+id,k,j,i) = tmp_dust_shape;
        // The D_tot comes out to be D_Z_init * Z_gas * Z_solar, i.e D_Z_init * Z_g
      }
    }
  });
  // Convert primitives to conserved
  if (is_hydro) {
    pmbp->phydro->peos->PrimToCons(w0, u0, is, ie, js, je, ks, ke);
  }


  return;
}

//! \fn void AddUserSrcs()
//! \brief Add User Source Terms
// NOTE source terms must all be computed using primitive (w0) and NOT conserved (u0) vars
void AddUserSrcs(Mesh *pm, const Real bdt) {
  // MeshBlockPack *pmbp = pm->pmb_pack;
  // bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
  // DvceArray5D<Real> &u0 = (is_mhd) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  // const DvceArray5D<Real> &w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  // const EOS_Data &eos_data = (is_mhd) ?
  //                 pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  // if (ptrml->adjust_temp_floor) AdjustTempTFloor(pm,bdt,u0,w0,eos_data);
  if (ptrml->dust_model){
    AddDustSource(pm, bdt);
  }
  return;
}

KOKKOS_INLINE_FUNCTION
Real sputtering(Real n_H, Real T, Real Z, Real grain_porosity){

  Real da_dt = -1.0; // um/Myr
  da_dt *= n_H/1.0; // cm^-3
  da_dt /= 1. + pow(2.0e6 / T, 2.5); // T in K
  da_dt *= Z; // in Z_solar
  da_dt *= pow(grain_porosity, -2./3.);

  return da_dt;
}

KOKKOS_INLINE_FUNCTION
Real accretion(Real n_H, Real T, Real Z, Real grain_porosity){

  Real da_dt = 1.8622e-4; // um/Myr
  da_dt *= n_H/1.0e3; // cm^-3
  da_dt /= sqrt(T/10.0); // T in K
  da_dt *= Z; // in Z_solar
  da_dt *= pow(grain_porosity, -2./3.);

  return da_dt;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
int searchsorted(const Real* a, const int n, const T& x) {
  if (n <= 0) return 0;
  int lo = 0;
  int hi = n; // insertion range [0, n]

  if (x <= a[lo]) return 0;
  if (x > a[n-1]) return hi;

  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (a[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

KOKKOS_INLINE_FUNCTION
Real dist_num_integral(const Real* d_arr, const Real* edges, const int n,
                       Real a1, Real a2, Real rho_gr, bool mass=false){
  // Integrate the distribution d_arr on edges between a1 and a2
  // These arrays are padded arrays

  int i1 = searchsorted(edges, n, a1)-1;
  int i2 = searchsorted(edges, n, a2)-1;

  Real K_dust = (4.0/3.0)*M_PI*rho_gr;

  Real integral = 0.0;   
  // Both are in the same bin
  if (i1==i2){

    if (mass) integral += 0.25 * d_arr[i1] * K_dust * (pow(a2, 4) - pow(a1, 4));
    else integral +=d_arr[i1]*(a2-a1);

  }

  // Add the left edge
  if (mass) integral += 0.25 * d_arr[i1] * K_dust * (pow(edges[i1+1], 4) - pow(a1, 4));
  else integral += d_arr[i1] * (edges[i1+1] - a1);
  

  // Add the right edge
  if (mass) integral += 0.25 * d_arr[i2] * K_dust * (pow(a2, 4) - pow(edges[i2], 4));
  else integral += d_arr[i2] * (a2 - edges[i2]);

  // Add the rest
  if (i2 > i1+1){

    if (mass) {
      for (int i=0; i<a2; i++){
        integral += 0.25 * d_arr[i] * K_dust * (pow(edges[i+1], 4) - pow(edges[i], 4));
      }
    }
    else {
      for (int i=0; i<a2; i++){
        integral += d_arr[i] * (edges[i+1] - edges[i]);
      }
    }

  }

  return integral;
}

KOKKOS_INLINE_FUNCTION
void rebin(Real* d_arr, const Real* edges, const int n, Real shift, Real rho_gr) {
  // We take in the arrays already padded with ghosts
  Real K_dust = (4.0/3.0)*M_PI*rho_gr;
  constexpr int kMaxRebinEntries = 34;  // Ndust_bins <= 32, plus two ghost entries
  if (n > kMaxRebinEntries) return;

  Real mass_in_bin[kMaxRebinEntries];
  Real shifted_edges[kMaxRebinEntries];
  for (int i=0; i<n; i++){
    mass_in_bin[i] = 0.0;
    shifted_edges[i] = edges[i]+shift;
  }

  for (int i=0; i<n-1; i++){
    mass_in_bin[i] = dist_num_integral(
      d_arr, shifted_edges, n, edges[i], edges[i+1], rho_gr, true);
  }

  for (int i=1; i<n-1; i++){
    d_arr[i] = mass_in_bin[i]/K_dust/(pow(edges[i+1], 4)-pow(edges[i], 4));
  }
  d_arr[0] = mass_in_bin[0];
  d_arr[n-1] = mass_in_bin[n-1];

  return;
}


//! \fn void AddDustSource()
//! \brief Apply dust model to each cell
//! Includes thermal sputtering, accretion, shattering and coagulation
void AddDustSource(Mesh *pm, const Real bdt){
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  // auto &size = pmbp->pmb->mb_size;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;
  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
  auto &u0 = (is_mhd) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  const auto &w0 = (is_mhd) ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  const EOS_Data &eos_data = (is_mhd) ?
                  pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  int &nfluid = (is_mhd) ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;

  // int nx1 = indcs.nx1;
  // int nx2 = indcs.nx2;
  // int nx3 = indcs.nx3;

  // const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  // const int nkji = nx3*nx2*nx1;
  // const int nji  = nx2*nx1;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;

  Real bta_time = bdt/pm->dt;
  Real use_e = eos_data.use_e;
  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;

  // Get code unit variables for temperature conversions
  Real KELVIN = 1.0;
  if (pmbp->punit != nullptr) {
    KELVIN = pmbp->punit->temperature_cgs();
  } else if (global_variable::my_rank == 0) {
    std::cout << "ERROR: <units> block missing..." << std::endl;
    std::exit(1);
  }

  int nmetal = nfluid; // index for metal
  // int ntracer = nfluid+1; // index for tracer for cloud
  int ndusti = nfluid+2; // first dust bin scalar
  const int Ndust_bins = ptrml->Ndust_bins; // Number of dust bins
  constexpr int kMaxDustBins = 32; // upper bound supported in device stack arrays
  if (Ndust_bins > kMaxDustBins) {
    if (global_variable::my_rank == 0) {
      std::cout << "ERROR: Ndust_bins=" << Ndust_bins
                << " exceeds kMaxDustBins=" << kMaxDustBins << std::endl;
    }
    std::exit(1);
  }


  Real Z_solar = ptrml->Z_solar;
  auto dust_bins = ptrml->dbins;

  Real size_d1 = 1.0; // in nm
  Real size_d2 = 10.0; // in nm

  Real sig_DL = 10; // in km/s
  Real sig_DS = 0.1; // in km/s
  Real F_coag = 0.5; // fudge factor for coagulation

  Real grain_porosity = ptrml->grain_porosity;
  Real rho_gr = ptrml->rho_gr;

  Real edges_pad[kMaxDustBins+3];
  for (int i=0; i<Ndust_bins+3; i++) edges_pad[i] = ptrml->dedges_pad(i);

  par_for("user_source", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dens = w0(m, IDN, k, j, i); // in cm^-3
    Real temp = (w0(m,IEN,k,j,i) * gm1) / dens * KELVIN;

    Real rho_di = 0.0;
    Real dust_rate[kMaxDustBins];
    for (int id = 0; id < kMaxDustBins; ++id) dust_rate[id] = 0.0;
    
    // To prevent NaNs
    Real Z_tol = 1e-20; 
    Real clamp_dust_tot = 1.0; // tmp var used for clamping later

    Real Z_local = w0(m,nmetal,k,j,i);
    Z_local = Kokkos::max(Z_tol, Z_local);

    int kk = k - ks;
    int jj = j - js;
    int ii = i - is;
    int cell_idx = m*nkji + kk*nji + jj*nx1 + ii;

    Real dust_arr[kMaxDustBins+2];

    // Total dust ratio across size bins
    dust_arr[0] = 0.0;
    dust_arr[Ndust_bins+1] = 0.0;
    for (int id=0; id<Ndust_bins; id++){
      dust_arr[id+1] = u0(m, ndusti+id, k, j, i);
    }

    Real rate_Z = 0.0, tmp_rate=0.0;
    Real dust_i=0, dust_in=0, dust_ip=0;

    Real tot_shift = sputtering(dens, temp, Z_local/Z_solar, grain_porosity);
    tot_shift += accretion(dens, temp, Z_local/Z_solar, grain_porosity);
    tot_shift *= bdt;

    Real delta_dmass = -dist_num_integral(
      dust_arr, edges_pad, n, ptrml->dust_a, ptrml->dust_b, rho_gr, true);
  
    rebin(dust_arr, edges_pad, Ndust_bins+2, tot_shift, rho_gr);

    delta_dmass += dist_num_integral(
      dust_arr, edges_pad, n, ptrml->dust_a, ptrml->dust_b, rho_gr, true);

    // Loop through the dust bins to calculate the rates
    for (int id=0; id<Ndust_bins; id++){

      // current dust amount in the bin
      dust_i = w0(m,ndusti+id,k,j,i);
      dust_i = Kokkos::max(Z_tol, dust_i);

      // Is there a next bin
      if (id<(Ndust_bins-1)) {
        dust_in = 1.0;
      }
      else dust_in = 0.0;

      // Is there a previous bin
      if (id>0) {
        dust_ip = 1.0;
      }
      else dust_ip = 0.0;

      // dust mass per unit vol.
      rho_di = dust_i * dens;


      // tmp_rate = 0.0;
      
      // //* Thermal sputtering
      // Real t_sp = 70.0; // in Myr
      // t_sp *= (1.e-3/dens);
      // t_sp *= 1. + pow((temp/2.0e6), -2.5);

      // tmp_rate = -rho_di / (t_sp * dust_bins(id));
      // dust_rate[id] += tmp_rate;
      // rate_Z += -tmp_rate;

      // //* Accretion
      // Real t_ac = 200.0; // in Myr
      // t_ac *= (20.0/dens);
      // t_ac *= pow((50.0/temp), 0.5);
      // t_ac *= Z_local/Z_solar; // in Z_sol 

      // t_ac /= 1-(dust_tot/Z_local);

      // tmp_rate = rho_di / (t_ac * dust_bins(id));
      // dust_rate[id] += tmp_rate;

      rate_Z += -tmp_rate;

      // Convert rate_Z to metal mass
      rate_Z *= 1.0;  //! Assuming Z_dust ~ 1

      // //* Shattering out of this bin, into previous one
      // // From Dubois 2024
      // Real t_shatt = 54.0; // in Myr
      // t_shatt *= (1.0/dens); // in cm^-3
      // t_shatt *= (s_i/3.0);  // in cm^-3
      // t_shatt *= (0.01/dust_i); 
      // t_shatt *= (10/sig_DL);  // in km/s

      // tmp_rate = rho_di / (t_shatt * dust_bins(id));
      // tmp_rate *= id>0;

      // dust_rate[id] += -tmp_rate;
      // dust_rate[id-1] += tmp_rate;


      // //* Coagulation out of this bin, into next one
      // // From Dubois 2024
      // Real t_coag = 0.27; // in Myr
      // t_coag *= (s_i/3.0);  // in cm^-3
      // t_coag *= (1.0e3/dens); // in cm^-3
      // t_coag *= (0.01/dust_i); 
      // t_coag *= (0.1/sig_DS);  // in km/s
      // t_coag *= F_coag;

      // tmp_rate = rho_di / (t_coag * dust_bins(id) / 0.05);
      // tmp_rate *= id<(Ndust_bins-1);

      // if (tmp_rate!=0.0){
      //   dust_rate[id] += -tmp_rate;
      //   dust_rate[id+1] += tmp_rate;
      // }

    }

    //! These are w0 * dens
    u0(m, nmetal, k, j, i) += bdt * rate_Z;

    Real dust_tot = 0.0;
    for (int id=0; id<Ndust_bins; id++){
      u0(m, ndusti+id, k, j, i) += bdt * dust_rate[id];
      dust_tot += u0(m, ndusti+id, k, j, i); // for later
    }

    // We should check that scalars are guaranteed to be in [0,1] after all source terms are added.
    // Clamp metallicity 
    u0(m, nmetal, k, j, i) = Kokkos::clamp(u0(m, nmetal, k, j, i), 0.0, dens);
    // Clamp total dust-to-gas ratio to [0,1]
    // Dust to metal ratio is not clamped, as all of gas metals can be locked in dust
    // Also dust might have been created in a metal-rich area and moved...
    clamp_dust_tot = Kokkos::clamp(dust_tot, 0.0, dens);

    // Scaled to sum to clamp_dust_tot 
    for (int id=0; id<Ndust_bins; id++){
      u0(m, ndusti+id, k, j, i) *= clamp_dust_tot/dust_tot;
    }

  });

  return;
}

// TODO: Update the history to work with N-bins of dust
void TurbulentHistory(HistoryData *pdata, Mesh *pm) {
  auto &w0_ = pm->pmb_pack->phydro->w0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  int &nhist_ = pdata->nhist;

  bool is_mhd = (pm->pmb_pack->pmhd != nullptr) ? true : false;
  const EOS_Data &eos_data = (is_mhd) ?
                  pm->pmb_pack->pmhd->peos->eos_data : pm->pmb_pack->phydro->peos->eos_data;
  int &nfluid_ = (is_mhd) ? pm->pmb_pack->pmhd->nmhd : pm->pmb_pack->phydro->nhydro;

  int nmetal_ = nfluid_; // index for metal
  // int ntracer_ = nfluid_+1; // index for tracer for cloud
  int ndusti_ = nfluid_+2; // first dust bin scalar
  const int Ndust_bins = ptrml->Ndust_bins; // Number of dust bins

  int count = 0;
  pdata->label[count] = "U^2"; count++;
  pdata->label[count] = "Mcold"; count++;
  pdata->label[count] = "T_sum"; count++;

  bool dust_model = ptrml->dust_model;

  int DUST_HIST = count;
  if (dust_model){
    pdata->label[count] = "Zgas"; count++;
    for (int id = 0; id < ptrml->Ndust_bins; ++id) {
      pdata->label[count] = "D" + std::to_string(id);
      count++;
    }
  }

  pdata->nhist = count;

  // loop over all MeshBlocks in this pack
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is; int nx1 = indcs.nx1;
  int js = indcs.js; int nx2 = indcs.nx2;
  int ks = indcs.ks; int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;

  Real gamma = eos_data.gamma;
  Real gm1 = gamma - 1.0;
  Real KELVIN = 1.0;
  if (pm->pmb_pack->punit != nullptr) {
    KELVIN = pm->pmb_pack->punit->temperature_cgs();
  } else if (global_variable::my_rank == 0) {
    std::cout << "ERROR: <units> block missing..." << std::endl;
    std::exit(1);
  }
  Real T_cold = ptrml->T_cold;

  array_sum::GlobalSum sum_this_mb;
  // store data into hdata array
  for (int n=0; n<NREDUCTION_VARIABLES; ++n) {
    sum_this_mb.the_array[n] = 0.0;
  }

  Kokkos::parallel_reduce("HistSums",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    // compute n,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real dens = w0_(m, IDN, k, j, i); 
    Real temp = (w0_(m,IEN,k,j,i) * gm1) / dens * KELVIN;

    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    Real dx_squared = size.d_view(m).dx1 * size.d_view(m).dx1;
    dens *= vol; // This is mass from now

    array_sum::GlobalSum hvars;
    hvars.the_array[0] += ((w0_(m,IVX,k,j,i)*w0_(m,IVX,k,j,i))
                        + (w0_(m,IVY,k,j,i)*w0_(m,IVY,k,j,i))
                        + (w0_(m,IVZ,k,j,i)*w0_(m,IVZ,k,j,i)))*vol;

    hvars.the_array[1] += (temp < 2.0*T_cold ? dens : 0);
    hvars.the_array[2] += temp;
    // Note that here dens is actually dens*vol

    if (dust_model){
      int dust_hist = DUST_HIST; // local, modifiable copy of starting dust history index
      hvars.the_array[dust_hist] += dens * w0_(m, nmetal_, k, j, i);
      ++dust_hist;

      for (int id=0; id < Ndust_bins; id++){
        hvars.the_array[dust_hist + id] += dens * w0_(m, ndusti_+id, k, j, i);
      }
    }
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

//===========================================================================//
//                              Refinement                                   //
//===========================================================================//

// Refine region based on density gradient threshold
void RefinementCondition(MeshBlockPack* pmbp) {
  Mesh *pmesh       = pmbp->pmesh;
  int nmb           = pmbp->nmb_thispack;
  int mbs           = pmesh->gids_eachrank[global_variable::my_rank];
  auto &refine_flag = pmesh->pmr->refine_flag;
  // auto &multi_d     = pmesh->multi_d;
  // auto &three_d     = pmesh->three_d;
  auto &indcs       = pmesh->mb_indcs;
  int &is = indcs.is, nx1 = indcs.nx1;
  int &js = indcs.js, nx2 = indcs.nx2;
  int &ks = indcs.ks, nx3 = indcs.nx3;
  const int nkji = nx3 * nx2 * nx1;
  const int nji  = nx2 * nx1;
  auto &u0       = pmbp->phydro->u0;
  auto &w0       = pmbp->phydro->w0;

  auto &ddens_thresh = ptrml->ddens_threshold;

  par_for_outer("UserRefineCond",DevExeSpace(), 0, 0, 0, (nmb-1),
  KOKKOS_LAMBDA(TeamMember_t tmember, const int m) {

    Real team_ddmax;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkji),
    [=](const int idx, Real& ddmax) {
      int k = (idx)/nji;
      int j = (idx - k*nji)/nx1;
      int i = (idx - k*nji - j*nx1) + is;
      j += js;
      k += ks;

      // Calculate density gradient
      Real d2 = (SQR(u0(m,IDN,k,j,i+1) - u0(m,IDN,k,j,i-1))
               + SQR(u0(m,IDN,k,j+1,i) - u0(m,IDN,k,j-1,i))
               + SQR(u0(m,IDN,k+1,j,i) - u0(m,IDN,k-1,j,i)));
      ddmax = fmax((sqrt(d2)/u0(m,IDN,k,j,i)), ddmax);

      // Calculate pressure gradient
      Real p2 = (SQR(w0(m,IEN,k,j,i+1) - w0(m,IEN,k,j,i-1))
               + SQR(w0(m,IEN,k,j+1,i) - w0(m,IEN,k,j-1,i))
	       + SQR(w0(m,IEN,k+1,j,i) - w0(m,IEN,k-1,j,i)));
      ddmax = fmax((sqrt(p2)/w0(m,IEN,k,j,i)), ddmax);
    },Kokkos::Max<Real>(team_ddmax));

    if (team_ddmax > ddens_thresh) {refine_flag.d_view(m+mbs) = 1;}
    if (team_ddmax < 0.1*ddens_thresh) {refine_flag.d_view(m+mbs) = -1;}

  });

  // sync host and device
  refine_flag.template modify<DevExeSpace>();
  refine_flag.template sync<HostMemSpace>();
}
