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

#include <Kokkos_Random.hpp>

//----------------------------------------------------------------------------------------
//! \struct pgen_trml
//! \brief Data structure holding all problem-specific parameters for the TRML simulation.
//!
//! This structure stores configuration parameters read from the input file and computed
//! derived quantities needed throughout the simulation. It is allocated once at
//! initialization and persists for the entire simulation run.

struct pgen_trml {
  ParameterInput *pin = nullptr; //!< Live parameter store used for restart-safe updates

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
  // COOLING/HEATING PARAMETERS
  // ====================================================================================

  bool cgm_cooling;         //!< Is CGMCooling enabled?
  Real hrate;               //!< heating rate

  // ====================================================================================
  // INITIAL CLOUD PARAMETERS
  // ====================================================================================

  bool cloud_flag;             //!< Flag to turn the cloud on/off
  Real LboxR;                  //!< L_box/R_cloud ratio
  Real transition_R;           //!< R_cloud/transition layer thickness 

  Real grain_porosity = 1.0;

  // ====================================================================================
  // DUST MODEL PARAMETERS
  // ====================================================================================

  bool dust_model;             //!< Flag to turn the dust model on/off 
  Real Z_gas;                  //!< Initial Gas metallicity in Z_solar
  Real D_Z_init;               //!< Initial dust-to-gas ratio
  Real dust_a;                 //!< Min dust size
  Real dust_b;                 //!< Max dust size
  Real dust_exp;               //!< Dust size distribution exponent
  int Ndust_bins;              //!< Number of dust size bins
  DvceArray1D<Real> dbins;     //!< Dust size bin centers (device view)
  DvceArray1D<Real> dedges;    //!< Dust size bin edges (device view)
  DvceArray1D<Real> dedges_pad;//!< Dust size bin edges with ghosts (device view)

  Real rho_gr;                 //!< Dust grain density in g/cc

  Real Z_solar;                //!< Solar metallicity

  int vturb_est_ncell;         //!< Number of cells over which vturb is estimated
  bool coag_shatt_flag;        //!< Flag for coagulation & shattering

  // ====================================================================================
  // SUPERNOVA PARAMETERS
  // ====================================================================================

  bool sn_injection;
  Real sn_energy;              //! SN injection energy, in ergs
  Real sn_ejecta_mass;         //! SN injection ejecta mass, in solar mass
  Real sn_inj_time;            //! SN injection time, code units
  Real sn_inj_dt;              //! SN injection interval, code units, negative means 1 SN
  Real sn_inj_radius;          //! SN injection radius, pc
  int sn_posn_type;            //! 0 = at center, 1 = random everytime

  Real sn_ejecta_Z = 0.1;          //! Metallicity of supernova ejecta
  Real sn_ejecta_dust_ratio = 0.2; //! Dust to metal ratio in supernova ejecta
  bool sn_done_flag = false;       //! If the it's a 1 SN run, and if restart after SN injected

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

static void SyncSNRestartState(const pgen_trml *p) {
  if (p == nullptr || p->pin == nullptr || !(p->sn_injection)) return;

  p->pin->SetReal("problem", "sn_inj_time", p->sn_inj_time);
  p->pin->SetBoolean("problem", "sn_done_flag", p->sn_done_flag);
  p->pin->SetInteger("problem", "sn_posn_type", p->sn_posn_type);
}

static void DumpInput(const pgen_trml *p) {
  if (global_variable::my_rank != 0) return;

  std::cout.precision(16);
  auto print = [&](const char *name, const auto &value) {
    std::cout << name << " = " << value << '\n';
  };

  print("ztop", p->ztop);
  print("zbot", p->zbot);
  print("gamma_adi", p->gamma_adi);
  print("dfloor", p->dfloor);
  print("pfloor", p->pfloor);
  print("tfloor", p->tfloor);
  print("rho_0", p->rho_0);
  print("cloud_flag", p->cloud_flag);

  print("cgm_cooling", p->cgm_cooling);
  if (p->cgm_cooling){
    print("hrate", p->hrate);
  }

  if (p->cloud_flag){
    print("chi", p->chi);
    print("LboxR",p->LboxR);
    print("transition_R",p->transition_R);
  }

  print("dust_model", p->dust_model);
  if (p->dust_model){
    print("Z_gas", p->Z_gas);
    print("D_Z_init", p->D_Z_init);
    print("dust_a", p->dust_a);
    print("dust_b", p->dust_b);
    print("dust_exp", p->dust_exp);
    print("Ndust_bins", p->Ndust_bins);
    print("rho_gr", p->rho_gr);
    print("Z_solar", p->Z_solar);
    print("vturb_est_ncell", p->vturb_est_ncell);
    print("coag_shatt_flag", p->coag_shatt_flag);
  }

  print("sn_injection", p->sn_injection);
  if (p->sn_injection){
    print("sn_energy", p->sn_energy);
    print("sn_ejecta_mass", p->sn_ejecta_mass);
    print("sn_inj_time", p->sn_inj_time);
    print("sn_inj_dt", p->sn_inj_dt);
    print("sn_inj_radius", p->sn_inj_radius);
    print("sn_ejecta_dust_ratio", p->sn_ejecta_dust_ratio);
    print("sn_ejecta_Z", p->sn_ejecta_Z);
    print("sn_posn_type", p->sn_posn_type);
  }

  print("T_cold", p->T_cold);
  print("T_hot", p->T_hot);
  print("T_floor", p->T_floor);
  print("ddens_threshold", p->ddens_threshold);

  return;
}


// Global pointer to the TRML parameters structure.
// Allocated in ProblemGenerator::UserProblem() and persists for the simulation lifetime.
// 
pgen_trml* input = new pgen_trml();

constexpr int kMaxDustBins = 32;  // Upper bound supported in device stack arrays

// For returning fragment size limits
struct Alims{
  Real amin;
  Real amax;
};

//! \brief Add user-defined source terms (cooling/heating) to conserved variables.
//! Called every RK substep via user_srcs_func.
void AddUserSrcs(Mesh *pm, const Real bdt);

//! \brief Apply dust model to each cell
//! Includes thermal sputtering, accretion, shattering and coagulation
void AddDustSource(Mesh *pm, const Real bdt);

//! \brief Apply SN injection
void InjectSN(Mesh *pm);

void TurbulentHistory(HistoryData *pdata, Mesh *pm);

void UserWorkInLoop(Mesh *pm);

void RefinementCondition(MeshBlockPack* pmbp);

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  input->pin = pin;

  // Enroll user functions 
  user_srcs_func = AddUserSrcs;
  user_hist_func = TurbulentHistory;
  user_ref_func = RefinementCondition;
  user_work_in_loop_func = UserWorkInLoop;

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
  input->gamma_adi         = eos.gamma;
  input->dfloor            = eos.dfloor;
  input->pfloor            = eos.pfloor;
  input->tfloor            = eos.tfloor;

  Real rho_0               = pin->GetReal("problem", "rho_0");
  input->rho_0             = rho_0;
  input->Z_gas             = pin->GetReal("problem", "Z_gas");  

  input->T_floor           = input->tfloor * KELVIN;
  input->T_hot             = pin->GetOrAddReal("problem", "T_hot", 1e6);
  input->T_cold            = pin->GetOrAddReal("problem", "T_cold", 1e4);
  Real T_hot               = input->T_hot;

  input->cgm_cooling        = pin->GetBoolean("hydro", "cgm_cooling");
  if (input->cgm_cooling){
    input->hrate              = pin->GetReal("hydro", "hrate");
  }

  input->cloud_flag        = pin->GetBoolean("problem", "cloud_flag");
  if (input->cloud_flag){
    input->chi               = input->T_hot/input->T_cold;
    input->LboxR              = pin->GetReal("problem", "LboxR");
    input->transition_R       = pin->GetReal("problem", "transition_R");
  }

  // get ztop and bottom
  input->ztop              = pin->GetReal("mesh", "x3max");
  input->zbot              = pin->GetReal("mesh", "x3min");


  // For dust model
  input->dust_model        = pin->GetBoolean("problem", "dust_model");

  if (input->dust_model){
    input-> D_Z_init         = pin->GetReal("problem", "D_Z_init");  
    input->dust_a            = pin->GetReal("problem", "dust_a"); // Min dust size
    input->dust_b            = pin->GetReal("problem", "dust_b"); // Max dust size
    input->dust_exp          = pin->GetReal("problem", "dust_exp"); // Dust size distribution exponent
    input->rho_gr            = pin->GetReal("problem", "rho_gr"); // Dust grain density

    input->vturb_est_ncell = pin->GetInteger("problem", "vturb_est_ncell"); // vturb estimation num of cells

    input->coag_shatt_flag = pin->GetBoolean("problem", "coag_shatt_flag");
  }

  // For SN injection
  input->sn_injection = pin->GetBoolean("problem", "sn_injection");

  if (input->sn_injection){
    input->sn_energy = pin->GetReal("problem", "sn_energy"); // in ergs
    input->sn_ejecta_mass = pin->GetReal("problem", "sn_ejecta_mass"); // in solar mass 

    input->sn_inj_time = pin->GetReal("problem", "sn_inj_time");
    input->sn_inj_dt = pin->GetReal("problem", "sn_inj_dt");
    input->sn_posn_type = pin->GetOrAddInteger("problem", "sn_posn_type", 0);

    input->sn_inj_radius = pin->GetReal("problem", "sn_inj_radius"); // in pc
    input->sn_ejecta_dust_ratio =
        pin->GetOrAddReal("problem", "sn_ejecta_dust_ratio", input->sn_ejecta_dust_ratio);
    input->sn_ejecta_Z = pin->GetOrAddReal("problem", "sn_ejecta_Z", input->sn_ejecta_Z);

    input->sn_done_flag = pin->GetOrAddBoolean("problem", "sn_done_flag", false);
    SyncSNRestartState(input);
  }

  input-> Z_solar          = 0.02;
  Real Z_gas               = input->Z_gas;
  Real D_Z_init            = input->D_Z_init;
  Real Z_solar             = input->Z_solar;


  // Number of dust bins: total scalars minus metal and cloud tracer entries.
  // Must be initialized even on restarts because source terms/history use it.
  const int Ndust_bins = nscalars - 2;
  if (Ndust_bins <= 0) {
    throw std::runtime_error("Ndust_bins <= 0. Check <hydro>/nscalars in the input (need metal + tracer + >=1 dust bin)");
  }
  input->Ndust_bins = Ndust_bins;

  // Create dust-bin midpoint array used by source terms.
  {
    const Real dust_a = input->dust_a;
    const Real dust_b = input->dust_b;

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
    dedges_ph(0) = dedges_ph(1) / r / 10.0;
    dedges_ph(Ndust_bins+2) = dedges_ph(Ndust_bins+1) * r * 10.0;
    for (int id = 0; id < Ndust_bins; ++id) {
      dbins_h(id) = 0.5*(dedges_h(id)+dedges_h(id+1));
    }
    
    Kokkos::deep_copy(dbins, dbins_h);
    Kokkos::deep_copy(dedges, dedges_h);
    Kokkos::deep_copy(dedges_pad, dedges_ph);

    input->dbins = dbins;
    input->dedges= dedges;
    input->dedges_pad= dedges_pad;
  }


  // Read the density gradient threshold for refinement
  input->ddens_threshold = pin->GetReal("problem", "ddens_max");

  // Print the input variables
  DumpInput(input);

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

  Real box = abs(input->ztop - input->zbot);
  Real radius = box/input->LboxR;
  Real smoothing_thickness = radius/input->transition_R;
  Real chi = input->chi;

  const bool cloud_flag = input->cloud_flag;

  int nmb1 = pmbp->nmb_thispack - 1;

  const int nmetal = nfluid;                 // index for metal
  const int ntracer = nfluid+1;              // index for tracer for cloud
  const int ndusti = nfluid+2;               // first dust bin num scalar

  printf("======================\n");
  printf("nmetal = %d, ntracer = %d, ndusti = %d\n", nmetal, ntracer, ndusti);
  printf("Ndust_bins = %d\n", Ndust_bins);
  printf("======================\n");
  printf("coag_shatt_flag: %d\n", static_cast<int>(input->coag_shatt_flag));
  printf("======================\n");

  // Set initial conditions
  if (global_variable::my_rank == 0) {
    std::cout << "Now initializing hydro/MHD variables." << "\n";
  }

  par_for("pgen_turb", DevExeSpace(),0,nmb1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {

    Real shape = 0.0;
    if (cloud_flag) {
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

      shape = 0.5 * (1.0 + std::tanh((radius - std::sqrt(x1v*x1v + x2v*x2v + x3v*x3v))
                                     / smoothing_thickness));
    }

    w0(m,IDN,k,j,i) = rho_0*(1.0 + (chi-1.0)*shape);
    w0(m,IVX,k,j,i) = 0.0;
    w0(m,IVY,k,j,i) = 0.0;
    w0(m,IVZ,k,j,i) = 0.0;
    if (eos.is_ideal) {
      w0(m,IEN,k,j,i) = T_hot * rho_0/ KELVIN  / gm1;
    }

    // add passive scalars
    if(nscalars>0){

      // First scalar has to be metallicity
      w0(m,nmetal,k,j,i) = Z_gas * Z_solar;

      // Cloud material tracer
      w0(m,ntracer,k,j,i) = 1.0 * shape;
      
      Real tmp_dust_shape = D_Z_init * Z_gas * Z_solar/ Real(Ndust_bins); //* shape / float(Ndust_bins);
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
  // if (input->adjust_temp_floor) AdjustTempTFloor(pm,bdt,u0,w0,eos_data);
  if (input->dust_model){
    AddDustSource(pm, bdt);
  }

  return;
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
Real mass_to_num(Real mass_in_bin, Real edges1, Real edges2, const Real rho_gr){
  // Convert mass change in bins to change in num distribution

    Real rho_gr_um = rho_gr * 1.0e-12;
    Real K_dust_um = (4.0/3.0)*M_PI*rho_gr_um;

    Real dist_change = 4.*mass_in_bin/K_dust_um;
    dist_change /= pow(edges2, 4.) - pow(edges1, 4.);

  return dist_change;
}

KOKKOS_INLINE_FUNCTION
Real dist_num_integral(const Real* d_arr, const Real* edges, const int n,
                       Real a1, Real a2, Real rho_gr, bool mass=false){
  // Integrate the distribution d_arr on edges between a1 and a2
  // These arrays are padded arrays

  int i1 = searchsorted(edges, n, a1)-1;
  int i2 = searchsorted(edges, n, a2)-1;

  Real rho_gr_um = rho_gr * 1.0e-12;
  Real K_dust_um = (4.0/3.0)*M_PI*rho_gr_um;

  Real integral = 0.0;   
  // Both are in the same bin
  if (i1==i2){

    if (mass) integral += 0.25 * d_arr[i1] * K_dust_um * (pow(a2, 4) - pow(a1, 4));
    else integral +=d_arr[i1]*(a2-a1);

    return integral;

  }

  // Add the left edge
  if (mass) integral += 0.25 * d_arr[i1] * K_dust_um * (pow(edges[i1+1], 4) - pow(a1, 4));
  else integral += d_arr[i1] * (edges[i1+1] - a1);
  

  // Add the right edge
  if (mass) integral += 0.25 * d_arr[i2] * K_dust_um * (pow(a2, 4) - pow(edges[i2], 4));
  else integral += d_arr[i2] * (a2 - edges[i2]);

  // Add the rest
  if (i2 > i1+1){

    if (mass) {
      for (int i=i1+1; i<i2; i++){
        integral += 0.25 * d_arr[i] * K_dust_um * (pow(edges[i+1], 4) - pow(edges[i], 4));
      }
    }
    else {
      for (int i=i1+1; i<i2; i++){
        integral += d_arr[i] * (edges[i+1] - edges[i]);
      }
    }

  }

  return integral;
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

KOKKOS_INLINE_FUNCTION
void rebin(Real* d_arr, const Real* edges, const int n, Real shift, Real rho_gr) {
  // We take in the arrays already padded with ghosts
  constexpr int kMaxRebinEntries = 35;  // Ndust_bins <= 32, plus three edge entries
  if (n > kMaxRebinEntries) return;

  Real mass_in_bin[kMaxRebinEntries];
  Real shifted_edges[kMaxRebinEntries];
  for (int i=0; i<n; i++){
    mass_in_bin[i] = 0.0;
    shifted_edges[i] = edges[i]+shift;
  }

  // Calculate mass in the unshifted bins
  for (int i=0; i<n-1; i++){
    mass_in_bin[i] = dist_num_integral(
      d_arr, shifted_edges, n, edges[i], edges[i+1], rho_gr, true);
  }

  // Calculate the new number dist
  for (int i=1; i<n-1; i++){
    d_arr[i] = mass_to_num(mass_in_bin[i],edges[i],edges[i+1], rho_gr);
  }
  d_arr[0] = mass_in_bin[0];
  d_arr[n-1] = mass_in_bin[n-1];

  return;
}

KOKKOS_INLINE_FUNCTION
Real v_grain(Real a, Real M_g, Real n_H, Real T, Real rho_gr){
  // a in um, M local Mach number, n_H in cm^-3
  // T in K, rho_gr in g/cc
  // From Eqn (18) Hirashita & Aoyama 2019
  return 1.1*pow(M_g, 1.5)*sqrt(a/0.1)*pow(T/1.0e4, 0.25)*pow(n_H, -0.25)*sqrt(rho_gr/3.5);
  // in km/s
}

KOKKOS_INLINE_FUNCTION
Real maxwell_tail_mean(Real v, Real vgr){
  constexpr Real vgr_tol = 1.0e-20;
  if (vgr<=vgr_tol) return 0.0;

  Real v0 = vgr * sqrt(2./3.);

  Real vmean = sqrt(8./M_PI) * v0;
  vmean *= 1. + 0.5*pow(v/v0, 2.);
  vmean *= exp(-0.5*pow(v/v0,2.));

  return vmean;
}

KOKKOS_INLINE_FUNCTION
Real maxwell_head_mean(Real v, Real vgr){
  constexpr Real vgr_tol = 1.0e-20;
  if (vgr<=vgr_tol) return 0.0;

  Real v0 = vgr * sqrt(2./3.);

  Real vmean = sqrt(8/M_PI) * v0;
  vmean *= 1. + 0.5*pow(v/v0, 2.);
  vmean *= exp(-0.5*pow(v/v0,2.));

  return (sqrt(8/M_PI)*v0 - vmean);
}

KOKKOS_INLINE_FUNCTION
void calc_interaction(Real bdt, Real intr_arr[][kMaxDustBins], const Real dist[],
                      const Real dbins[], const Real edges[], int nbin,
                      Real M_g, Real rho_gr, Real n_H, Real T, Real vol_cc, bool shatter) {
  // Assuming bdt is in seconds

  Real F_stick = 10.;

  Real gamma_Si = 2.7;
  Real gamma_C = 1.2;
  Real gamma_d = 0.5*(gamma_C + gamma_Si);

  // Young's modulus
  Real E_Si = 5.4e11;
  Real E_C = 3.4e10; 
  Real E_d = 0.5*(E_C + E_Si);

  // Shattering threshold velocity
  Real v_shatt_Si = 2.7;  // km/s
  Real v_shatt_C = 1.2;  // km/s
  Real v_shatt = 0.5 * (v_shatt_C + v_shatt_Si);  // km/s

  Real um_cgs = 1.0e-4;

  Real vc = 2.14;
  vc *= F_stick * pow(gamma_d, 5./6.) / pow(E_d, 1./3.) / sqrt(rho_gr);

  for (int i=0; i<nbin; i++){

    Real dela1i = edges[i+1] - edges[i];
    Real dela2i = pow(edges[i+1], 2.) - pow(edges[i], 2.);
    Real dela3i = pow(edges[i+1], 3.) - pow(edges[i], 3.);

    for (int j=0; j<nbin; j++){

      //* Can be cached if not memory-limited
      Real dela1j = edges[j+1] - edges[j];
      Real dela2j = pow(edges[j+1], 2.) - pow(edges[j], 2.);
      Real dela3j = pow(edges[j+1], 3.) - pow(edges[j], 3.);

      Real vproc = 0.0; 

      if (shatter){
        vproc = maxwell_tail_mean(v_shatt,
                                  v_grain(dbins[i], M_g, n_H, T, rho_gr)); // km/s
      }     
      else {
        vproc = (pow(dbins[i], 3.) + pow(dbins[j], 3.)) /
                pow(dbins[i] + dbins[j], 3.);
        vproc = sqrt(vproc)* pow((dbins[i] + dbins[j]) / dbins[i] / dbins[j], 5./6.);
        vproc = maxwell_head_mean(vproc,
                                  v_grain(dbins[i], M_g, n_H, T, rho_gr)); // km/s
      }

      intr_arr[i][j] = dela3i*dela1j/3. + dela2i*dela2j/2. +
                       dela1i*dela3j/3.; // um^4

      intr_arr[i][j] *= M_PI * vproc / vol_cc; // in km/s / cc
      intr_arr[i][j] *= dist[i] * dist[j]; // (#/um)^2

      // Result has a unit of um^2 km/s / cc
      intr_arr[i][j] *= 1.0e-3*bdt; // now dimensionless

    }
  }

  return;
}

KOKKOS_INLINE_FUNCTION
void add_coagulate(Real intr_arr[][kMaxDustBins], Real mass_in_bin[],
                   const Real dbins[], const Real edges[], int nbin,
                   Real rho_gr) {
  //! Assuming that intr_arr is already in the correct units
  //TODO: Check the units

  // Maybe pass this to avoid recalculating this?
  // Dust density in g/um^3
  const Real rho_gr_um = rho_gr * 1.0e-12;
  Real K_dust_um = (4.0/3.0)*M_PI*rho_gr_um;
  
  for (int i=0; i<nbin; i++){
    for (int j=0; j<nbin; j++){

          Real ai3 = pow(dbins[i], 3.);
          Real aj3 = pow(dbins[j], 3.);
          Real a_coag = pow(ai3+aj3, 1. / 3.);

          int k = searchsorted(edges, nbin + 1, a_coag)-1;

          if ((k < nbin) && (k>0)){ 
            // Remove coagulated grains
            mass_in_bin[i] -= intr_arr[i][j] * K_dust_um * ai3;
            mass_in_bin[j] -= intr_arr[i][j] * K_dust_um * aj3;
            // Add back coagulated product grain
            mass_in_bin[k] += intr_arr[i][j] * K_dust_um * (ai3 + aj3);
          }

    }
  }
  return;
}

KOKKOS_INLINE_FUNCTION
Real sigma_fn(Real M, const Real s){

  Real M_inv = 1./M;
  return 0.3 * pow(s + M_inv - 0.11, 0.13)/(s + M_inv - 1.);
}

// TODO: Come up with a way to abstract away the dust constants
KOKKOS_INLINE_FUNCTION
Real M_shocked(Real Mrel, Real Mproj, const Real rho_gr, const Real sigma_1,
               const Real s, const Real R, const Real M_1) {

  constexpr Real Mrel_tol = 1.0e-6;
  if (Mrel<=Mrel_tol) return 0.0;

  Real sigma_r = sigma_fn(Mrel/(1.+R), s);

  Real Msh = (1.+2.*R) / pow(1.+R, 9./16.) / pow(sigma_r, 1./9.);
  Msh *= pow(Mrel/sigma_1/M_1, 8./9.);

  return (Msh * Mproj); 
}

KOKKOS_INLINE_FUNCTION
Alims a_lim_frag(Real Mfrac, const Real P1, const Real rho_gr, const Real Pv, const Real z_const){
  // size limits for the fragments of the target

  // Dust density in g/um^3
  const Real rho_gr_um = rho_gr * 1.0e-12;

  Alims afrac_lim{0.0, 0.0};

  afrac_lim.amax = 3./4./M_PI/rho_gr_um * (1. + z_const);
  afrac_lim.amax /= 4. * pow(z_const, 3.) * (z_const - 2.);
  afrac_lim.amax *= Mfrac;

  afrac_lim.amin = afrac_lim.amax * pow(P1/Pv, 1.47);
   
  return afrac_lim;
}

KOKKOS_INLINE_FUNCTION
Real frag_norm(Real Mfrag, Alims alim, Real expo, const Real rho_gr) {
  // Normalization constant for the shattered fragment distribution

  const Real rho_gr_um = rho_gr * 1.0e-12;
  Real K_dust_um = (4.0/3.0)*M_PI*rho_gr_um;
  Real mex = expo + 4.;

  return Mfrag * mex / K_dust_um / (pow(alim.amax, mex) - pow(alim.amin, mex));
}

KOKKOS_INLINE_FUNCTION
Real integrate_frag(Real Anorm, Real a1, Real a2, Real expo, const Real rho_gr){
  // Integrate and return the fragment mass in a bin
  // amin and amax are sim dust range

  const Real rho_gr_um = rho_gr * 1.0e-12;
  Real K_dust_um = (4.0/3.0)*M_PI*rho_gr_um;
  Real mex = expo + 4.;

  return Anorm * K_dust_um * (pow(a2, mex) - pow(a1, mex)) / mex;
}

KOKKOS_INLINE_FUNCTION
void deposit_frag(Real mass_in_bin[], const Real edges[], int nbin, Alims alim,
                  Real Anorm, Real expo, const Real rho_gr) {
  // Deposit fragment mass in grain distribution

  int k = searchsorted(edges, nbin + 1, alim.amin)-1;
  int l = searchsorted(edges, nbin + 1, alim.amax)-1;

  if (l>=nbin){
    Kokkos::printf("### FATAL ERROR in %s at line %d\n", __FILE__, __LINE__);
    Kokkos::printf("Shattered fragments bigger than amax! "
                   "amin=%e amax=%e last_edge=%e nbin=%d\n",
                   static_cast<double>(alim.amin), static_cast<double>(alim.amax),
                   static_cast<double>(edges[nbin]), nbin);
    Kokkos::abort("deposit_frag overflowed dust-bin bounds");
  }

  if (l<0){ // Whole fragment distribution out of bounds
    mass_in_bin[0] += integrate_frag(Anorm, alim.amin, alim.amax, expo, rho_gr);
    return;
  }

  Real aleft = alim.amin;
  if (k<0){ // Just frag_amin out of bounds
    mass_in_bin[0] += integrate_frag(Anorm, alim.amin, edges[0], expo, rho_gr);
    k = 0;
    aleft = edges[0];
  }

  // if in the same bin
  if (k == l) {
    mass_in_bin[k] += integrate_frag(Anorm, aleft, alim.amax, expo, rho_gr);
    return;
  }

  // left edge
  mass_in_bin[k] += integrate_frag(Anorm, aleft, edges[k+1], expo, rho_gr);

  for (int i = k + 1; i < l; ++i) {
    mass_in_bin[i] += integrate_frag(Anorm, edges[i], edges[i+1], expo, rho_gr);
  }

  // right edge
  mass_in_bin[l] += integrate_frag(Anorm, edges[l], alim.amax, expo, rho_gr);

  return;
}

KOKKOS_INLINE_FUNCTION
void add_shatters(Real intr_arr[][kMaxDustBins], Real mass_in_bin[],
                  const Real dbins[], const Real edges[], int nbin, Real M_g,
                  const Real rho_gr, Real n_H, Real T) {
  
  // Dust property consts
  constexpr Real R = 1.0;

  // Critical Pressure
  constexpr Real P1_Si = 3e11;  // dyn/cm^2
  constexpr Real P1_C = 4e10;  // dyn/cm^2
  constexpr Real P1 = (P1_Si + P1_C) / 2.;  // dyn/cm^2

  // Critical Pressure of vapourisation
  constexpr Real Pv_Si = 5.4e12; // dyn/cm^2
  constexpr Real Pv_C = 5.8e12; // dyn/cm^2
  constexpr Real Pv = (Pv_Si + Pv_C) / 2.;  // dyn/cm^2

  constexpr Real v_shatt_Si = 2.7;  // km/s
  constexpr Real v_shatt_C = 1.2;  // km/s
  constexpr Real v_shatt = 0.5 * (v_shatt_C + v_shatt_Si);  // km/s

  // Dust property const
  constexpr Real s_Si = 1.2;
  constexpr Real s_C = 1.9;
  constexpr Real s = (s_Si + s_C) / 2.;  // dimensionless

  // Const in reln for excavation flow
  constexpr Real z_const = 3.4;  // dimensionless

  // Grain sound speed
  // Tielens+ 1994, Table 1
  // https://ui.adsabs.harvard.edu/abs/1994ApJ...431..321T/abstract
  constexpr Real c0_Si = 5.0;  // km/s
  constexpr Real c0_C = 1.8;  // km/s
  constexpr Real c0 = (c0_Si + c0_C) / 2.;  // km/s
  constexpr Real c0_cgs = c0 * 1.0e5; // cm/s

  constexpr Real shatter_expo = -3.3;
  constexpr Real mex = shatter_expo + 4.;

  const Real phi_1 = P1 / rho_gr/ c0_cgs/ c0_cgs;
  const Real M_1 = 2. * phi_1 / (1. + sqrt(1. + 4.*s*phi_1));
  const Real sigma_1 = sigma_fn(M_1, s);

  const Real rho_gr_um = rho_gr * 1.0e-12;
  const Real K_dust_um = (4.0/3.0)*M_PI*rho_gr_um;

  constexpr Real Mfrag_tol = 1.0e-20;

  for (int i=0; i<nbin; i++){
    for (int j=0; j<nbin; j++){

      Real ai3 = pow(dbins[i], 3.);
      Real aj3 = pow(dbins[j], 3.);
      Real Mi = K_dust_um * ai3;
      Real Mj = K_dust_um * aj3;

      // Real Mproj = Kokkos::min(Mi, Mj);
      // Real Mtarg = Kokkos::max(Mi, Mj);

      Real v_grain_i = v_grain(dbins[i], M_g, n_H, T, rho_gr);
      Real vrel = maxwell_tail_mean(v_shatt, v_grain_i);
      Real Msh = M_shocked(vrel/c0, Mj, rho_gr, sigma_1, s, R, M_1);

      Real Mfrag = (Msh > 0.5*Mi) ? Mi : 0.4 * Msh;
      Real Mleft = Mi - Mfrag;

      if ((Mfrag/Mi)>= Mfrag_tol){
        Real a_left = pow(Mleft/K_dust_um, 1./3.);
        
        // Take care of parent grains, and left over target (left)
        int k = searchsorted(edges, nbin + 1, a_left)-1;
        if ((k < nbin) && (k>0)){ 
          // Remove shattered grains
          mass_in_bin[i] -= intr_arr[i][j] * Mi;
          mass_in_bin[j] -= intr_arr[i][j] * Mj;
          // Add back left product grain twice
          // Once for (i,j) and once of (j,i)
          mass_in_bin[k] += 2.*intr_arr[i][j] * Mleft;
        }

        //* Take care of fragments from M1
        // Factor of 2 for Mfrag, is to include (j,i) in this calculation
        Alims frag_alim = a_lim_frag(2.*Mfrag, P1, rho_gr, Pv, z_const);
        Real frag_norm =  Mfrag * mex / K_dust_um / (pow(frag_alim.amax, mex) - pow(frag_alim.amin, mex));

        deposit_frag(mass_in_bin, edges, nbin, frag_alim, frag_norm, shatter_expo, rho_gr);
      }

      // // Take care of fragments of the projectile
      // Alims proj_alim = a_lim_proj(Mproj, Mtarg, vrel, rho_gr, sigma_1, s, R,
      //                              c0, M_1);
      // Real proj_norm = frag_norm(Mproj, proj_alim, shatter_expo, rho_gr);

      // deposit_frag(mass_in_bin, edges, nbin, proj_alim, proj_norm, shatter_expo, rho_gr);

    }
  }


  return;
}

//! \fn void AddDustSource()
//! \brief Apply dust model to each cell
//! Includes thermal sputtering, accretion, shattering and coagulation
void AddDustSource(Mesh *pm, const Real bdt){
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
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
  Real MYR = 1.0; // how many code times in 1 Myr
  Real SEC = 1.0; // how many code units in a sec
  Real KM_S = 1.0; // how many code vel in 1 km/s
  Real CM3 = 1.0; // how many code volumes in 1 cc
  if (pmbp->punit != nullptr) {
    KELVIN = pmbp->punit->temperature_cgs();
    MYR = pmbp->punit->myr();
    SEC = pmbp->punit->s();
    CM3 = pow(pmbp->punit->cm(), 3.);
  } else if (global_variable::my_rank == 0) {
    std::cout << "ERROR: <units> block missing..." << std::endl;
    std::exit(1);
  }

  int nmetal = nfluid; // index for metal
  // int ntracer = nfluid+1; // index for tracer for cloud
  int ndusti = nfluid+2; // first dust bin scalar
  const int Ndust_bins = input->Ndust_bins; // Number of dust bins
  if (Ndust_bins > kMaxDustBins) {
    if (global_variable::my_rank == 0) {
      std::cout << "ERROR: Ndust_bins=" << Ndust_bins
                << " exceeds kMaxDustBins=" << kMaxDustBins << std::endl;
    }
    std::exit(1);
  }

  Real Z_solar = input->Z_solar;
  auto dust_bins = input->dbins;

  const Real gamma_adi = input->gamma_adi;
  const int vturb_est_ncell = input->vturb_est_ncell;

  const Real grain_porosity = input->grain_porosity;
  const Real rho_gr = input->rho_gr;
  const Real dust_a = input->dust_a;
  const Real dust_b = input->dust_b;
  const int n_edges = Ndust_bins + 3;

  const bool coag_shatt_flag = input->coag_shatt_flag;

  auto dedges_pad_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                           input->dedges_pad);
  auto dbins_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                           input->dbins);
  Real edges_pad[kMaxDustBins+3];
  Real edges[kMaxDustBins+1];
  Real dbins[kMaxDustBins];
  for (int i=0; i<n_edges; i++) edges_pad[i] = dedges_pad_h(i);
  for (int i=0; i<=Ndust_bins; i++) edges[i] = dedges_pad_h(i+1);
  for (int i=0; i<Ndust_bins; i++) dbins[i] = dbins_h(i);

  par_for("user_source", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

    Real dens = w0(m, IDN, k, j, i); // in cm^-3
    Real temp = (w0(m,IEN,k,j,i) * gm1) / dens * KELVIN;
    const Real M = 0.0;  // Placeholder until a local Mach estimate is implemented.

    // cell volume in cm^3
    Real vol_cc = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3 / CM3;

    // To prevent NaNs
    constexpr Real Z_tol = 1e-20; 
    constexpr Real D_tol = 1e-20; 
    Real clamp_dust_tot = 1.0; // tmp var used for clamping later

    Real Z_local = w0(m,nmetal,k,j,i);
    Z_local = Kokkos::max(Z_tol, Z_local);

    int kk = k - ks;
    int jj = j - js;
    int ii = i - is;
    int cell_idx = m*nkji + kk*nji + jj*nx1 + ii;

    Real dust_arr[kMaxDustBins+2];

    //! Should I use w0 scalar instead?
    // u0 is equiv to rho_d
    //*========== Initialisation =============//

    // Total dust ratio across size bins
    dust_arr[0] = 0.0;
    dust_arr[Ndust_bins+2] = 0.0; 
    // two additional ghost bins to catch overflow and underflow

    // For tolerance check 
    Real D_check = 0.0;

    for (int id=1; id<(Ndust_bins+1); id++){
      dust_arr[id] = mass_to_num(u0(m, ndusti+id-1, k, j, i), 
                        edges_pad[id], edges_pad[id+1], rho_gr);
      D_check += u0(m, ndusti+id-1, k, j, i);
    }

    // If too little dust, skip the dust evolution
    if (D_check <= D_tol) return;

    // Initial dust mass
    Real dust_tot = dist_num_integral(
      dust_arr, edges_pad, n_edges, dust_a, dust_b, rho_gr, true);

    Real delta_dmass = -dust_tot;

    //*========== Sputtering & Accretion =============//

    // Calculate shift from sputtering and accretion
    Real tot_shift = sputtering(dens, temp, Z_local/Z_solar, grain_porosity);
    tot_shift += accretion(dens, temp, Z_local/Z_solar, grain_porosity); // in um/Myr
    tot_shift *= bdt/MYR; // bdt converted to Myr


    // Rebin after shifting
    rebin(dust_arr, edges_pad, n_edges, tot_shift, rho_gr);

    //*=========== Turbulent velocity estimation ==============//
    // Use Welford's online method to calculate the std dev of velocity

    Real vx_mean = 0.0, vy_mean = 0.0, vz_mean = 0.0;
    Real M2x = 0.0, M2y = 0.0, M2z = 0.0;
    int n_count = 0;

    for (int ia=-vturb_est_ncell; ia<=vturb_est_ncell; ia++){
      for (int ja=-vturb_est_ncell; ja<=vturb_est_ncell; ja++){
        for (int ka=-vturb_est_ncell; ka<=vturb_est_ncell; ka++){

          Real vx = w0(m, IVX, k+ka, j+ja, i+ia); 
          Real vy = w0(m, IVY, k+ka, j+ja, i+ia); 
          Real vz = w0(m, IVZ, k+ka, j+ja, i+ia); 

          n_count++;

          Real deltax = vx - vx_mean;
          Real deltay = vy - vy_mean;
          Real deltaz = vz - vz_mean;

          vx_mean += deltax / n_count;
          vy_mean += deltay / n_count;
          vz_mean += deltaz / n_count;

          Real delta2x = vx - vx_mean;
          Real delta2y = vy - vy_mean;
          Real delta2z = vz - vz_mean;

          M2x += deltax * delta2x;
          M2y += deltay * delta2y;
          M2z += deltaz * delta2z;

        }
      }
    }

    Real vturb_g = sqrt((M2x+M2y+M2z)/Real(n_count+1)); // at nbd, in code units
    
    // Assuming Kolmogorov, scale from nbd to cell lengthscale
    vturb_g *= pow(1/Real(n_count), 1./9.) / KM_S ; // in code units

    // Gas mach number
    Real cs_g = sqrt(gamma_adi*temp/KELVIN) * 1.0e-5; 
    Real M_g = vturb_g / cs_g;

    //*========== Coagulation & Shattering =============//

    // No coagulation or shattering without any changes if turbulent velocity is almost zero
    constexpr Real M_g_tol = 1e-6; // km/s
    if ((M_g >= M_g_tol) && coag_shatt_flag) { 

      Real intr_arr[kMaxDustBins][kMaxDustBins];  // Interaction array

      //* Coagulation
      // Calculate the interaction frequency
      calc_interaction(bdt/SEC, intr_arr, dust_arr + 1, dbins, edges, Ndust_bins, M_g,
                      rho_gr, dens, temp, vol_cc, false);

      // Calculate mass redistribution as mass_in_bin
      Real mass_in_bin[kMaxDustBins];   // Mass change
      for (int id=0; id<Ndust_bins; ++id) mass_in_bin[id] = 0.0;

      add_coagulate(intr_arr, mass_in_bin, dbins, edges, Ndust_bins, rho_gr);

      //* Shattering
      // Calculate the interaction frequency
      calc_interaction(bdt/SEC, intr_arr, dust_arr + 1, dbins, edges, Ndust_bins, M_g,
                      rho_gr, dens, temp, vol_cc, true);

      // Calculate mass redistribution as mass_in_bin
      add_shatters(intr_arr, mass_in_bin, dbins, edges, Ndust_bins, M_g,
                  rho_gr, dens, temp);

      // Add the changes to dust_arr from coagulation and shattering
      for (int id=1; id<(Ndust_bins+1); ++id){ 
        dust_arr[id] += mass_to_num(mass_in_bin[id-1], edges[id-1], edges[id], rho_gr);
      }
    }

    //*============= Reassignment =================//

    // Mass after rebinning
    dust_tot = dist_num_integral(
      dust_arr, edges_pad, n_edges, dust_a, dust_b, rho_gr, true);
    delta_dmass += dust_tot;

    // Kokkos::printf("delta_dmass: %.50lf um \n", delta_dmass);

    // Calculate metal mass change
    Real rate_Z = -delta_dmass;
    rate_Z *= 1.0;  //! Assuming Z_dust ~ 1

    // Kokkos::printf("rate_Z: %.50lf um \n", rate_Z);

    const Real rho_gr_um = rho_gr * 1.0e-12;
    const Real K_dust_um = (4.0/3.0)*M_PI*rho_gr_um;

    // Loop through the dust bins
    for (int id=0; id<Ndust_bins; id++){
      u0(m, ndusti+id, k, j, i) = 0.25 * K_dust_um * dust_arr[id+1];
      u0(m, ndusti+id, k, j, i) *= pow(edges_pad[id+2], 4.)-pow(edges_pad[id+1], 4.);
    }

    //! These are w0 * dens
    u0(m, nmetal, k, j, i) += rate_Z;

    // We should check that scalars are guaranteed to be in [0,1] after all source terms are added.
    // Clamp metallicity 
    u0(m, nmetal, k, j, i) = Kokkos::clamp(u0(m, nmetal, k, j, i), 0.0, dens);
    // Clamp total dust-to-gas ratio to [0,1]
    // Dust to metal ratio is not clamped, as all of gas metals can be locked in dust
    // Also dust might have been created in a metal-rich area and moved...
    clamp_dust_tot = Kokkos::clamp(dust_tot, 0.0, dens);

    // Scaled to sum to clamp_dust_tot 
    for (int id=0; id<Ndust_bins; id++){
      if (dust_tot <= D_tol) u0(m, ndusti+id, k, j, i) = 0.0;
      else u0(m, ndusti+id, k, j, i) *= clamp_dust_tot/dust_tot;
    }

  });

  return;
}

//! \fn void AddSN()
//! \brief Apply SN injection
void InjectSN(Mesh *pm){

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;
  bool is_mhd = (pmbp->pmhd != nullptr) ? true : false;
  auto &u0 = (is_mhd) ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  int &nfluid = (is_mhd) ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;
  int &nscalars = (is_mhd) ? pmbp->pmhd->nscalars : pmbp->phydro->nscalars;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;


  // Get code unit variables for temperature conversions
  Real PC = 1.0; // how many code volumes in 1 pc
  Real MSOL = 1.0; // how many code mass in 1 solar mass
  Real ERG = 1.0; // how many code energies in 1 erg
  if (pmbp->punit != nullptr) {
    PC = pmbp->punit->pc();
    MSOL = pmbp->punit->msun();
    ERG = pmbp->punit->erg();
  } else if (global_variable::my_rank == 0) {
    std::cout << "ERROR: <units> block missing..." << std::endl;
    std::exit(1);
  }

  const int nmetal = nfluid; // index for metal
  const int ntracer = nfluid+1; // index for tracer for cloud
  const int ndusti = nfluid+2; // first dust bin scalar
  const int Ndust_bins = input->Ndust_bins; // Number of dust bins

  const Real sn_inj_radius = input->sn_inj_radius*PC; // code length
  const Real sn_inj_vol = 4./3. * M_PI * pow(sn_inj_radius, 3.);

  const Real sn_energy_den = input->sn_energy*ERG/sn_inj_vol; // code energy
  const Real sn_ejecta_den= input->sn_ejecta_mass*MSOL/sn_inj_vol; // code mass

  const Real sn_ejecta_Z = input->sn_ejecta_Z;
  const Real sn_ejecta_dust_ratio = input->sn_ejecta_dust_ratio;

  static Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>
      sn_rand_pool(1234567);

  Real sn_x = 0.5*(pm->mesh_size.x1min + pm->mesh_size.x1max);
  Real sn_y = 0.5*(pm->mesh_size.x2min + pm->mesh_size.x2max);
  Real sn_z = 0.5*(pm->mesh_size.x3min + pm->mesh_size.x3max);

  if (input->sn_posn_type == 1){
    // Generate a random SN injection position
    if (global_variable::my_rank == 0) {
      auto rand_gen = sn_rand_pool.get_state();

      Real x1min = pm->mesh_size.x1min + sn_inj_radius;
      Real x1max = pm->mesh_size.x1max - sn_inj_radius;
      Real x2min = pm->mesh_size.x2min + sn_inj_radius;
      Real x2max = pm->mesh_size.x2max - sn_inj_radius;
      Real x3min = pm->mesh_size.x3min + sn_inj_radius;
      Real x3max = pm->mesh_size.x3max - sn_inj_radius;

      if (x1min > x1max) {
        Real xmid = 0.5*(pm->mesh_size.x1min + pm->mesh_size.x1max);
        x1min = xmid;
        x1max = xmid;
      }
      if (x2min > x2max) {
        Real ymid = 0.5*(pm->mesh_size.x2min + pm->mesh_size.x2max);
        x2min = ymid;
        x2max = ymid;
      }
      if (x3min > x3max) {
        Real zmid = 0.5*(pm->mesh_size.x3min + pm->mesh_size.x3max);
        x3min = zmid;
        x3max = zmid;
      }

      Real rand = rand_gen.frand();
      sn_x = (1.0 - rand)*x1min + rand*x1max;
      rand = rand_gen.frand();
      sn_y = (1.0 - rand)*x2min + rand*x2max;
      rand = rand_gen.frand();
      sn_z = (1.0 - rand)*x3min + rand*x3max;

      sn_rand_pool.free_state(rand_gen);
    }
  }

// Communicate the newly generated position to all MPI ranks
#if MPI_PARALLEL_ENABLED
  MPI_Bcast(&sn_x, 1, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sn_y, 1, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sn_z, 1, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
#endif


  par_for("user_source", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

    Real x1min = size.d_view(m).x1min;
    Real x1max = size.d_view(m).x1max;

    Real x2min = size.d_view(m).x2min;
    Real x2max = size.d_view(m).x2max;

    Real x3min = size.d_view(m).x3min;
    Real x3max = size.d_view(m).x3max;

    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    Real dx = x1v - sn_x;
    Real dy = x2v - sn_y;
    Real dz = x3v - sn_z;

    if ((dx*dx + dy*dy + dz*dz) < sn_inj_radius*sn_inj_radius){
      u0(m, IDN, k, j, i) += sn_ejecta_den; 
      u0(m, IEN, k, j, i) += sn_energy_den;


      // add passive scalars
      if(nscalars>0){
        // First scalar has to be metallicity
        u0(m,nmetal,k,j,i) += sn_ejecta_Z * sn_ejecta_den * (1.-sn_ejecta_dust_ratio);

        // Cloud material tracer
        u0(m,ntracer,k,j,i) += 1.0 * sn_ejecta_den;
        
        Real tmp_dust_shape = sn_ejecta_Z * sn_ejecta_den * sn_ejecta_dust_ratio
                              / Real(Ndust_bins);
        for (int id=0; id<Ndust_bins; id++){
          // Dust
          u0(m,ndusti+id,k,j,i) += tmp_dust_shape;
        }
      }

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
  const int Ndust_bins = input->Ndust_bins; // Number of dust bins

  int count = 0;
  pdata->label[count] = "U^2"; count++;
  pdata->label[count] = "Mcold"; count++;
  pdata->label[count] = "T_sum"; count++;

  bool dust_model = input->dust_model;

  int DUST_HIST = count;
  if (dust_model){
    pdata->label[count] = "Zgas"; count++;
    for (int id = 0; id < input->Ndust_bins; ++id) {
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
  Real T_cold = input->T_cold;

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

// ! \fn void UserWorkInLoop()
// ! \brief Function called in hydro or mhd tasks in "after_timeintegrator" stage
void UserWorkInLoop(Mesh *pm) {
  // Frame tracking


  if (!input->sn_injection || input->sn_done_flag) return;
  if ((pm->time < input->sn_inj_time) || input->sn_done_flag) return;

  //! If the simulation with (sn_inj_dt<0, i.e. single SN) is restarted
  //! after the SN has been injected, set sn_done_flag in input to 0 !! 
  if (input->sn_inj_dt < 0.) input->sn_done_flag = true;


  InjectSN(pm);  // Inject supernova
  input->sn_inj_time += input->sn_inj_dt;
  SyncSNRestartState(input);

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

  auto &ddens_thresh = input->ddens_threshold;

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

// ----------------------------------------------------------------------------------------
