//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mass_removal_test.cpp
//  \brief Problem generator for testing mass removal
#include <iostream> // cout

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
  // Real velocity;               //!< Shear velocity amplitude (hot phase moves at +v, cold at -v)
  // Real t_shear;                //!< Shear timescale = 1/velocity (eddy turnover time)

  // ====================================================================================
  // INITIAL CLOUD PARAMETERS
  // ====================================================================================
  Real LboxR;                  //!< L_box/R_cloud ratio
  Real transition_R;           //!< R_cloud/transition layer thickness 

  // // ====================================================================================
  // // COOLING FUNCTION PARAMETERS
  // // ====================================================================================
  // Real xi;                     //!< Dimensionless ratio t_shear/t_cool (key parameter!)
  //                              //!< Controls whether cooling or mixing dominates
  // Real t_cool_start;           //!< Simulation time at which to turn on radiative cooling
  // Real beta_lo;                //!< Power-law slope of cooling function below T_peak
  // Real beta_hi;                //!< Power-law slope of cooling function above T_peak
  // Real alpha_heat;             //!< Derived heating coefficient exponent
  // Real heat_coefficient;       //!< Derived heating amplitude coefficient
  // bool balanced_heating;       //!< Use legacy heating normalization when true
  // Real dt_cutoff;              //!< Minimum allowed timestep (prevents runaway)
  // Real cfl_cool;               //!< Cooling CFL number (limits dT/T per timestep)

  // ====================================================================================
  // DUST MODEL PARAMETERS
  // ====================================================================================

  bool dust_model;             //! Flag to turn the dust model on/off 
  Real Z_gas;                  //! Initial Gas metallicity in Z_solar
  Real D_Z_init;               //! Initial dust-to-gas ratio
  Real Z_solar;                //! Solar metallicity

  // ====================================================================================
  // TEMPERATURE THRESHOLDS
  // ====================================================================================
  Real T_cold;                 //!< Cold phase temperature (= pgas_0/rho_0/contrast)
  Real T_hot;                  //!< Hot phase temperature (= pgas_0/rho_0)
  Real T_floor;                //!< Temperature floor for adjusting cold gas

  // //!< Temperature bands for phase classification (logarithmically spaced)
  // Real T_peak_lo;              //!< Lower bound of T_peak band (for diagnostics)
  // Real T_peak_hi;              //!< Upper bound of T_peak band (for diagnostics)
  // Real T_hot_lo;               //!< Lower temperature considered "hot"
  // Real T_cold_hi;              //!< Upper temperature considered "cold"

  // // ====================================================================================
  // // MHD PARAMETERS
  // // ====================================================================================
  // Real alpha_magdens;          //!< Power-law index for B-field scaling with density
  //                              //!< (used when B_dens_ratio = true)

  // ====================================================================================
  // ADAPTIVE MESH REFINEMENT THRESHOLDS
  // ====================================================================================
  Real ddens_threshold; 

  // Real t_start_refine;          //!< Start applying refinement only after this time
  // Real density_ratio_threshold; //!< Refine if density gradient exceeds this
  // Real vel2_rms_threshold;      //!< Refine if velocity RMS exceeds this
  // Real T_max_threshold;         //!< Refine if T_max in block exceeds this
  // Real T_min_threshold;         //!< Refine if T_min in block is below this
  // std::vector<Real> refine_level_times; //!< Time-based max refinement schedule (t:num_levels)
  // std::vector<int> refine_level_counts; //!< Num levels allowed at each scheduled time
  // int max_num_levels_input;      //!< Max num_levels from input file
  // int max_level_floor;           //!< Lowest allowed logical max level (from existing mesh)

  // // ====================================================================================
  // // COOLING CONTROL FLAGS AND PARAMETERS
  // // ====================================================================================
  // int  ndiag;                  //!< Output diagnostics every ndiag timesteps (-1 = off)
  // bool heating_on;             //!< Enable background heating
  // bool cooling_below_T_cold;   //!< Enable cooling for T < T_cold
  // bool cool_subc;              //!< Enable cooling subcycling (multiple cooling steps per hydro step)
  // Real cool_subfac;            //!< Subcycling factor (dt_cool = cool_subfac * t_cool)
  // bool use_temp_floor_cool;    //!< Enable temperature floor for cooling
  // Real T_floor_cool;           //!< Minimum temperature for cooling (cooling off below this)
  // bool use_temp_ceiling;       //!< Enable maximum temperature cap
  // Real T_ceiling;              //!< Maximum allowed temperature in simulation
  // bool adjust_temp_floor;      //!< If T < T_floor, average with neighbors
  // bool use_dens_ceiling;       //!< Enable density ceiling for cooling
  // Real dens_ceiling;           //!< Disable cooling above this density

  // // ====================================================================================
  // // DIAGNOSTIC TRACKING VARIABLES
  // // ====================================================================================
  // Real tot_mass = 0.0;         //!< Total mass in domain (for conservation tracking)
  // Real tot_coolrate = 0.0;     //!< Total cooling rate (energy loss per unit time)

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
void AddDustSource(Mesh *pm, const Real bdt);

void TurbulentHistory(HistoryData *pdata, Mesh *pm);

void RefinementCondition(MeshBlockPack* pmbp);

//----------------------------------------------------------------------------------------
//  \brief Problem Generator for mass removal

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
  ptrml-> Z_solar          = 0.0134;  
  Real Z_gas               = ptrml->Z_gas;
  Real D_Z_init            = ptrml->D_Z_init;
  Real Z_solar             = ptrml->Z_solar;

  ptrml->LboxR              = pin->GetReal("problem", "LboxR");
  ptrml->transition_R       = pin->GetReal("problem", "transition_R");


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
  int &nscalars = (is_mhd) ?
                pmbp->pmhd->nscalars : pmbp->phydro->nscalars;
  int &nfluid = (is_mhd) ? pmbp->pmhd->nmhd : pmbp->phydro->nhydro;
  Real gm1 = eos.gamma - 1.0;
  auto &size = pmbp->pmb->mb_size;

  Real box = abs(ptrml->ztop - ptrml->zbot);
  Real radius = box/ptrml->LboxR;
  Real smoothing_thickness = radius/ptrml->transition_R;
  Real chi = ptrml->chi;

  int nmb1 = pmbp->nmb_thispack - 1;


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
    if(nscalars>0){
      w0(m,nfluid,k,j,i) = 1.0 * shape;

      // Dust
      w0(m,nfluid+1,k,j,i) = Z_gas * Z_solar * (0.1 + 0.9 * shape);
      w0(m,nfluid+2,k,j,i) = 0.5 * D_Z_init * Z_gas * Z_solar * shape;
      w0(m,nfluid+3,k,j,i) = 0.5 * D_Z_init * Z_gas * Z_solar * shape;
      // The D_tot comes out to be D_Z_init * Z_gas * Z_solar, i.e D_Z_init * Z_g
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

  Real Z_solar             = ptrml->Z_solar;

  Real size_d1 = 1.0; // in nm
  Real size_d2 = 10.0; // in nm

  Real sig_DL = 10; // in km/s
  Real sig_DS = 0.1; // in km/s
  Real F_coag = 0.5; // fudge factor for coagulation

  Real s_i = 3.0; // Dust grain density


  par_for("user_source", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dens = w0(m, IDN, k, j, i); // in cm^-3
    Real temp = (w0(m,IEN,k,j,i) * gm1) / dens * KELVIN;
    
    Real Z_local = w0(m, nfluid+1 , k, j, i)/Z_solar; // in solar metallicity
    Real dust_1 = w0(m, nfluid+2 , k, j, i);
    Real dust_2 = w0(m, nfluid+3, k, j, i);

    // To prevent NaNs
    Real Z_tol = 1e-20;
    Z_local = Kokkos::max(Z_tol, Z_local);
    dust_1 = Kokkos::max(Z_tol, dust_1);
    dust_2 = Kokkos::max(Z_tol, dust_2);

    Real dust_tot = dust_1 + dust_2;
    Real clamp_dust_tot = 1.0; // tmp var used for clamping later

    Real rho_d1 = dust_1 * dens;
    Real rho_d2 = dust_2 * dens;
    Real rate_d1 = 0.0;
    Real rate_d2 = 0.0;
    Real rate_Z = 0.0;
    Real tmp_rate = 0.0;

    
    //* Thermal sputtering
    Real t_sp = 70.0; // in Myr
    t_sp *= (1.e-3/dens);
    t_sp *= 1. + pow((temp/2.0e6), -2.5);


    tmp_rate = -rho_d1 / (t_sp * size_d1);
    rate_d1 += tmp_rate;
    rate_Z += -tmp_rate;


    tmp_rate = -rho_d2 / (t_sp * size_d2);
    rate_d2 += tmp_rate;
    rate_Z += -tmp_rate;


    //* Accretion
    Real t_ac = 200.0; // in Myr
    t_ac *= (20.0/dens);
    t_ac *= pow((50.0/temp), 0.5);
    t_ac *= Z_local; // in Z_sol 

    t_ac /= 1-(dust_tot/Z_local/Z_solar);


    tmp_rate = rho_d1 / (t_ac * size_d1);
    rate_d1 += tmp_rate;
    rate_Z += -tmp_rate;


    tmp_rate = rho_d2 / (t_ac * size_d2);
    rate_d2 += tmp_rate;
    rate_Z += -tmp_rate;


    // Convert rate_Z to metal mass
    rate_Z *= 1.0;  //! Assuming Z_dust ~ 1

    //* Shattering
    // From Dubois 2024
    Real t_shatt = 54.0; // in Myr
    t_shatt *= (1.0/dens); // in cm^-3
    t_shatt *= (s_i/3.0);  // in cm^-3
    t_shatt *= (0.01/dust_2); 
    t_shatt *= (10/sig_DL);  // in km/s

    tmp_rate = rho_d2 / (t_shatt * size_d2);

    rate_d1 += tmp_rate;
    rate_d2 += -tmp_rate;


    //* Coagulation
    // From Dubois 2024
    Real t_coag = 0.27; // in Myr
    t_coag *= (s_i/3.0);  // in cm^-3
    t_coag *= (1.0e3/dens); // in cm^-3
    t_coag *= (0.01/dust_2); 
    t_coag *= (0.1/sig_DS);  // in km/s
    t_coag *= F_coag;

    tmp_rate = rho_d1 / (t_coag * size_d1 / 0.05);

    rate_d1 += -tmp_rate;
    rate_d2 += tmp_rate;


    //! These are w0 * dens
    u0(m, nfluid+1, k, j, i) += bdt * rate_Z;
    u0(m, nfluid+2, k, j, i) += bdt * rate_d1;
    u0(m, nfluid+3, k, j, i) += bdt * rate_d2;

    // We should check that scalars are guaranteed to be in [0,1] after all source terms are added.
    // Clamp metallicity 
    u0(m, nfluid+1, k, j, i) = Kokkos::clamp(u0(m, nfluid+1, k, j, i), 0.0, dens);

    // Clamp total dust-to-gas ratio to [0,1]
    // Dust to metal ratio is not clamped, as all of gas metals can be locked in dust
    // Also dust might have been created in a metal-rich area and moved...
    dust_tot = u0(m, nfluid+2, k, j, i) + u0(m, nfluid+3, k, j, i);
    clamp_dust_tot = Kokkos::clamp(dust_tot, 0.0, dens);

    // Scaled to sum to clamp_dust_tot 
    u0(m, nfluid+2, k, j, i) *= clamp_dust_tot/dust_tot;
    u0(m, nfluid+3, k, j, i) *= clamp_dust_tot/dust_tot;


  });

  return;
}

void TurbulentHistory(HistoryData *pdata, Mesh *pm) {

  int count = 0;
  pdata->label[count] = "U^2"; count++;
  pdata->label[count] = "Mcold"; count++;
  pdata->label[count] = "T_sum"; count++;

  bool dust_model = ptrml->dust_model;

  int DUST_HIST = count;
  if (dust_model){
    pdata->label[count] = "Zgas"; count++;
    pdata->label[count] = "Dstot"; count++;
    pdata->label[count] = "Dltot"; count++;
  }

  pdata->nhist = count;

  auto &w0_ = pm->pmb_pack->phydro->w0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  int &nhist_ = pdata->nhist;

  bool is_mhd = (pm->pmb_pack->pmhd != nullptr) ? true : false;
  const EOS_Data &eos_data = (is_mhd) ?
                  pm->pmb_pack->pmhd->peos->eos_data : pm->pmb_pack->phydro->peos->eos_data;
  int &nfluid_ = (is_mhd) ? pm->pmb_pack->pmhd->nmhd : pm->pmb_pack->phydro->nhydro;

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
      hvars.the_array[DUST_HIST] += dens * w0_(m, nfluid_+1, k, j, i);
      hvars.the_array[DUST_HIST+1] += dens * w0_(m, nfluid_+2, k, j, i);
      hvars.the_array[DUST_HIST+2] += dens * w0_(m, nfluid_+3, k, j, i);
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
  auto &multi_d     = pmesh->multi_d;
  auto &three_d     = pmesh->three_d;
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

