//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mass_removal_test.cpp
//  \brief Problem generator for a convective box
#include <iostream> // cout

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"

void UserSource(Mesh* pm, const Real bdt);
void UserBoundary(Mesh* pm);

// Input parameters that need to be accessible in the user functions
static Real T_top, T_bot, rho_bot, Hbox;

//----------------------------------------------------------------------------------------
//  \brief Problem Generator for mass removal

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;

  int nscalars = pmbp->phydro->nscalars;
  int nhydro = pmbp->phydro->nhydro;

  // Enroll user functions 
  user_srcs_func = UserSource;
  user_bcs_func = UserBoundary;

  if (restart) return;

  // Get acceleration from input file
  if (!pin->DoesParameterExist("hydro", "const_accel_val")) {
    std::stringstream msg;
    msg << "hydro/const_accel_val is required in the input file.";
    throw std::runtime_error(msg.str());
  }
  Real g_accel = pin->GetReal("hydro", "const_accel_val");

  // Get top and bottom temperatures from input file
  if (!pin->DoesParameterExist("problem", "T_top") || 
      !pin->DoesParameterExist("problem", "T_bot")) {
    std::stringstream msg;
    msg << "problem/T_top and problem/T_bot are required in the input file.";
    throw std::runtime_error(msg.str());
  }
  T_top = pin->GetReal("problem", "T_top");
  T_bot = pin->GetReal("problem", "T_bot");

  // Get bottom density from input file
  if (!pin->DoesParameterExist("problem", "rho_bot")) {
    std::stringstream msg;
    msg << "problem/rho_bot is required in the input file.";
    throw std::runtime_error(msg.str());
  }
  rho_bot = pin->GetReal("problem", "rho_bot");
  Hbox = pin->GetReal("mesh", "x1max") - pin->GetReal("mesh", "x1min");

  // Read density perturbation amplitude and wavenumber from input file
  if (!pin->DoesParameterExist("problem", "density_perturbation_amp") ||
      !pin->DoesParameterExist("problem", "density_perturbation_k")) {
    std::stringstream msg;
    msg << "problem/density_perturbation_amp and problem/density_perturbation_k are required in the input file.";
    throw std::runtime_error(msg.str());
  }
  Real rho_perturb_amp = pin->GetReal("problem", "density_perturbation_amp");
  Real rho_perturb_k = pin->GetReal("problem", "density_perturbation_k");

  Real inv_Hbox = 1.0/Hbox;
  Real k_seed = 2.0*Kokkos::numbers::pi*inv_Hbox*rho_perturb_k;

  const Real T_bot_device = T_bot;
  const Real rho_bot_device = rho_bot;

  Real gradT = (T_top - T_bot)/Hbox;

  printf("T_bot_device in pgen: %f\n", T_bot_device);

  // Real den0 =  pow(T_bot/T_ref, 1. - g_accel/gradT);

  // printf("Convection problem parameters:\n");
  // printf("  Rayleigh number = %e\n", Ra);
  // printf("  Gravity = %e\n", g_accel);
  // printf("  Temperature gradient = %e\n", gradT);
  // printf("  delta T = %e\n", deltaT);
  // printf("  Reference Temperature = %e\n", T_ref);
  // printf("  Density at bottom = %e\n", rho_bot);
  // printf("  Density at top = %e\n", den0);
  // printf("  Tbot/Tref = %e\n", T_bot/T_ref);
  // printf("  g_accel/gradT = %e\n", g_accel/gradT);



  // Capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  // Initialize Hydro variables -------------------------------
  auto &u0 = pmbp->phydro->u0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;


  // Set initial conditions
  par_for("pgen_turb", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    Real R2 = x1v*x1v + x2v*x2v + x3v*x3v;

    Real temp = T_bot_device + gradT*x1v;
    Real den = rho_bot_device * pow(T_bot_device/temp, 1. - g_accel/gradT);

    // printf("x1=%f, temp=%f, den=%f\n", x1v, temp, den);

    u0(m,IDN,k,j,i) = den*(1.0 + 
      rho_perturb_amp * 
      std::cos(k_seed*x1v) * 
      std::cos(k_seed*x2v) * 
      std::cos(k_seed*x3v)
    );
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = den * temp/ gm1 +
       0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
       SQR(u0(m,IM3,k,j,i)))/den;
    
    Real dye_conc = 0.0;
    if (R2 < 0.1) { // RHS should be r^2
      dye_conc = 1.0;
    }

    u0(m, nhydro  , k, j, i) = dye_conc * den; // first scalar
    u0(m, nhydro+1, k, j, i) = dye_conc * den; // second scalar
    u0(m, nhydro+2, k, j, i) = dye_conc * den; // third scalar
    u0(m, nhydro+3, k, j, i) = dye_conc * den; // fourth scalar
    u0(m, nhydro+4, k, j, i) = 0.5 * den;      // fifth scalar - mean gradient forcing
  });
  


  return;
}

KOKKOS_INLINE_FUNCTION
// Advance x in [0,1] using a triangular heating profile.
Real heatStep(Real x, Real amplitude, Real dt)
{
  // Clamp input to [0,1] just in case
  x = Kokkos::clamp(x, 0.0, 1.0);

  // Piecewise-linear heating rate
  Real rate = amplitude * -1.0 * sin(x*2*Kokkos::numbers::pi);

  // Update and clamp to 1.0
  Real x_new = x + rate * dt;
  if (x_new > 1.0) x_new = 1.0;

  return x_new;
}

void UserSource(Mesh* pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  int nscalars = pmbp->phydro->nscalars;
  int nhydro = pmbp->phydro->nhydro;



  par_for("user_source", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real density = w0(m, IDN, k, j, i);
    Real scalar_conc_1 = w0(m, nhydro  , k, j, i);
    Real scalar_conc_2 = w0(m, nhydro+1, k, j, i);
    Real scalar_conc_3 = w0(m, nhydro+2, k, j, i);
    Real scalar_conc_4 = w0(m, nhydro+3, k, j, i);

    // Apply heating to first scalar
    Real new_scalar_conc_1 = heatStep(scalar_conc_1, 0.0, bdt);
    u0(m, nhydro  , k, j, i) += (new_scalar_conc_1 - scalar_conc_1) * density;

    // Apply heating to second scalar
    Real new_scalar_conc_2 = heatStep(scalar_conc_2, 0.1, bdt);
    u0(m, nhydro+1, k, j, i) += (new_scalar_conc_2 - scalar_conc_2) * density;

    // Apply heating to third scalar
    Real new_scalar_conc_3 = heatStep(scalar_conc_3, 1.0, bdt);
    u0(m, nhydro+2, k, j, i) += (new_scalar_conc_3 - scalar_conc_3) * density;

    // Apply heating to fourth scalar
    Real new_scalar_conc_4 = heatStep(scalar_conc_4, 10.0, bdt);
    u0(m, nhydro+3, k, j, i) += (new_scalar_conc_4 - scalar_conc_4) * density;

    // Appy Mean Gradient Forcing to fifth scalar
    // Passive scalar source term = u.G
    // let G = (ds/dx, 0, 0) , then u.G = v_x * ds/dx
    Real G = 0.5; // We set the x gradient to be small
    Real vx = w0(m, IVX, k, j, i);
    // Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    u0(m, nhydro+4, k, j, i) += vx * G * bdt;

    // We should check that scalars are guaranteed to be in [0,1] after all source terms are added.
    u0(m, nhydro  , k, j, i) = Kokkos::clamp(u0(m, nhydro  , k, j, i), 0.0, u0(m,IDN,k,j,i));
    u0(m, nhydro+1, k, j, i) = Kokkos::clamp(u0(m, nhydro+1, k, j, i), 0.0, u0(m,IDN,k,j,i));
    u0(m, nhydro+2, k, j, i) = Kokkos::clamp(u0(m, nhydro+2, k, j, i), 0.0, u0(m,IDN,k,j,i));
    u0(m, nhydro+3, k, j, i) = Kokkos::clamp(u0(m, nhydro+3, k, j, i), 0.0, u0(m,IDN,k,j,i));
    u0(m, nhydro+4, k, j, i) = Kokkos::clamp(u0(m, nhydro+4, k, j, i), 0.0, u0(m,IDN,k,j,i));


  });

  return;
}

//----------------------------------------------------------------------------------------
// User boundary condition mimicking reflective BCs. Activate by setting ix?/ox?_bc=user.
void UserBoundary(Mesh* pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pmbp->phydro == nullptr) return;

  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  auto &indcs = pm->mb_indcs;
  int ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;
  int &ie = indcs.ie;
  int &js = indcs.js;
  int &je = indcs.je;
  int &ks = indcs.ks;
  int &ke = indcs.ke;

  auto &u0 = pmbp->phydro->u0;
  int nvar = u0.extent_int(1);
  int nmb1 = pmbp->nmb_thispack - 1;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  const Real T_bot_device = T_bot;
  const Real T_top_device = T_top;

  printf("T_bot_device in boundary: %f\n", T_bot_device);
  printf("T_top_device in boundary: %f\n", T_top_device);

  // Handle X1 faces
  par_for("user_boundary_x1", DevExeSpace(), 0, nmb1, 0, (nvar-1), 0, (n3-1), 0, (n2-1),
  KOKKOS_LAMBDA(int m, int n, int k, int j) {
    if (mb_bcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
      for (int i = 0; i < ng; ++i) {
        int ghost = is - i - 1;
        int mirr = is + i;
        u0(m, n, k, j, ghost) = (n == IVX) ? -u0(m, n, k, j, mirr)
                                           :  u0(m, n, k, j, mirr);

        if (n == IEN) {
          Real den_mirr = u0(m, IDN, k, j, mirr);
          Real mx_mirr = u0(m, IM1, k, j, mirr);
          Real my_mirr = u0(m, IM2, k, j, mirr);
          Real mz_mirr = u0(m, IM3, k, j, mirr);
          Real kin = 0.5*(SQR(mx_mirr) + SQR(my_mirr) + SQR(mz_mirr))/den_mirr;
          u0(m, IEN, k, j, ghost) = den_mirr * T_bot_device/ gm1 + kin;
        }
      }
    }
    if (mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
      for (int i = 0; i < ng; ++i) {
        int ghost = ie + i + 1;
        int mirr = ie - i;
        u0(m, n, k, j, ghost) = (n == IVX) ? -u0(m, n, k, j, mirr)
                                           :  u0(m, n, k, j, mirr);

        if (n == IEN) {
          Real den_mirr = u0(m, IDN, k, j, mirr);
          Real mx_mirr = u0(m, IM1, k, j, mirr);
          Real my_mirr = u0(m, IM2, k, j, mirr);
          Real mz_mirr = u0(m, IM3, k, j, mirr);
          Real kin = 0.5*(SQR(mx_mirr) + SQR(my_mirr) + SQR(mz_mirr))/den_mirr;
          u0(m, IEN, k, j, ghost) = den_mirr * T_top_device/ gm1 + kin;
        }
      }
    }
  });

  // Handle X2 faces (if they exist)
  if (n2 > 1) {
    par_for("user_boundary_x2", DevExeSpace(), 0, nmb1, 0, (nvar-1), 0, (n3-1), 0, (n1-1),
    KOKKOS_LAMBDA(int m, int n, int k, int i) {
      if (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::user) {
        for (int j = 0; j < ng; ++j) {
          int ghost = js - j - 1;
          int mirr = js + j;
          u0(m, n, k, ghost, i) = (n == IVY) ? -u0(m, n, k, mirr, i)
                                             :  u0(m, n, k, mirr, i);
        }
      }
      if (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::user) {
        for (int j = 0; j < ng; ++j) {
          int ghost = je + j + 1;
          int mirr = je - j;
          u0(m, n, k, ghost, i) = (n == IVY) ? -u0(m, n, k, mirr, i)
                                             :  u0(m, n, k, mirr, i);
        }
      }
    });
  }

  // Handle X3 faces (if they exist)
  if (n3 > 1) {
    par_for("user_boundary_x3", DevExeSpace(), 0, nmb1, 0, (nvar-1), 0, (n2-1), 0, (n1-1),
    KOKKOS_LAMBDA(int m, int n, int j, int i) {
      if (mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user) {
        for (int k = 0; k < ng; ++k) {
          int ghost = ks - k - 1;
          int mirr = ks + k;
          u0(m, n, ghost, j, i) = (n == IVZ) ? -u0(m, n, mirr, j, i)
                                             :  u0(m, n, mirr, j, i);
        }
      }
      if (mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::user) {
        for (int k = 0; k < ng; ++k) {
          int ghost = ke + k + 1;
          int mirr = ke - k;
          u0(m, n, ghost, j, i) = (n == IVZ) ? -u0(m, n, mirr, j, i)
                                             :  u0(m, n, mirr, j, i);
        }
      }
    });
  }
}
