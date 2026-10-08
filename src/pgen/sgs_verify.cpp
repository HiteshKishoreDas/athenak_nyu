//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sgs_verify.cpp
//! \brief Deterministic periodic initial conditions for HD SGS verification.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "driver/driver.hpp"
#include "globals.hpp"
#include "pgen/pgen.hpp"
#include "srcterms/turb_driver.hpp"

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (pin->GetOrAddBoolean("problem", "user_source", false)) {
    user_srcs = true;
    user_srcs_func = [](Mesh *, const Real) {};
  }
  if (restart) return;
  auto *pack = pmy_mesh_->pmb_pack;
  if (pack->phydro == nullptr || !pack->phydro->peos->eos_data.is_ideal) {
    std::cerr << "sgs_verify requires ideal-gas hydro" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // Test-only input: nx, ny, nvar; then y/x ordered conserved fields and K.
  // E excludes K. SGS-off files must provide zero K rather than lose energy.
  if (pin->DoesParameterExist("problem", "initial_state")) {
    const int nx = pin->GetInteger("mesh", "nx1");
    const int ny = pin->GetInteger("mesh", "nx2");
    const int nvar = pack->phydro->nhydro + pack->phydro->nscalars + 1;
    auto fail = []() {
      std::cerr << "sgs_verify: invalid 2D initial_state file/layout/reservoir"
                << std::endl;
      std::exit(EXIT_FAILURE);
    };
    if (pin->GetInteger("mesh", "nx3") != 1 || pmy_mesh_->multilevel ||
        (pack->phydro->isgs >= 0 && !pack->phydro->sgs_pgen_initialization)) fail();
    std::ifstream input(pin->GetString("problem", "initial_state"));
    int file_nx = 0, file_ny = 0, file_nvar = 0;
    if (!(input >> file_nx >> file_ny >> file_nvar) || file_nx != nx ||
        file_ny != ny || file_nvar != nvar) fail();
    HostArray3D<Real> host("sgs_initial_host", ny, nx, nvar);
    for (int j=0; j<ny; ++j) {
      for (int i=0; i<nx; ++i) {
        for (int n=0; n<nvar; ++n) {
          if (!(input >> host(j,i,n)) || !std::isfinite(host(j,i,n))) fail();
        }
        Real rho = host(j,i,IDN);
        Real kinetic = 0.0;
        for (int n=IM1; n<=IM3; ++n) kinetic += SQR(host(j,i,n));
        if (rho <= 0.0 || host(j,i,IEN) <= kinetic/(2.0*rho) ||
            host(j,i,nvar-1) < 0.0 ||
            (pack->phydro->isgs < 0 && host(j,i,nvar-1) != 0.0)) fail();
      }
    }
    std::string extra;
    if (input >> extra) fail();
    DvceArray3D<Real> data("sgs_initial_data", ny, nx, nvar);
    Kokkos::deep_copy(data, host);
    const auto ind = pmy_mesh_->mb_indcs;
    const auto domain = pmy_mesh_->mesh_size;
    auto size = pack->pmb->mb_size;
    auto u = pack->phydro->u0;
    const int count = pack->phydro->nvars;
    par_for("pgen_sgs_initial_state", DevExeSpace(), 0, pack->nmb_thispack-1,
            ind.ks, ind.ke, ind.js, ind.je, ind.is, ind.ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real x = CellCenterX(i-ind.is, ind.nx1,
                           size.d_view(m).x1min, size.d_view(m).x1max);
      Real y = CellCenterX(j-ind.js, ind.nx2,
                           size.d_view(m).x2min, size.d_view(m).x2max);
      int ix = static_cast<int>(nx*(x-domain.x1min)/(domain.x1max-domain.x1min));
      int iy = static_cast<int>(ny*(y-domain.x2min)/(domain.x2max-domain.x2min));
      for (int n=0; n<count; ++n) u(m,n,k,j,i) = data(iy,ix,n);
    });
    return;
  }
  const std::string mode = pin->GetOrAddString("problem", "mode", "uniform");
  int which = -1;
  const char *modes[] = {"uniform", "shear", "diffusion", "wave", "contact",
                         "collision", "box"};
  for (int n=0; n<7; ++n) {
    if (mode == modes[n]) which = n;
  }
  if (which < 0) {
    std::cerr << "sgs_verify: unknown problem/mode " << mode << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const Real rho0 = pin->GetOrAddReal("problem", "density0", 1.0);
  const Real p0 = pin->GetOrAddReal("problem", "pressure0", 1.0);
  const Real k0 = pin->GetOrAddReal("problem", "k0", 0.01);
  const Real amp = pin->GetOrAddReal("problem", "amplitude", 0.2);
  const Real v0 = pin->GetOrAddReal("problem", "velocity", 0.3);
  if (!std::isfinite(rho0) || rho0 <= 0.0 || !std::isfinite(p0) || p0 <= 0.0 ||
      !std::isfinite(k0) || k0 < 0.0 || !std::isfinite(amp) || amp < 0.0 ||
      amp >= 1.0 || !std::isfinite(v0)) {
    std::cerr << "sgs_verify: require density0,pressure0 > 0; k0 >= 0; "
              << "0 <= amplitude < 1; finite velocity" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (which == 2 && p0 <= (2.0/3.0)*rho0*k0*amp) {
    std::cerr << "sgs_verify diffusion: pressure0 too small for pressure balance"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const auto ind = pmy_mesh_->mb_indcs;
  const auto mesh_size = pmy_mesh_->mesh_size;
  const auto size = pack->pmb->mb_size;
  auto u = pack->phydro->u0;
  const int nhyd = pack->phydro->nhydro;
  const int nscal = pack->phydro->nscalars;
  const int isgs = pack->phydro->isgs;
  const bool profile = isgs >= 0 && pack->phydro->sgs_pgen_initialization;
  const Real gamma = pack->phydro->peos->eos_data.gamma;
  const Real gm1 = gamma-1.0;
  const Real twopi = 2.0*acos(-1.0);
  const bool multi_d = pmy_mesh_->multi_d;
  const bool three_d = pmy_mesh_->three_d;
  const bool mixed_wave = pin->GetOrAddBoolean("problem", "wave_linear_mixed", false);
  const Real longitudinal_b = 2.0/3.0 + 4.0*pack->phydro->sgs_coeff.c2/3.0;
  if (mixed_wave && (which != 3 || multi_d || isgs < 0 ||
      pack->phydro->sgs_coeff.c1 != 0.0 || pack->phydro->sgs_coeff.ceps != 0.0 ||
      pack->phydro->sgs_coeff.ckappa != 0.0)) {
    std::cerr << "mixed wave oracle requires 1D SGS wave, c1=ceps=ckappa=0"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (which == 1 && !multi_d) {
    std::cerr << "sgs_verify shear requires at least two dimensions" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  par_for("pgen_sgs_verify", DevExeSpace(), 0, pack->nmb_thispack-1,
          ind.ks, ind.ke, ind.js, ind.je, ind.is, ind.ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(i-ind.is, ind.nx1,
                               size.d_view(m).x1min, size.d_view(m).x1max);
    const Real y = CellCenterX(j-ind.js, ind.nx2,
                               size.d_view(m).x2min, size.d_view(m).x2max);
    const Real z = CellCenterX(k-ind.ks, ind.nx3,
                               size.d_view(m).x3min, size.d_view(m).x3max);
    const Real ax = twopi*(x-mesh_size.x1min)/(mesh_size.x1max-mesh_size.x1min);
    const Real ay = twopi*(y-mesh_size.x2min)/(mesh_size.x2max-mesh_size.x2min);
    const Real az = twopi*(z-mesh_size.x3min)/(mesh_size.x3max-mesh_size.x3min);
    Real rho = rho0, p = p0, K = rho0*k0;
    Real vx = v0, vy = 0.0, vz = 0.0;
    if (which == 1) {
      vx = v0*sin(ay);
    } else if (which == 2) {
      vx = 0.0;
      K *= 1.0+amp*cos(ax);
      // Initial effective pressure is constant: P + 2K/3 = p0 + 2rho0*k0/3.
      p -= (2.0/3.0)*(K-rho0*k0);
    } else if (which == 3) {
      const Real perturbation = amp*sin(ax+(mixed_wave ? 0.123 : 0.0));
      rho *= 1.0+perturbation;
      p *= pow(rho/rho0, gamma);
      K *= pow(rho/rho0, mixed_wave ? 1.0+longitudinal_b : 5.0/3.0);
      const Real acoustic = mixed_wave ? longitudinal_b*(1.0+longitudinal_b) :
                                        10.0/9.0;
      vx = sqrt(gamma*p0/rho0+acoustic*k0)*perturbation;
    } else if (which == 4) {
      rho *= 1.0+amp*cos(ax);
      // Constant thermal pressure and K density give a moving material contact.
    } else if (which == 5) {
      vx = -v0*sin(ax);
    } else if (which == 6) {
      rho *= 1.0+amp*cos(ax)*(multi_d ? cos(ay) : 1.0);
      vx = v0*(multi_d ? sin(ay) : sin(ax));
      vy = multi_d ? v0*(three_d ? sin(az) : sin(ax)) : 0.0;
      vz = three_d ? v0*sin(ax) : 0.0;
      K = rho*k0*(1.0+0.5*amp*cos(ax));
    }
    u(m,IDN,k,j,i) = rho;
    u(m,IM1,k,j,i) = rho*vx;
    u(m,IM2,k,j,i) = rho*vy;
    u(m,IM3,k,j,i) = rho*vz;
    u(m,IEN,k,j,i) = p/gm1+0.5*rho*(vx*vx+vy*vy+vz*vz);
    for (int n=0; n<nscal; ++n) {
      // Different positive tracer profiles expose an aliased SGS/passive slot.
      const Real scalar = 0.2+0.05*n+(mixed_wave ? 0.0 : 0.1*sin((n+1)*ax));
      u(m,nhyd+n,k,j,i) = rho*scalar;
    }
    if (profile) u(m,isgs,k,j,i) = K;
  });
  if (pin->GetOrAddBoolean("problem", "forcing_oracle", false)) {
    auto *turb = pack->pturb;
    if (which != 0 || turb == nullptr || !turb->stage_consistent ||
        !turb->IsDrivingActive()) {
      std::cerr << "forcing_oracle requires uniform mode and active stage forcing"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    DvceArray5D<Real> saved_u("oracle_u", u.extent(0), u.extent(1),
                              u.extent(2), u.extent(3), u.extent(4));
    auto w = pack->phydro->w0;
    DvceArray5D<Real> saved_w("oracle_w", w.extent(0), w.extent(1),
                              w.extent(2), w.extent(3), w.extent(4));
    auto force = turb->force;
    DvceArray5D<Real> saved_force("oracle_force", force.extent(0), force.extent(1),
                                  force.extent(2), force.extent(3), force.extent(4));
    Kokkos::deep_copy(saved_u, u);
    Kokkos::deep_copy(saved_w, w);
    Kokkos::deep_copy(saved_force, force);
    Real saved_dt = pmy_mesh_->dt;
    Real saved_work = turb->integrated_work;
    Real saved_mom[3];
    for (int d=0; d<3; ++d) saved_mom[d] = turb->integrated_momentum[d];
    const Real test_dt = 0.01;
    const Real accel[3] = {0.7, -0.2, 0.3};
    pmy_mesh_->dt = test_dt;
    pack->phydro->peos->ConsToPrim(u, w, false,
                        ind.is, ind.ie, ind.js, ind.je, ind.ks, ind.ke);
    par_for("oracle_force", DevExeSpace(), 0, pack->nmb_thispack-1,
            ind.ks, ind.ke, ind.js, ind.je, ind.is, ind.ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      for (int d=0; d<3; ++d) force(m,d,k,j,i) = accel[d];
    });
    Driver oracle_driver(pin, pmy_mesh_, 0.0, nullptr);
    turb->StageForcing(&oracle_driver, 1);
    pack->phydro->peos->ConsToPrim(u, w, false,
                        ind.is, ind.ie, ind.js, ind.je, ind.ks, ind.ke);
    const int nvars = pack->phydro->nvars;
    par_for("oracle_rk2_blend", DevExeSpace(), 0, pack->nmb_thispack-1,
            ind.ks, ind.ke, ind.js, ind.je, ind.is, ind.ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      for (int n=0; n<nvars; ++n) {
        u(m,n,k,j,i) = 0.5*saved_u(m,n,k,j,i)+0.5*u(m,n,k,j,i);
      }
    });
    turb->StageForcing(&oracle_driver, 2);
    auto before = Kokkos::create_mirror_view_and_copy(HostMemSpace(), saved_u);
    auto after = Kokkos::create_mirror_view_and_copy(HostMemSpace(), u);
    const Real expected_work = rho0*(v0*accel[0]*test_dt+
                   0.5*(accel[0]*accel[0]+accel[1]*accel[1]+accel[2]*accel[2])*
                   test_dt*test_dt);
    Real max_error = 0.0;
    for (int m=0; m<pack->nmb_thispack; ++m) {
      for (int k=ind.ks; k<=ind.ke; ++k) {
        for (int j=ind.js; j<=ind.je; ++j) {
          for (int i=ind.is; i<=ind.ie; ++i) {
            for (int n=0; n<nvars; ++n) {
              Real expected = before(m,n,k,j,i);
              if (n >= IM1 && n <= IM3) expected += rho0*accel[n-IM1]*test_dt;
              if (n == IEN) expected += expected_work;
              Real value = after(m,n,k,j,i);
              if (!std::isfinite(value)) {
                max_error = std::numeric_limits<Real>::infinity();
              }
              max_error = std::max(max_error, fabs(value-expected));
            }
          }
        }
      }
    }
    Real volume = (mesh_size.x1max-mesh_size.x1min)*(mesh_size.x2max-mesh_size.x2min)*
                  (mesh_size.x3max-mesh_size.x3min);
    max_error = std::max(max_error,
                    fabs(turb->integrated_work-saved_work-expected_work*volume));
    if (!std::isfinite(turb->integrated_work)) {
      max_error = std::numeric_limits<Real>::infinity();
    }
    for (int d=0; d<3; ++d) {
      max_error = std::max(max_error,
          fabs(turb->integrated_momentum[d]-saved_mom[d]-rho0*accel[d]*test_dt*volume));
      if (!std::isfinite(turb->integrated_momentum[d])) {
        max_error = std::numeric_limits<Real>::infinity();
      }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &max_error, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
    Kokkos::deep_copy(u, saved_u);
    Kokkos::deep_copy(w, saved_w);
    Kokkos::deep_copy(force, saved_force);
    pmy_mesh_->dt = saved_dt;
    turb->integrated_work = saved_work;
    for (int d=0; d<3; ++d) turb->integrated_momentum[d] = saved_mom[d];
    if (!std::isfinite(max_error) || max_error > 1.e-12) {
      std::cerr << "forcing-oracle failed: max error=" << max_error << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (global_variable::my_rank == 0) {
      std::cout << "forcing-oracle passed: max error=" << max_error
                << ", tolerance=1e-12" << std::endl;
    }
  }
}
