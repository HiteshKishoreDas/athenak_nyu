//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2026 Athena code contributors
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_sgs.cpp
//! \brief Uniform-grid Schmidt & Federrath (2011), equations (7)-(14).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "srcterms/turb_driver.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace hydro {
namespace {
void SGSFailure(const std::string &message) {
  std::cerr << "### FATAL ERROR: SGS: " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

int GlobalFailures(int local) {
#if MPI_PARALLEL_ENABLED
  int global;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return global;
#else
  return local;
#endif
}

KOKKOS_INLINE_FUNCTION
Real ClosureDifference(Real left, Real right) {
  Real difference = right-left;
  // Roundoff convention for closure gradients only. Production uses actual
  // differences, preserving compatible work even in a boosted frame.
  Real roundoff = 8.0*std::numeric_limits<Real>::epsilon()*
                  (fabs(left)+fabs(right));
  return (fabs(difference) <= roundoff) ? 0.0 : difference;
}

KOKKOS_INLINE_FUNCTION
bool Admissible(const Real u[5], Real K, const EOS_Data &eos) {
  for (int n=0; n<5; ++n) {
    if (!Kokkos::isfinite(u[n])) return false;
  }
  if (!Kokkos::isfinite(K) || K < 0.0 || u[IDN] <= 0.0 || u[IDN] < eos.dfloor) {
    return false;
  }
  Real e = u[IEN] - 0.5*(SQR(u[IM1])+SQR(u[IM2])+SQR(u[IM3]))/u[IDN];
  Real p = (eos.gamma-1.0)*e;
  return Kokkos::isfinite(e) && e > 0.0 && p >= eos.pfloor &&
         p/u[IDN] >= eos.tfloor && p/pow(u[IDN], eos.gamma) >= eos.sfloor;
}
}  // namespace

void Hydro::ConfigureSGS(ParameterInput *pin) {
  if (!pin->DoesBlockExist("sgs")) return;
  if (!pin->GetOrAddBoolean("sgs", "enabled", true)) return;
  if (pin->GetOrAddString("sgs", "model", "schmidt_federrath") !=
      "schmidt_federrath") SGSFailure("model must be schmidt_federrath");
  if (!peos->eos_data.is_ideal || pmy_pack->pcoord->is_special_relativistic ||
      pmy_pack->pcoord->is_general_relativistic) {
    SGSFailure("requires Newtonian ideal-gas HD");
  }
  for (const char *block : {"mhd", "radiation", "ion-neutral", "adm", "z4c",
                            "shearing_box", "particles", "frame_tracking",
                            "hydro_srcterms"}) {
    if (pin->DoesBlockExist(block)) SGSFailure(std::string("unsupported block: ")+block);
  }
  for (const char *option : {"ism_cooling", "cgm_cooling", "rel_cooling",
                             "const_accel", "fofc"}) {
    if (pin->DoesParameterExist("hydro", option) &&
        pin->GetBoolean("hydro", option)) {
      SGSFailure(std::string("unsupported hydro option: ")+option);
    }
  }
  if (pin->DoesBlockExist("cooling") &&
      pin->GetOrAddBoolean("cooling", "enabled", false)) {
    SGSFailure("unsupported cooling requires separate qualification");
  }
  for (const char *option : {"viscosity", "conductivity", "tdep_conductivity"}) {
    if (pin->DoesParameterExist("hydro", option)) {
      SGSFailure(std::string("unsupported hydro option: ")+option);
    }
  }
  // ponytail: no SGS mesh transfer; qualify AMR energy exchange before adding hooks.
  if (pmy_pack->pmesh->multilevel) SGSFailure("AMR/SMR requires separate qualification");
#if SINGLE_PRECISION_ENABLED
  SGSFailure("single precision is not qualified");
#endif
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL)
  SGSFailure("GPU execution is not qualified");
#endif
  if (pin->GetString("time", "evolution") != "dynamic" ||
      pin->GetOrAddString("time", "integrator", "rk2") != "rk2" ||
      pin->GetString("hydro", "rsolver") != "llf") {
    SGSFailure("requires dynamic evolution, rk2 integrator and llf solver");
  }
  std::string recon = pin->GetOrAddString("hydro", "reconstruct", "plm");
  if (recon != "plm" && recon != "dc" && recon != "ppm4") {
    SGSFailure("requires plm, dc or ppm4 reconstruction");
  }
  auto *pm = pmy_pack->pmesh;
  int ndim = 1 + static_cast<int>(pm->multi_d) + static_cast<int>(pm->three_d);
  if (pm->mb_indcs.ng < 2) SGSFailure("requires at least two ghost zones");
  for (int f=0; f<2*ndim; ++f) {
    if (pm->mesh_bcs[f] != BoundaryFlag::periodic) {
      SGSFailure("only periodic boundaries are qualified");
    }
  }
  Real dx = pm->mesh_size.dx1;
  if (dx <= 0.0 || (pm->multi_d && fabs(pm->mesh_size.dx2/dx-1.0)>1.e-12) ||
      (pm->three_d && fabs(pm->mesh_size.dx3/dx-1.0)>1.e-12)) {
    SGSFailure("requires cubic cells in active dimensions");
  }
  if (!(pm->cfl_no > 0.0 && pm->cfl_no <= 0.2)) {
    SGSFailure("requires 0 < time/cfl_number <= 0.2");
  }
  if (!(peos->eos_data.gamma > 1.0) || !std::isfinite(peos->eos_data.gamma) ||
      nscalars < 0) SGSFailure("invalid gamma or scalar count");
  sgs_coeff.c1 = pin->GetOrAddReal("sgs", "c1", 0.02);
  sgs_coeff.c2 = pin->GetOrAddReal("sgs", "c2", 0.7);
  sgs_coeff.ceps = pin->GetOrAddReal("sgs", "ceps", 1.5);
  sgs_coeff.ckappa = pin->GetOrAddReal("sgs", "ckappa", 0.65);
  std::string init = pin->GetOrAddString("sgs", "initialization", "uniform");
  if (init != "uniform" && init != "pgen") SGSFailure("unknown initialization policy");
  sgs_pgen_initialization = (init == "pgen");
  sgs_seed = pin->GetOrAddReal("sgs", "seed", 1.e-6);
  for (Real x : {sgs_coeff.c1, sgs_coeff.c2, sgs_coeff.ceps,
                  sgs_coeff.ckappa, sgs_seed}) {
    if (!std::isfinite(x) || x < 0.0) SGSFailure("coefficients/seed must be finite >= 0");
  }
  if (sgs_coeff.c2 > 1.0) SGSFailure("c2 must be <= 1");
  isgs = nvars++;
  sgs_produced = pin->GetOrAddReal("sgs", "integrated_production", 0.0);
  sgs_dissipated = pin->GetOrAddReal("sgs", "integrated_dissipation", 0.0);
  sgs_limited = pin->GetOrAddReal("sgs", "integrated_limited", 0.0);
  sgs_diff_limited = pin->GetOrAddReal("sgs", "integrated_diffusion_limited", 0.0);
  sgs_fallbacks = pin->GetOrAddInteger("sgs", "fallbacks", 0);
  if (global_variable::my_rank == 0) {
    std::cout << "SGS schmidt_federrath: K slot=" << isgs << ", seed k=" << sgs_seed
              << ", c1=" << sgs_coeff.c1 << ", c2=" << sgs_coeff.c2
              << ", ceps=" << sgs_coeff.ceps << ", ckappa=" << sgs_coeff.ckappa
              << ", Delta=dx=" << dx << std::endl;
  }
}

void Hydro::InitializeSGS() {
  auto ind = pmy_pack->pmesh->mb_indcs;
  auto u = u0;
  int idx = isgs;
  Real seed = sgs_seed;
  par_for("sgs_seed", DevExeSpace(), 0, pmy_pack->nmb_thispack-1,
          ind.ks, ind.ke, ind.js, ind.je, ind.is, ind.ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u(m,idx,k,j,i) = seed*u(m,IDN,k,j,i);
  });
}

// One shared traction per face. Normal derivatives use adjacent cell differences;
// transverse derivatives average centered cell derivatives. No tensor-volume storage.
void Hydro::SGSFluxes() {
  auto w = w0;
  auto ind = pmy_pack->pmesh->mb_indcs;
  auto size = pmy_pack->pmb->mb_size;
  auto coeff = sgs_coeff;
  int idx = isgs;
  int ndim = 1 + pmy_pack->pmesh->multi_d + pmy_pack->pmesh->three_d;
  Real dt = pmy_pack->pmesh->dt;
  for (int dir=0; dir<ndim; ++dir) {
    auto f = (dir == 0) ? uflx.x1f : ((dir == 1) ? uflx.x2f : uflx.x3f);
    auto t = (dir == 0) ? sgs_face.x1f : ((dir == 1) ? sgs_face.x2f : sgs_face.x3f);
    int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("sgs_faces", DevExeSpace(), 0, pmy_pack->nmb_thispack-1,
            ind.ks, ind.ke+dk, ind.js, ind.je+dj, ind.is, ind.ie+di,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      int il = i-di, jl = j-dj, kl = k-dk;
      Real dx = size.d_view(m).dx1;
      Real rho_l = w(m,IDN,kl,jl,il), rho_r = w(m,IDN,k,j,i);
      Real K_l = rho_l*w(m,idx,kl,jl,il), K_r = rho_r*w(m,idx,k,j,i);
      Real rho = 0.5*(rho_l+rho_r), K = 0.5*(K_l+K_r);
      Real a[3][3] = {};
      for (int v=0; v<3; ++v) {
        for (int d=0; d<ndim; ++d) {
          int qi = (d == 0), qj = (d == 1), qk = (d == 2);
          if (d == dir) {
            a[v][d] = ClosureDifference(w(m,IVX+v,kl,jl,il),w(m,IVX+v,k,j,i))/dx;
          } else {
            a[v][d] = (ClosureDifference(w(m,IVX+v,k-qk,j-qj,i-qi),
                                          w(m,IVX+v,k+qk,j+qj,i+qi)) +
                       ClosureDifference(w(m,IVX+v,kl-qk,jl-qj,il-qi),
                                          w(m,IVX+v,kl+qk,jl+qj,il+qi)))/(4.0*dx);
          }
        }
      }
      Real tau[3][3];
      sgs::Stress(rho, K, dx, a, coeff, tau);
      Real production = 0.0, work = 0.0;
      for (int v=0; v<3; ++v) {
        Real vl = w(m,IVX+v,kl,jl,il), vr = w(m,IVX+v,k,j,i);
        production += tau[v][dir]*(vr-vl)/(2.0*dx);
        work += tau[v][dir]*0.5*(vl+vr);
      }
      Real theta = 1.0;
      if (production < 0.0 && dt > 0.0) {
        theta = fmin(1.0, fmin(K_l,K_r)/(8.0*ndim*dt*(-production)));
      }
      for (int v=0; v<3; ++v) {
        t(m,v,k,j,i) = theta*tau[v][dir];
        f(m,IM1+v,k,j,i) -= t(m,v,k,j,i);
      }
      t(m,3,k,j,i) = (1.0-theta)*fabs(production);
      f(m,IEN,k,j,i) -= theta*work;
      // SGS diffusion is a represented-energy flux, not a paired E source.
      Real J = -sgs::Diffusivity(rho, K, dx, coeff)*
                (w(m,idx,k,j,i)-w(m,idx,kl,jl,il))/dx;
      Real raw_J = J;
      if (J != 0.0 && dt > 0.0) {
        Real donor = (J > 0.0) ? K_l : K_r;
        J *= fmin(1.0, donor*dx/(8.0*ndim*dt*fabs(J)));
      }
      t(m,4,k,j,i) = fabs(raw_J-J)/dx;
      f(m,idx,k,j,i) += J;
    });
  }
}

void Hydro::SGSRates() {
  auto w = w0;
  auto ind = pmy_pack->pmesh->mb_indcs;
  auto size = pmy_pack->pmb->mb_size;
  auto coeff = sgs_coeff;
  auto r = sgs_rates;
  auto t1 = sgs_face.x1f, t2 = sgs_face.x2f, t3 = sgs_face.x3f;
  int idx = isgs;
  int ndim = 1 + pmy_pack->pmesh->multi_d + pmy_pack->pmesh->three_d;
  par_for("sgs_rates", DevExeSpace(), 0, pmy_pack->nmb_thispack-1,
          ind.ks, ind.ke, ind.js, ind.je, ind.is, ind.ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real sigma = 0.0, limited = 0.0, diff_limited = 0.0;
    Real dx = size.d_view(m).dx1;
    for (int dir=0; dir<ndim; ++dir) {
      auto t = (dir == 0) ? t1 : ((dir == 1) ? t2 : t3);
      int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
      for (int v=0; v<3; ++v) {
        Real uc = w(m,IVX+v,k,j,i);
        sigma += (t(m,v,k+dk,j+dj,i+di)*(w(m,IVX+v,k+dk,j+dj,i+di)-uc) +
                  t(m,v,k,j,i)*(uc-w(m,IVX+v,k-dk,j-dj,i-di)))/(2.0*dx);
      }
      limited += t(m,3,k,j,i) + t(m,3,k+dk,j+dj,i+di);
      diff_limited += 0.5*(t(m,4,k,j,i) + t(m,4,k+dk,j+dj,i+di));
    }
    r(m,0,k,j,i) = sigma;
    r(m,1,k,j,i) = sgs::Dissipation(w(m,IDN,k,j,i),
                    w(m,IDN,k,j,i)*w(m,idx,k,j,i), dx, coeff);
    r(m,2,k,j,i) = limited;
    r(m,3,k,j,i) = diff_limited;
  });
}

int Hydro::SGSTrial() {
  auto ind = pmy_pack->pmesh->mb_indcs;
  auto u = u0;
  auto size = pmy_pack->pmb->mb_size;
  auto f1 = uflx.x1f, f2 = uflx.x2f, f3 = uflx.x3f;
  auto rates = sgs_rates;
  auto eos = peos->eos_data;
  int idx = isgs;
  bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  Real dt = pmy_pack->pmesh->dt;
  // Check the full-dt FE substep of SSP RK2 (also in the second stage).
  auto *turb = pmy_pack->pturb;
  bool forcing = turb != nullptr && turb->stage_consistent && turb->IsDrivingActive();
  auto force = forcing ? turb->force : DvceArray5D<Real>();
  auto w = w0;
  int bad = 0;
  Kokkos::parallel_reduce("sgs_trial", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      {0,ind.ks,ind.js,ind.is}, {pmy_pack->nmb_thispack,ind.ke+1,ind.je+1,ind.ie+1}),
  KOKKOS_LAMBDA(int m, int k, int j, int i, int &sum) {
    Real trial[5], K;
    for (int n=0; n<=5; ++n) {
      int var = (n == 5) ? idx : n;
      Real div = (f1(m,var,k,j,i+1)-f1(m,var,k,j,i))/size.d_view(m).dx1;
      if (md) div += (f2(m,var,k,j+1,i)-f2(m,var,k,j,i))/size.d_view(m).dx2;
      if (td) div += (f3(m,var,k+1,j,i)-f3(m,var,k,j,i))/size.d_view(m).dx3;
      Real val = u(m,var,k,j,i)-dt*div;
      Real exchange = rates(m,0,k,j,i)-rates(m,1,k,j,i);
      if (n == IEN) val -= dt*exchange;
      if (n == 5) {
        K = val+dt*exchange;
      } else {
        trial[n] = val;
      }
    }
    if (forcing) {
      Real rho = w(m,IDN,k,j,i);
      for (int v=0; v<3; ++v) {
        Real acceleration = force(m,v,k,j,i);
        trial[IM1+v] += dt*rho*acceleration;
        trial[IEN] += dt*rho*w(m,IVX+v,k,j,i)*acceleration;
      }
    }
    if (!Admissible(trial, K, eos)) ++sum;
  }, bad);
  return GlobalFailures(bad);
}

void Hydro::SGSSources(Driver *d, int stage) {
  auto ind = pmy_pack->pmesh->mb_indcs;
  auto u = u0;
  auto r = sgs_rates;
  auto size = pmy_pack->pmb->mb_size;
  int idx = isgs;
  Real dt = d->beta[stage-1]*pmy_pack->pmesh->dt;
  par_for("sgs_exchange", DevExeSpace(), 0, pmy_pack->nmb_thispack-1,
          ind.ks, ind.ke, ind.js, ind.je, ind.is, ind.ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real exchange = dt*(r(m,0,k,j,i)-r(m,1,k,j,i));
    u(m,idx,k,j,i) += exchange;
    u(m,IEN,k,j,i) -= exchange;
  });
  Real prod=0.0, diss=0.0, limited=0.0, diff_limited=0.0;
  Kokkos::parallel_reduce("sgs_budget", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      {0,ind.ks,ind.js,ind.is}, {pmy_pack->nmb_thispack,ind.ke+1,ind.je+1,ind.ie+1}),
  KOKKOS_LAMBDA(int m, int k, int j, int i, Real &p, Real &e, Real &l, Real &jl) {
    Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    p += vol*r(m,0,k,j,i);
    e += vol*r(m,1,k,j,i);
    l += vol*r(m,2,k,j,i);
    jl += vol*r(m,3,k,j,i);
  }, prod, diss, limited, diff_limited);
#if MPI_PARALLEL_ENABLED
  Real budget[4] = {prod, diss, limited, diff_limited}, total[4];
  MPI_Allreduce(budget, total, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  prod = total[0];
  diss = total[1];
  limited = total[2];
  diff_limited = total[3];
#endif
  // Final-solution quadrature weights differ from stage update coefficients.
  Real weight = d->SourceTermHistoryWeight(stage)*pmy_pack->pmesh->dt;
  sgs_produced += weight*prod;
  sgs_dissipated += weight*diss;
  sgs_limited += weight*limited;
  sgs_diff_limited += weight*diff_limited;
}

void Hydro::CheckSGSState() {
  auto ind = pmy_pack->pmesh->mb_indcs;
  auto u = u0;
  auto eos = peos->eos_data;
  int idx = isgs;
  int bad = 0;
  // Check active cells before the EOS gets any chance to apply a floor.
  Kokkos::parallel_reduce("sgs_admissible", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      {0,ind.ks,ind.js,ind.is}, {pmy_pack->nmb_thispack,ind.ke+1,ind.je+1,ind.ie+1}),
  KOKKOS_LAMBDA(int m, int k, int j, int i, int &sum) {
    Real state[5];
    for (int n=0; n<5; ++n) state[n] = u(m,n,k,j,i);
    if (!Admissible(state,u(m,idx,k,j,i),eos)) ++sum;
  }, bad);
  if (GlobalFailures(bad)>0) {
    SGSFailure("inadmissible state; SGS runs never apply floors");
  }
}

void Hydro::SGSNewTimeStep() {
  auto ind = pmy_pack->pmesh->mb_indcs;
  auto w = w0;
  auto size = pmy_pack->pmb->mb_size;
  auto eos = peos->eos_data;
  auto c = sgs_coeff;
  int idx = isgs;
  int ndim = 1 + pmy_pack->pmesh->multi_d + pmy_pack->pmesh->three_d;
  Real sgs_dt = std::numeric_limits<Real>::max();
  Kokkos::parallel_reduce("sgs_dt", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      {0,ind.ks,ind.js,ind.is}, {pmy_pack->nmb_thispack,ind.ke+1,ind.je+1,ind.ie+1}),
  KOKKOS_LAMBDA(int m, int k, int j, int i, Real &dtmin) {
    Real dx = size.d_view(m).dx1, specific = w(m,idx,k,j,i);
    Real rho = w(m,IDN,k,j,i), K = rho*specific;
    Real wave_k = specific, rate = 0.0;
    // Frozen-coefficient row sum uses receiving density, including contrasts.
    for (int dir=0; dir<ndim; ++dir) {
      int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
      for (int side=-1; side<=1; side+=2) {
        Real rho_n = w(m,IDN,k+side*dk,j+side*dj,i+side*di);
        Real K_n = rho_n*w(m,idx,k+side*dk,j+side*dj,i+side*di);
        Real rho_f = 0.5*(rho+rho_n), K_f = 0.5*(K+K_n);
        wave_k = fmax(wave_k,K_f/rho);
        Real mu = 2.0*c.c1*dx*sqrt(2.0)*sqrt(rho_f)*sqrt(K_f);
        Real kappa = sgs::Diffusivity(rho_f,K_f,dx,c);
        rate += (2.0*mu+kappa)/(rho*dx*dx);
      }
    }
    // An identically zero SGS field imposes no extra timestep constraint.
    if (wave_k == 0.0) return;
    Real cs = sqrt(eos.gamma*(eos.gamma-1.0)*w(m,IEN,k,j,i)/rho+6.0*wave_k);
    Real speed = ndim*cs;
    for (int v=0; v<ndim; ++v) speed += fabs(w(m,IVX+v,k,j,i));
    Real local = dx/speed;
    if (rate > 0.0) local = fmin(local,1.0/rate);
    if (c.ceps > 0.0 && specific > 0.0) {
      local = fmin(local,dx/(2.0*c.ceps*sqrt(specific)));
    }
    dtmin = fmin(dtmin,local);
  }, Kokkos::Min<Real>(sgs_dt));
  dtnew = std::min(dtnew,sgs_dt);
}
}  // namespace hydro
