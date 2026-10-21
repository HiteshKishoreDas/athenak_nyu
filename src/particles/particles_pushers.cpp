//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particle_pushers.cpp
//  \brief

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "particles.hpp"
#include "particles/fluid_velocity.hpp"
#include "hydro/hydro.hpp"

namespace particles {
//----------------------------------------------------------------------------------------
//! \fn  void Particles::ParticlesPush
//  \brief

TaskStatus Particles::Push(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is;
  int js = indcs.js;
  int ks = indcs.ks;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &pi = prtcl_idata;
  auto &pr = prtcl_rdata;
  auto dt_ = (pmy_pack->pmesh->dt);
  auto gids = pmy_pack->gids;

  switch (pusher) {
    case ParticlesPusher::drift:

      par_for("part_update",DevExeSpace(),0,(nprtcl_thispack-1),
      KOKKOS_LAMBDA(const int p) {
        int m = pi(PGID,p) - gids;
        int ip = (pr(IPX,p) - mbsize.d_view(m).x1min)/mbsize.d_view(m).dx1 + is;
        pr(IPX,p) += 0.5*dt_*pr(IPVX,p);

        if (multi_d) {
          int jp = (pr(IPY,p) - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 + js;
          pr(IPY,p) += 0.5*dt_*pr(IPVY,p);
        }

        if (three_d) {
          int kp = (pr(IPZ,p) - mbsize.d_view(m).x3min)/mbsize.d_view(m).dx3 + ks;
          pr(IPZ,p) += 0.5*dt_*pr(IPVZ,p);
        }
      });

    break;
  case ParticlesPusher::lagrangian_tracer: {
    if (stage != 1 && stage != 2) Kokkos::abort("invalid Lagrangian RK2 stage");
    auto w = pmy_pack->phydro->w0;
    const auto ind = indcs;
    const bool three = three_d;
    const int nmb = pmy_pack->nmb_thispack;
    par_for("velocity_following_tracers",DevExeSpace(),0,nprtcl_thispack-1,
    KOKKOS_LAMBDA(const int p) {
      const int m = pi(PGID,p)-gids;
      if (m < 0 || m >= nmb) Kokkos::abort("invalid tracer leaf ownership");
      Real x[3]={pr(IPX,p),pr(IPY,p),pr(IPZ,p)}, v[3];
      if (FluidVelocity(w,m,x,mbsize.d_view(m),ind,three,v)) {
        Kokkos::abort("invalid Lagrangian velocity interpolation stencil");
      }
      for (int d=0; d<(three ? 3 : 2); ++d) {
        // Predictor velocity lives in the native real payload, so interstage MPI
        // migration carries it. The increment form handles periodic wraps without
        // averaging coordinates on opposite sides of the domain.
        if (stage == 1) {
          pr(2*d+1,p)=v[d];
          pr(2*d,p) += dt_*v[d];
        } else {
          pr(2*d,p) += 0.5*dt_*(v[d]-pr(2*d+1,p));
          pr(2*d+1,p)=v[d];
        }
        if (!Kokkos::isfinite(pr(2*d,p))) Kokkos::abort("nonfinite tracer position");
      }
    });
    break;
  }
  case ParticlesPusher::lagrangian_mc:
    return PushLagrangianMC(pdriver, stage);
  default:
    break;
  }

  return TaskStatus::complete;
}
} // namespace particles
