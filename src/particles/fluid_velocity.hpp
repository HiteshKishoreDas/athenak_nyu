// Copyright (C) 2026 the AthenaK contributors. BSD-3-Clause.
#ifndef PARTICLES_FLUID_VELOCITY_HPP_
#define PARTICLES_FLUID_VELOCITY_HPP_
#include "athena.hpp"
#include "mesh/mesh.hpp"

namespace particles {
// Cell-centered multilinear primitive velocity, shared by tracers and curve nodes.
// Return a diagnostic instead of reading outside the filled ghost stencil.
KOKKOS_INLINE_FUNCTION
int FluidVelocity(const DvceArray5D<Real> &w, int m, const Real x[3],
                  const RegionSize &s, const RegionIndcs &ind, bool three, Real v[3]) {
  for (int d=0; d<3; ++d) {
    v[d]=0;
    if (!Kokkos::isfinite(x[d])) return 2;
  }
  Real f[3] = {(x[0]-s.x1min)/s.dx1+ind.is-0.5,
                (x[1]-s.x2min)/s.dx2+ind.js-0.5,
                three ? (x[2]-s.x3min)/s.dx3+ind.ks-0.5 : 0};
  int cell[3] = {static_cast<int>(Kokkos::floor(f[0])),
                 static_cast<int>(Kokkos::floor(f[1])),
                 three ? static_cast<int>(Kokkos::floor(f[2])) : 0};
  if (cell[0]<0 || cell[0]+1>=static_cast<int>(w.extent(4)) ||
      cell[1]<0 || cell[1]+1>=static_cast<int>(w.extent(3)) ||
      cell[2]<0 || (three && cell[2]+1>=static_cast<int>(w.extent(2)))) return 1;
  for (int kz=0; kz<(three ? 2 : 1); ++kz) {
    for (int jy=0; jy<2; ++jy) {
      for (int ix=0; ix<2; ++ix) {
        int i=cell[0]+ix, j=cell[1]+jy, k=cell[2]+kz;
        Real weight=(ix ? f[0]-cell[0] : 1-f[0]+cell[0])*
                    (jy ? f[1]-cell[1] : 1-f[1]+cell[1]);
        if (three) weight *= kz ? f[2]-cell[2] : 1-f[2]+cell[2];
        if (!(w(m,IDN,k,j,i)>0) || !Kokkos::isfinite(w(m,IDN,k,j,i))) return 2;
        for (int d=0; d<(three ? 3 : 2); ++d) {
          Real u=w(m,IVX+d,k,j,i);
          if (!Kokkos::isfinite(u)) return 2;
          v[d] += weight*u;
        }
      }
    }
  }
  return 0;
}
}  // namespace particles
#endif  // PARTICLES_FLUID_VELOCITY_HPP_
