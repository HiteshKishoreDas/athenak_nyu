// Copyright (C) 2026 the AthenaK contributors
// Licensed under the 3-clause BSD License (the "LICENSE").
#ifndef SHOCKS_OBSERVATION_HPP_
#define SHOCKS_OBSERVATION_HPP_
#include "shocks/detector.hpp"
#include "shocks/encounter.hpp"

namespace shocks {
KOKKOS_INLINE_FUNCTION
Observation ObserveSignal(const Signal &s, const Real x[3], Real time,
                          uint64_t generation) {
  Observation o;
  o.time=time; o.generation=generation; o.valid=s.valid; o.inside=s.inside;
  o.rho=s.rho; o.pressure=s.pressure; o.temperature=s.temperature;
  o.mach_pressure=s.mach_pressure; o.mach_temperature=s.mach_temperature;
  for (int d=0; d<3; ++d) {o.x[d]=x[d]; o.normal[d]=s.normal[d];}
  o.upstream[0]=s.rho_up; o.upstream[1]=s.pressure_up;
  o.upstream[2]=s.temperature_up;
  o.downstream[0]=s.rho_down; o.downstream[1]=s.pressure_down;
  o.downstream[2]=s.temperature_down;
  return o;
}
}  // namespace shocks
#endif  // SHOCKS_OBSERVATION_HPP_
