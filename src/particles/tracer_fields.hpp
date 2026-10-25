#ifndef PARTICLES_TRACER_FIELDS_HPP_
#define PARTICLES_TRACER_FIELDS_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file tracer_fields.hpp
//! \brief parser and host-side evaluators for lagrangian_mc tracer fields

#include <string>
#include <vector>

#include "athena.hpp"

namespace particles {

enum class TracerFieldKind {
  density, pressure, temperature, entropy, internal_energy,
  v1, v2, v3, vmag, sound_speed, mach, scalar,
  b1, b2, b3, bmag, magnetic_pressure, beta, alfven_speed,
  shock_valid, shock_region, shock_mach_pressure, shock_mach_temperature,
  shock_mach_pressure_valid, shock_mach_temperature_valid
};

KOKKOS_INLINE_FUNCTION
bool IsShockField(TracerFieldKind kind) {
  return kind>=TracerFieldKind::shock_valid && kind<=TracerFieldKind::shock_mach_temperature_valid;
}

// Signal and completed-step Observation share these evidence fields.
template<typename Sample> KOKKOS_INLINE_FUNCTION
Real EvaluateShockField(TracerFieldKind kind, const Sample &sample) {
  switch (kind) {
    case TracerFieldKind::shock_valid: return sample.valid;
    case TracerFieldKind::shock_region: return sample.inside;
    case TracerFieldKind::shock_mach_pressure: return sample.mach_pressure;
    case TracerFieldKind::shock_mach_temperature: return sample.mach_temperature;
    case TracerFieldKind::shock_mach_pressure_valid:
      return sample.valid && Kokkos::isfinite(sample.mach_pressure);
    case TracerFieldKind::shock_mach_temperature_valid:
      return sample.valid && Kokkos::isfinite(sample.mach_temperature);
    default: return 0;
  }
}

struct TracerField {
  TracerFieldKind kind = TracerFieldKind::density;
  int scalar_index = -1;
  std::string name = "density";
};

std::vector<TracerField> ParseTracerFieldList(const std::string &field_list,
                                              bool has_mhd, int nscalars,
                                              const std::string &context);

TracerField ParseTracerFieldName(const std::string &field_name, bool has_mhd,
                                 int nscalars, const std::string &context);

std::vector<std::string> TracerFieldNames(const std::vector<TracerField> &fields);

Real EvaluateTracerFieldHost(const TracerField &field, const HostArray5D<Real> &w0,
                             const HostArray5D<Real> &bcc, bool has_mhd,
                             Real gamma, int nfluid, int m, int k, int j, int i);

} // namespace particles

#endif // PARTICLES_TRACER_FIELDS_HPP_
