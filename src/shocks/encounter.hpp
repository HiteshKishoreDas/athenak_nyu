//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#ifndef SHOCKS_ENCOUNTER_HPP_
#define SHOCKS_ENCOUNTER_HPP_

#include <cstdint>
#include <limits>
#include "athena.hpp"

namespace shocks {
enum Phase { outside, active, pending_exit };
enum Flags : uint32_t {
  left_censored = 1, invalid_observation = 2, merged_flicker = 4,
  no_plateau = 8, cadence_limited = 16, temperature_invalid = 32,
  right_censored = 64, domain_exit = 128, direction_ambiguous = 256,
  downstream_reentry = 512
};
enum Status { completed = 0, partial = 1 };

// All quantities use Athena code units. No physical crossing time is inferred.
// Images count global periodic wraps; x is the wrapped observation coordinate.
struct Observation {
  Real time = -1, x[3] = {}, rho = 0, pressure = 0, temperature = 0;
  Real mach_pressure = 0, mach_temperature = 0, normal[3] = {};
  Real upstream[3] = {}, downstream[3] = {};
  uint64_t generation = 0;
  int64_t image[3] = {};
  bool valid = false, inside = false;
};
struct EncounterState {
  uint64_t particle_id = 0, sequence = 0;
  uint64_t last_generation = UINT64_MAX;
  int64_t image[3] = {};
  int phase = outside, clean = 0;
  bool suppressing_reentry = false;
  uint64_t suppressed_reentries = 0;
  uint32_t flags = cadence_limited;
  Observation previous, first, last, post;
};
struct Event {
  uint64_t particle_id = 0, sequence = 0;
  uint32_t status = completed, flags = 0;
  Observation pre, first, last, post, representative;
};

KOKKOS_INLINE_FUNCTION
void FillEvent(const EncounterState &s, Event &e, uint32_t reason = 0) {
  e.particle_id = s.particle_id;
  e.sequence = s.sequence;
  e.flags = s.flags | reason | no_plateau;
  e.status = (e.flags & (left_censored | right_censored | domain_exit)) ? partial
                                                                    : completed;
  e.pre = s.previous;
  e.first = s.first;
  e.last = s.last;
  e.post = s.post;
  e.representative = s.first;  // first inside sample, never a peak Mach selector
}

// Two successive valid outside observations rearm an episode. Reentry sooner
// merges it and records ambiguity. Fixed endpoints are not certified plateaus.
// Retry of a generation is idempotent. Invalid data never closes an encounter.
KOKKOS_INLINE_FUNCTION
bool UpdateEncounter(EncounterState &s, const Observation &o, Event &e) {
  if (s.last_generation == o.generation) return false;
  s.last_generation = o.generation;
  if (!o.valid) {
    s.flags |= invalid_observation;
    s.clean = 0;
    return false;
  }
  if (o.inside) {
    if (s.phase == outside) {
      // A valid previous material state must be closer to this shock's
      // simultaneous upstream endpoint than to its downstream endpoint.
      // This ordering has no fitted tolerance. Downstream numerical lobes
      // are retained as suppressed reentry evidence, not another passage.
      if (s.previous.valid && o.upstream[0] > 0 && o.upstream[1] > 0 &&
          o.downstream[0] > 0 && o.downstream[1] > 0) {
        Real up = Kokkos::fabs(Kokkos::log(s.previous.rho)-Kokkos::log(o.upstream[0])) +
                  Kokkos::fabs(Kokkos::log(s.previous.pressure)-
                               Kokkos::log(o.upstream[1]));
        Real down = Kokkos::fabs(Kokkos::log(s.previous.rho)-
                                 Kokkos::log(o.downstream[0])) +
                    Kokkos::fabs(Kokkos::log(s.previous.pressure)-
                                 Kokkos::log(o.downstream[1]));
        if (!(up < down)) {
          s.flags |= downstream_reentry;
          if (!s.suppressing_reentry) {
            if (s.suppressed_reentries == UINT64_MAX) {
              s.flags |= invalid_observation;
            } else {
              ++s.suppressed_reentries;
            }
          }
          s.suppressing_reentry = true;
          return false;
        }
      }
      s.suppressing_reentry = false;
      s.flags &= ~downstream_reentry;
      s.flags |= s.previous.valid ? 0u : left_censored;
      s.first = o;
    } else if (s.phase == pending_exit) {
      s.flags |= merged_flicker;
    }
    if (s.phase != outside && s.last.valid) {
      Real dot = 0;
      for (int d=0; d<3; ++d) dot += s.last.normal[d]*o.normal[d];
      if (!(dot > 0)) s.flags |= direction_ambiguous;
    }
    s.phase = active;
    s.clean = 0;
    s.last = o;
    return false;
  }
  if (s.phase == outside) {
    s.suppressing_reentry = false;
    s.previous = o;
    return false;
  }
  if (s.phase == active) {
    s.phase = pending_exit;
    s.post = o;  // preserve first outside, not the later rearming observation
    s.clean = 1;
    return false;
  }
  if (++s.clean < 2) return false;
  FillEvent(s, e);
  // Caller checks sequence saturation before accepting the observation.
  ++s.sequence;
  s.phase = outside;
  s.clean = 0;
  s.flags = cadence_limited;
  s.previous = o;
  s.first = Observation{};
  s.last = Observation{};
  s.post = Observation{};
  return true;
}
}  // namespace shocks
#endif  // SHOCKS_ENCOUNTER_HPP_
