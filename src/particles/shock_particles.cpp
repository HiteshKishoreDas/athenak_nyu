//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "eos/eos.hpp"
#include "srcterms/srcterms.hpp"
#include "srcterms/turb_driver.hpp"
#include "pgen/pgen.hpp"
#include "particles.hpp"
#include "shock_version.hpp"
#include "shocks/store.hpp"
#include "outputs/tracer_binary.hpp"

namespace {
using shocks::Digest;
using shocks::FileDigest;
using shocks::RankName;
using shocks::CommitFile;
using shocks::PutObservation;
using shocks::PutEvent;
[[noreturn]] void ShockError(const std::string &message) {shocks::Fail(message);}

struct CheckpointHeader {
  char magic[8] = {'A','K','S','H','O','C','K','2'};
  uint32_t version = 2, endian = 0x12345678, real_bytes = sizeof(Real);
  uint32_t state_bytes = sizeof(shocks::EncounterState), int_bytes = sizeof(int);
  int rank = 0, ranks = 0, gids = 0, nmb = 0, nrdata = 0, nidata = 0;
  uint64_t generation = 0, count = 0, contract_bytes = 0;
  uint64_t fluid_digest = 0, history_digest = 0;
  uint64_t sampled = 0, invalid = 0, complete = 0, partial = 0;
  int highwater = 0;
  Real time = 0;
};
static_assert(std::is_trivially_copyable<shocks::EncounterState>::value,
              "Migration/checkpoint state must have no pointers");

template<class T>
void AppendBytes(std::string &data, const T *p, size_t count) {
  data.append(reinterpret_cast<const char*>(p), sizeof(T)*count);
}
}  // namespace

namespace particles {
void Particles::ConfigureShockTracking(ParameterInput *pin) {
  if (pin->DoesBlockExist("shock_tracking")) {
    shock_enabled = pin->GetOrAddBoolean("shock_tracking", "enabled", false);
    if (!shock_enabled &&
        pin->GetOrAddBoolean("shock_tracking", "history_present", false)) {
      ShockError("Shock history cannot be disabled on restart");
    }
  }
  if (shock_enabled) {
    if (pusher != ParticlesPusher::lagrangian_tracer || !IsLagrangianTracer()) {
      ShockError("Shock tracking requires particle_type=lagrangian_tracer and "
                 "pusher=lagrangian_tracer; MC histories are unsupported");
    }
    if (pin->GetOrAddString("time", "integrator", "rk2") != "rk2" ||
        (pin->DoesBlockExist("cooling") &&
         pin->GetOrAddBoolean("cooling", "enabled", false))) {
      ShockError("Shock Lagrangian qualification requires RK2 and cooling disabled");
    }
    for (auto block=pin->block.begin(); block != pin->block.end(); ++block) {
      if (block->block_name.rfind("output",0) == 0 &&
          pin->DoesParameterExist(block->block_name,"file_type") &&
          pin->GetString(block->block_name,"file_type") == "rst" &&
          pin->GetOrAddBoolean(block->block_name,"single_file_per_node",false)) {
        ShockError("Shock Lagrangian restart does not support per-node files");
      }
    }
    std::set<int> schedule_ids;
    for (const auto &seed : seed_schedules_) {
      if (!schedule_ids.insert(seed.id).second) {
        ShockError("Shock restart requires unique tracer seed schedule IDs");
      }
    }
    auto &c = shock_config;
    c.width = pin->GetOrAddInteger("shock_tracking", "width", 1);
    c.min_mach = pin->GetOrAddReal("shock_tracking", "min_mach", 1.1);
    c.mach_max = pin->GetOrAddReal("shock_tracking", "mach_max", 1000);
    c.grad_floor = pin->GetOrAddReal("shock_tracking", "grad_floor", 1e-30);
    c.div_limit = pin->GetOrAddReal("shock_tracking", "div_limit", 0);
    c.density_floor = pin->GetOrAddReal("shock_tracking", "density_floor", 0);
    c.temperature_floor = pin->GetOrAddReal("shock_tracking", "temperature_floor", 0);
    shock_capacity = pin->GetOrAddInteger("shock_tracking", "event_capacity", 4096);
    shock_lineage = pin->GetString("shock_tracking", "lineage");
    shock_directory = pin->GetOrAddString("shock_tracking", "directory", "shock_events");
    shock_trace = pin->GetOrAddBoolean("shock_tracking", "trace", false);
    shock_binary = tracer_output::Binary(pin);
    bool ok = c.width >= 1 && c.width <= pmy_pack->pmesh->mb_indcs.ng &&
              shock_capacity > 0 && !shock_lineage.empty() &&
              shock_lineage.size() <= 128 && !shock_directory.empty();
    for (char ch : shock_lineage) {
      ok = ok && ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ||
                  (ch >= '0' && ch <= '9') || ch == '-' || ch == '_');
    }
    ok = ok && shocks::ConfigValid(c) && c.div_limit <= 0;
    if (!ok) {
      ShockError("Invalid shock_tracking configuration: finite floors and "
                 "nonpositive divergence threshold, "
                 "width in [1,nghost], positive capacity and safe lineage required");
    }
    c.gamma = pmy_pack->phydro ? pmy_pack->phydro->peos->eos_data.gamma : 0;
    std::ostringstream contract;
    contract.precision(std::numeric_limits<Real>::max_digits10);
    contract << "v1-entry-upstream-logdistance-velocity-multilinear-ssprk2-temperature-"
             << "nearest-reference6f5ec8bc "
             << SHOCK_IMPLEMENTATION_HASH << ' ' << sizeof(Real) << ' '
             << c.gamma << ' ' << c.width << ' ' << c.min_mach
             << ' ' << c.mach_max << ' ' << c.grad_floor << ' ' << c.div_limit
             << ' ' << c.density_floor << ' ' << c.temperature_floor << ' '
             << pin->GetOrAddString("time", "integrator", "rk2") << ' '
             << pin->GetOrAddString("hydro", "reconstruct", "plm") << ' '
             << pin->GetOrAddString("hydro", "rsolver", "llf") << ' '
             << pin->GetOrAddReal("time", "cfl_number", 0.3);
    auto &ms = pmy_pack->pmesh->mesh_size;
    auto &mi = pmy_pack->pmesh->mesh_indcs;
    auto &bi = pmy_pack->pmesh->mb_indcs;
    contract << ' ' << mi.nx1 << ' ' << mi.nx2 << ' ' << mi.nx3 << ' '
             << bi.nx1 << ' ' << bi.nx2 << ' ' << bi.nx3 << ' '
             << ms.x1min << ' ' << ms.x1max << ' ' << ms.x2min << ' ' << ms.x2max
             << ' ' << ms.x3min << ' ' << ms.x3max;
    for (auto bc : pmy_pack->pmesh->mesh_bcs) contract << ' ' << static_cast<int>(bc);
    if (pmy_pack->phydro) {
      auto eos = pmy_pack->phydro->peos->eos_data;
      const char *names[] = {"dfloor","pfloor","tfloor","sfloor"};
      int index = 0;
      for (Real floor : {eos.dfloor,eos.pfloor,eos.tfloor,eos.sfloor}) {
        if (!std::isfinite(floor) || floor < 0) {
          ShockError("Shock Lagrangian requires finite nonnegative native EOS floors");
        }
        pin->SetRealPrecise("hydro",names[index++],floor);
        contract << ' ' << floor;
      }
      contract << ' ' << pmy_pack->phydro->use_fofc
               << ' ' << pmy_pack->phydro->nscalars;
    }
    contract << " rng " << random_seed << " schedules " << seed_schedules_.size();
    for (const auto &seed : seed_schedules_) {
      // Native GetOrAddReal formats defaults with limited precision. Persist
      // resolved immutable values exactly so authentic restarts retain the contract.
      for (const auto &entry : std::initializer_list<std::pair<const char*,Real>>{
          {"start_time",seed.start_time},{"end_time",seed.end_time},
          {"cadence",seed.cadence},{"x1min",seed.x1min},{"x1max",seed.x1max},
          {"x2min",seed.x2min},{"x2max",seed.x2max},{"x3min",seed.x3min},
          {"x3max",seed.x3max},{"center1",seed.center1},{"center2",seed.center2},
          {"center3",seed.center3},{"radius",seed.radius},
          {"slab_min",seed.slab_min},{"slab_max",seed.slab_max}}) {
        if (!std::isfinite(entry.second)) {
          ShockError("Shock tracer schedule contains nonfinite values");
        }
        pin->SetRealPrecise(seed.block_name,entry.first,entry.second);
      }
      if ((seed.has_target_min && !std::isfinite(seed.target_min)) ||
          (seed.has_target_max && !std::isfinite(seed.target_max))) {
        ShockError("Shock tracer target contains nonfinite thresholds");
      }
      if (seed.has_target_min) {
        pin->SetRealPrecise(seed.block_name,"target_min",seed.target_min);
      }
      if (seed.has_target_max) {
        pin->SetRealPrecise(seed.block_name,"target_max",seed.target_max);
      }
      contract << ' ' << seed.id << ' ' << seed.start_time << ' ' << seed.end_time
               << ' ' << seed.cadence << ' ' << seed.count_per_event << ' ' << seed.seed
               << ' ' << static_cast<int>(seed.weight)
               << ' ' << static_cast<int>(seed.region)
               << ' ' << seed.x1min << ' ' << seed.x1max
               << ' ' << seed.x2min << ' ' << seed.x2max
               << ' ' << seed.x3min << ' ' << seed.x3max
               << ' ' << seed.center1 << ' ' << seed.center2 << ' ' << seed.center3
               << ' ' << seed.radius << ' ' << seed.slab_axis << ' ' << seed.slab_min
               << ' ' << seed.slab_max << ' ' << seed.has_target
               << ' ' << seed.target_field.name << ' ' << seed.has_target_min
               << ' ' << seed.has_target_max << ' ' << seed.target_min
               << ' ' << seed.target_max;
    }
    if (pmy_pack->pturb) {
      // Bind the complete resolved static forcing configuration, including RNG seed.
      // Modal amplitudes/RNG progress already persist in the native fluid restart.
      std::vector<std::pair<std::string,std::string>> forcing;
      for (const auto &block : pin->block) {
        if (block.block_name != "turb_driving") continue;
        for (const auto &line : block.line) {
          forcing.emplace_back(line.param_name,line.param_value);
        }
      }
      std::sort(forcing.begin(),forcing.end());
      contract << " turb-driver-config-v1 " << forcing.size();
      for (const auto &entry : forcing) {
        contract << ' ' << entry.first.size() << ':' << entry.first
                 << ' ' << entry.second.size() << ':' << entry.second;
      }
    }
    // Regrid controls affect subsequent trajectory sampling and must persist.
    for (const auto &block : pin->block) {
      if (block.block_name != "mesh_refinement" &&
          block.block_name.rfind("amr_criterion",0) != 0) continue;
      std::vector<std::pair<std::string,std::string>> settings;
      for (const auto &line : block.line) {
        settings.emplace_back(line.param_name,line.param_value);
      }
      std::sort(settings.begin(),settings.end());
      contract << ' ' << block.block_name << ' ' << settings.size();
      for (const auto &entry : settings) {
        contract << ' ' << entry.first.size() << ':' << entry.first
                 << ' ' << entry.second.size() << ':' << entry.second;
      }
    }
    for (const char *name : {"refine_time","coarsen_time","refine_x1max"}) {
      Real value = pin->GetOrAddReal("problem",name,1.0e30);
      if (!std::isfinite(value)) ShockError("nonfinite timed AMR callback parameter");
      pin->SetRealPrecise("problem",name,value);
      contract << ' ' << name << ' ' << value;
    }
    if (shock_binary) contract << " tracer-binary-v1";
    shock_contract = contract.str();
    if (shock_contract.size()+shock_lineage.size()+1 > 16384) {
      ShockError("Shock restart contract exceeds supported metadata size");
    }
    Kokkos::realloc(shock_state, nprtcl_thispack);
    Kokkos::realloc(shock_events, shock_capacity);
    const auto &bi2 = pmy_pack->pmesh->mb_indcs;
    Kokkos::realloc(shock_signal, pmy_pack->nmb_thispack,
                    bi2.nx3 > 1 ? bi2.nx3 + 2*bi2.ng : 1,
                    bi2.nx2 + 2*bi2.ng, bi2.nx1 + 2*bi2.ng);
    pin->SetBoolean("shock_tracking", "history_present", true);
  }
  if (shock_enabled) ValidateShockPhysics();
}

void Particles::ValidateShockPhysics() {
  if (std::string(DevExeSpace::name()) != "Serial" || sizeof(Real) != 8) {
    ShockError("qualified development configuration requires Serial CPU/double; "
               "threaded/GPU/single precision need separate runtime qualification");
  }
  auto pm = pmy_pack->pmesh;
  auto h = pmy_pack->phydro;
  auto coord = pmy_pack->pcoord;
  bool ok = h && h->peos->eos_data.is_ideal && h->peos->eos_data.use_e &&
            std::isfinite(h->peos->eos_data.gamma) && h->peos->eos_data.gamma > 1;
  ok = ok && !pm->one_d &&
       pm->nmb_packs_thisrank == 1 && !pmy_pack->pmhd && !pmy_pack->pionn &&
       !pmy_pack->prad && !pmy_pack->pz4c && !pmy_pack->pdyngr && !pmy_pack->padm &&
       !pmy_pack->pframe_tracker && IsLagrangianTracer() &&
       !coord->is_special_relativistic &&
       !coord->is_general_relativistic && !coord->is_dynamical_relativistic;
  if (pmy_pack->pturb) {
    ok = ok && pmy_pack->pturb->RestartMetadata().StaticFieldsFinite();
  }
  if (h) {
    auto s = h->psrc;
    ok = ok && !h->pvisc && !h->pcond && !h->porb_u && !h->psbox_u &&
         !h->use_fofc;
    if (s) ok = ok && !s->const_accel && !s->rel_cooling && !s->rad_beam;
  }
  if (pm->pgen) ok = ok && !pm->pgen->user_srcs && !pm->pgen->user_srcs_func;
  for (int d=0; d<(pm->three_d ? 6 : 4); ++d) {
    ok = ok && (pm->mesh_bcs[d] == BoundaryFlag::periodic ||
                pm->mesh_bcs[d] == BoundaryFlag::outflow);
  }
  auto &b = pm->mb_indcs;
  ok = ok && b.ng >= 1 && b.nx1 >= 4 && b.nx2 >= 4 &&
       (!pm->three_d || b.nx3 >= 4);
  if (shock_enabled) ok = ok && shocks::ConfigValid(shock_config);
  if (!ok) ShockError("requires 2D/3D Newtonian ideal Cartesian hydro, "
                      "periodic/outflow boundaries, no other sources/transport, one pack; "
                      "native forcing parameters must be finite");
}

void Particles::ResizeShockState(int nold) {
  if (!shock_enabled) return;
  Kokkos::resize(shock_state, nprtcl_thispack);
  auto s = shock_state;
  auto pi = prtcl_idata;
  par_for("initialize_encounters", DevExeSpace(), nold, nprtcl_thispack-1,
  KOKKOS_LAMBDA(int p) {
    s(p) = shocks::EncounterState{};
    s(p).particle_id = static_cast<uint64_t>(pi(PTAG,p));
  });
}

void Particles::PrepareShockStore(bool restarting) {
  try {
    std::filesystem::create_directories(shock_directory);
  } catch (const std::exception &e) {
    ShockError(std::string("cannot create event directory: ") + e.what());
  }
  std::string manifest = "AKSHOCK-EVENTS 1\n" + shock_lineage + "\n" +
                         shock_contract + "\n";
  CommitFile(shock_directory + "/manifest" + RankName(), manifest, restarting);
  // FNV-1a chained checksum covers accidental corruption, not hostile collisions.
  // Every generation has an immutable chunk, including empty generations.
  // Missing committed history cannot silently become a fresh output stream.
  if (restarting && shock_generation != UINT64_MAX) {
    uint64_t digest = UINT64_C(14695981039346656037);
    for (uint64_t g=0; g<=shock_generation; ++g) {
      const std::string path = shock_directory + "/cycle" + std::to_string(g) +
                               RankName() + tracer_output::Extension(shock_binary);
      if (!std::filesystem::exists(path)) {
        ShockError("restart requires complete event chunks through "
                   "checkpoint generation");
      }
      std::ifstream in(path, std::ios::binary);
      std::string data((std::istreambuf_iterator<char>(in)), {});
      if (!in || in.bad()) ShockError("cannot validate event chunk: " + path);
      digest = Digest(digest,data.data(),data.size());
    }
    if (digest != shock_history_digest) ShockError("committed event history corrupted");
  }
  if (global_variable::my_rank == 0) {
    std::cout << "Experimental shock history: " << shock_contract << "\nlineage="
              << shock_lineage << " capacity=" << shock_capacity
              << "; first inside representative; code units; no certified plateaus; "
              << "endpoint cadence limited; scientific review pending" << std::endl;
  }
}

void Particles::ObserveShocks(Real time, uint64_t generation) {
  if (!shock_enabled || generation == shock_generation) return;
  if (shock_generation != UINT64_MAX && generation != shock_generation+1) {
    ShockError("diagnostic generation gap");
  }
  Kokkos::Timer timer;
  auto pm = pmy_pack->pmesh;
  auto b = pm->mb_indcs;
  auto size = pmy_pack->pmb->mb_size.d_view;
  auto w = pmy_pack->phydro->w0;
  auto config = shock_config;
  const int nk = b.nx3 > 1 ? b.nx3+2*b.ng : 1;
  const int nj = b.nx2+2*b.ng, ni = b.nx1+2*b.ng;
  if (shock_signal.extent(0) != static_cast<size_t>(pmy_pack->nmb_thispack) ||
      shock_signal.extent(1) != static_cast<size_t>(nk) ||
      shock_signal.extent(2) != static_cast<size_t>(nj) ||
      shock_signal.extent(3) != static_cast<size_t>(ni)) {
    Kokkos::realloc(shock_signal,pmy_pack->nmb_thispack,nk,nj,ni);
  }
  auto signal = shock_signal;
  const bool td = pm->three_d;
  bool physical[6];
  for (int d=0; d<6; ++d) physical[d] = pm->mesh_bcs[d] != BoundaryFlag::periodic;
  auto ms = pm->mesh_size;
  Real low[3] = {ms.x1min,ms.x2min,ms.x3min};
  Real high[3] = {ms.x1max,ms.x2max,ms.x3max};
  par_for("live_unreduced_shock_detector", DevExeSpace(), 0,
          pmy_pack->nmb_thispack-1, b.ks,b.ke,b.js,b.je,b.is,b.ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    auto s = shocks::DetectCell(w,m,k,j,i,size(m).dx1,size(m).dx2,size(m).dx3,td,config);
    int axis = 0;
    for (int d=1; d<(td ? 3 : 2); ++d) {
      if (Kokkos::fabs(s.normal[d]) > Kokkos::fabs(s.normal[axis])) axis = d;
    }
    Real position[3] = {size(m).x1min+(i-b.is+0.5)*size(m).dx1,
                        size(m).x2min+(j-b.js+0.5)*size(m).dx2,
                        size(m).x3min+(k-b.ks+0.5)*size(m).dx3};
    Real spacing[3] = {size(m).dx1,size(m).dx2,size(m).dx3};
    if ((physical[2*axis] && position[axis]-config.width*spacing[axis] < low[axis]) ||
        (physical[2*axis+1] && position[axis]+config.width*spacing[axis] >= high[axis])) {
      s.valid = false;  // extrapolated outflow ghosts are not physical endpoints
      s.inside = false;
    }
    signal(m,k,j,i) = s;
  });
  auto state = shock_state;
  const bool trace = shock_trace;
  DvceArray1D<shocks::Observation> observations;
  if (trace) Kokkos::realloc(observations,nprtcl_thispack);
  auto events = shock_events;
  auto pr = prtcl_rdata;
  auto pi = prtcl_idata;
  const int gids = pmy_pack->gids, capacity = shock_capacity;
  Kokkos::View<int, DevMemSpace> count("shock_event_count"), bad("shock_state_error");
  Kokkos::View<int, DevMemSpace> invalid("shock_invalid_samples");
  par_for("observe_fluid_encounters", DevExeSpace(), 0, nprtcl_thispack-1,
  KOKKOS_LAMBDA(int p) {
    if (state(p).sequence == UINT64_MAX || state(p).suppressed_reentries == UINT64_MAX) {
      Kokkos::atomic_inc(&bad());
      return;
    }
    int m = pi(PGID,p)-gids;
    int i = static_cast<int>(Kokkos::floor((pr(IPX,p)-size(m).x1min)/size(m).dx1));
    int j = static_cast<int>(Kokkos::floor((pr(IPY,p)-size(m).x2min)/size(m).dx2));
    int k = td ? static_cast<int>(Kokkos::floor(
                                 (pr(IPZ,p)-size(m).x3min)/size(m).dx3)) : 0;
    shocks::Observation o;
    o.time = time;
    o.generation = generation;
    for (int d=0; d<3; ++d) {
      o.x[d] = pr(2*d,p);
      o.image[d] = state(p).image[d];
    }
    if (i >= 0 && i < b.nx1 && j >= 0 && j < b.nx2 && k >= 0 && k < b.nx3) {
      auto s = signal(m,k+b.ks,j+b.js,i+b.is);
      o.valid = s.valid;
      o.inside = s.inside;
      o.rho = s.rho;
      o.pressure = s.pressure;
      o.temperature = s.temperature;
      o.mach_pressure = s.mach_pressure;
      o.mach_temperature = s.mach_temperature;
      for (int d=0; d<3; ++d) o.normal[d] = s.normal[d];
      o.upstream[0] = s.rho_up;
      o.upstream[1] = s.pressure_up;
      o.upstream[2] = s.temperature_up;
      o.downstream[0] = s.rho_down;
      o.downstream[1] = s.pressure_down;
      o.downstream[2] = s.temperature_down;
      if (s.inside && !Kokkos::isfinite(s.mach_temperature)) {
        state(p).flags |= shocks::temperature_invalid;
      }
    }
    if (!o.valid) Kokkos::atomic_inc(&invalid());
    if (trace) observations(p) = o;
    shocks::Event e;
    if (shocks::UpdateEncounter(state(p), o, e)) {
      int slot = Kokkos::atomic_fetch_add(&count(), 1);
      if (slot < capacity) events(slot) = e;
    }
  });
  auto hc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), count);
  auto hb = Kokkos::create_mirror_view_and_copy(HostMemSpace(), bad);
  auto hi = Kokkos::create_mirror_view_and_copy(HostMemSpace(), invalid);
  if (shock_sampled > UINT64_MAX-static_cast<uint64_t>(nprtcl_thispack) ||
      shock_invalid > UINT64_MAX-hi()) ShockError("diagnostic counter exhausted");
  shock_sampled += nprtcl_thispack;
  shock_invalid += hi();
  if (hb()) ShockError("event sequence exhausted");
  if (hc() > capacity || shock_pending_exits.size() + hc() >
                         static_cast<size_t>(capacity)) {
    ShockError("event buffer overflow; no checkpoint written, replay prior checkpoint "
               "with larger event_capacity");
  }
  auto he = Kokkos::create_mirror_view_and_copy(HostMemSpace(), events);
  for (int n=0; n<hc(); ++n) shock_pending_exits.push_back(he(n));
  shock_highwater = std::max(shock_highwater,
                             static_cast<int>(shock_pending_exits.size()));
  for (const auto &e : shock_pending_exits) {
    if (e.status == shocks::completed) {
      if (shock_complete == UINT64_MAX) ShockError("completed counter exhausted");
      ++shock_complete;
    } else {
      if (shock_partial == UINT64_MAX) ShockError("partial counter exhausted");
      ++shock_partial;
    }
  }
  DrainShockEvents(shock_pending_exits, generation);
  shock_pending_exits.clear();
  shock_generation = generation;
  shock_time = time;
  if (shock_trace) {
    auto hs = Kokkos::create_mirror_view_and_copy(HostMemSpace(), state);
    auto ho = Kokkos::create_mirror_view_and_copy(HostMemSpace(), observations);
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<Real>::max_digits10);
    if (!shock_binary) for (int p=0; p<nprtcl_thispack; ++p) {
      out << hs(p).particle_id << ',' << hs(p).sequence << ',' << hs(p).phase << ','
          << hs(p).flags << ',';
      PutObservation(out, ho(p));
      out << '\n';
    }
    std::string data=out.str();
    if (shock_binary) {
      tracer_output::Table table("particle_trace",tracer_output::TraceColumns(false),
          nprtcl_thispack,generation,time,"position=wrapped\nlineage="+shock_lineage+
          "\ncontract="+shock_contract+"\n");
      for (int p=0; p<nprtcl_thispack; ++p) {
        table.Put(hs(p).particle_id); table.Put(hs(p).sequence);
        table.Put(static_cast<uint32_t>(hs(p).phase)); table.Put(hs(p).flags);
        tracer_output::Observation(table,ho(p));
      }
      data=table.Finish();
    }
    CommitFile(shock_directory + "/trace" + std::to_string(generation) + RankName() +
                 tracer_output::Extension(shock_binary),data,true);
  }
  shock_seconds += timer.seconds();
}

void Particles::DrainShockEvents(const std::vector<shocks::Event> &events,
                                uint64_t generation) {
  auto ordered = events;
  std::sort(ordered.begin(), ordered.end(), [](const auto &a, const auto &b) {
    return a.particle_id < b.particle_id ||
           (a.particle_id == b.particle_id && a.sequence < b.sequence);
  });
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<Real>::max_digits10)
      << "# AKSHOCK-EVENTS 1 " << shock_lineage << ' ' << generation << ' '
      << ordered.size() << '\n' << "# " << shock_contract << '\n';
  if (!shock_binary) for (const auto &e : ordered) PutEvent(out,e);
  out << "# END " << ordered.size() << '\n';
  std::string data=out.str();
  if (shock_binary) {
    tracer_output::Table table("particle_events",tracer_output::EventColumns(false),
        ordered.size(),generation,pmy_pack->pmesh->time,
        "position=wrapped\nlineage="+shock_lineage+"\ncontract="+shock_contract+"\n");
    for (const auto &e : ordered) tracer_output::Event(table,e);
    data=table.Finish();
  }
  CommitFile(shock_directory + "/cycle" + std::to_string(generation) + RankName() +
               tracer_output::Extension(shock_binary), data, true);
  shock_history_digest = Digest(shock_history_digest,data.data(),data.size());
}

void Particles::EmitShockExits(const std::vector<int> &indices) {
  if (!shock_enabled) return;
  auto s = Kokkos::create_mirror_view_and_copy(HostMemSpace(), shock_state);
  auto pr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  for (int p : indices) {
    if (p < 0 || p >= nprtcl_thispack) ShockError("invalid destruction index");
    if (s(p).phase == shocks::outside) continue;
    if (s(p).sequence == UINT64_MAX || shock_pending_exits.size() >=
                                       static_cast<size_t>(shock_capacity)) {
      ShockError("domain-exit event capacity/sequence exhausted");
    }
    shocks::Event e;
    shocks::FillEvent(s(p), e, shocks::right_censored | shocks::domain_exit);
    // The endpoint lies outside the physical domain: retain its time/position,
    // with invalid fluid evidence, rather than inventing an exit crossing time.
    e.post = shocks::Observation{};
    e.post.time = pmy_pack->pmesh->time+pmy_pack->pmesh->dt;
    e.post.generation = pmy_pack->pmesh->ncycle+1;
    for (int d=0; d<3; ++d) {
      e.post.x[d] = pr(2*d,p);
      e.post.image[d] = s(p).image[d];
    }
    shock_pending_exits.push_back(e);
  }
}

void Particles::PreflightShockRestart(const std::string &fname) {
  int exists = std::filesystem::exists(fname) ||
               std::filesystem::exists(fname + ".shock" + RankName());
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE,&exists,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
  if (exists) {
    ShockError("existing paired checkpoint: " + fname +
               "; restart to a new job/basename before writing checkpoints");
  }
}

void Particles::WriteShockRestart(const std::string &fname) {
  if (!shock_pending_exits.empty() || shock_generation !=
      static_cast<uint64_t>(pmy_pack->pmesh->ncycle) ||
      shock_time != pmy_pack->pmesh->time) {
    ShockError("checkpoint requires drained synchronized final-step state");
  }
  CheckpointHeader h;
  h.rank = global_variable::my_rank;
  h.ranks = global_variable::nranks;
  h.gids = pmy_pack->gids;
  h.nmb = pmy_pack->nmb_thispack;
  h.nrdata = nrdata;
  h.nidata = nidata;
  h.generation = shock_generation;
  h.time = shock_time;
  h.count = nprtcl_thispack;
  h.history_digest = shock_history_digest;
  h.sampled = shock_sampled;
  h.invalid = shock_invalid;
  h.complete = shock_complete;
  h.partial = shock_partial;
  h.highwater = shock_highwater;
  int fd = open(fname.c_str(),O_RDONLY);
  if (fd < 0 || fsync(fd) || close(fd)) ShockError("fluid checkpoint flush failed");
  h.fluid_digest = FileDigest(fname);
  const std::string contract = shock_lineage + "\n" + shock_contract;
  h.contract_bytes = contract.size();
  auto s = Kokkos::create_mirror_view_and_copy(HostMemSpace(), shock_state);
  std::string data;
  AppendBytes(data, &h, 1);
  data += contract;
  if (nprtcl_thispack > 0) AppendBytes(data,s.data(),nprtcl_thispack);
  uint64_t checksum = Digest(UINT64_C(14695981039346656037),data.data(),data.size());
  AppendBytes(data,&checksum,1);
  CommitFile(fname + ".shock" + RankName(), data, false);
}

void Particles::ReadShockRestart(const std::string &fname) {
  const std::string path = fname + ".shock" + RankName();
  std::ifstream in(path, std::ios::binary);
  CheckpointHeader h;
  in.read(reinterpret_cast<char*>(&h), sizeof(h));
  if (!in || std::memcmp(h.magic,"AKSHOCK2",8) || h.version != 2 ||
      h.endian != 0x12345678 || h.real_bytes != sizeof(Real) ||
      h.state_bytes != sizeof(shocks::EncounterState) || h.int_bytes != sizeof(int) ||
      h.rank != global_variable::my_rank || h.ranks != global_variable::nranks ||
      h.gids != pmy_pack->gids || h.nmb != pmy_pack->nmb_thispack ||
      h.nrdata != nrdata || h.nidata != nidata ||
      h.count != static_cast<uint64_t>(nprtcl_thispack) ||
      h.contract_bytes > 16384 || h.time != pmy_pack->pmesh->time ||
      h.generation != static_cast<uint64_t>(pmy_pack->pmesh->ncycle)) {
    ShockError("missing/incompatible shock sidecar: " + path);
  }
  uint64_t record_bytes = sizeof(shocks::EncounterState);
  const auto bytes = std::filesystem::file_size(path);
  if (bytes < sizeof(h) + h.contract_bytes + sizeof(uint64_t) ||
      bytes - sizeof(h) - h.contract_bytes - sizeof(uint64_t) != h.count*record_bytes) {
    ShockError("truncated or malformed shock sidecar: " + path);
  }
  in.seekg(bytes-sizeof(uint64_t));
  uint64_t stored_checksum;
  in.read(reinterpret_cast<char*>(&stored_checksum),sizeof(stored_checksum));
  if (!in || FileDigest(path,bytes-sizeof(uint64_t)) != stored_checksum) {
    ShockError("shock sidecar checksum mismatch");
  }
  in.seekg(sizeof(h));
  std::string contract(h.contract_bytes, '\0');
  in.read(contract.data(), contract.size());
  if (contract != shock_lineage + "\n" + shock_contract) {
    ShockError("restart lineage/detector/solver/grid contract changed");
  }
  if (FileDigest(fname) != h.fluid_digest) ShockError("fluid/sidecar checksum mismatch");
  Kokkos::realloc(shock_state, nprtcl_thispack);
  auto pr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_rdata);
  auto pi = Kokkos::create_mirror_view_and_copy(HostMemSpace(), prtcl_idata);
  auto s = Kokkos::create_mirror_view(shock_state);
  std::set<uint64_t> ids;
  for (int p=0; p<nprtcl_thispack; ++p) {
    in.read(reinterpret_cast<char*>(&s(p)), sizeof(shocks::EncounterState));
    bool valid_position = true;
    int m = pi(PGID,p)-h.gids;
    if (m >= 0 && m < h.nmb) {
      auto ms = pmy_pack->pmb->mb_size.h_view(m);
      Real low[3] = {ms.x1min,ms.x2min,ms.x3min};
      Real high[3] = {ms.x1max,ms.x2max,ms.x3max};
      for (int d=0; d<(pmy_pack->pmesh->three_d ? 3 : 2); ++d) {
        valid_position = valid_position && std::isfinite(pr(2*d,p)) &&
                         pr(2*d,p) >= low[d] && pr(2*d,p) < high[d];
      }
    } else {
      valid_position = false;
    }
    if (!in || !valid_position || pi(PGID,p) < h.gids || pi(PGID,p) >= h.gids + h.nmb ||
        pi(PTAG,p) < 0 || s(p).particle_id != static_cast<uint64_t>(pi(PTAG,p)) ||
        !ids.insert(s(p).particle_id).second || s(p).phase < shocks::outside ||
        s(p).phase > shocks::pending_exit || s(p).last_generation != h.generation) {
      ShockError("invalid particle state in sidecar");
    }
  }
  Kokkos::deep_copy(shock_state, s);
  shock_sampled = h.sampled;
  shock_invalid = h.invalid;
  shock_complete = h.complete;
  shock_partial = h.partial;
  shock_highwater = h.highwater;
  shock_history_digest = h.history_digest;
  shock_generation = h.generation;
  shock_time = h.time;
}

void Particles::WriteShockSuspended() {
  // Open episodes persist without consuming their identity. They are censored at
  // this run boundary, and can resume from the paired fluid checkpoint.
  auto s = Kokkos::create_mirror_view_and_copy(HostMemSpace(), shock_state);
  std::vector<shocks::Event> records;
  uint64_t suppressed = 0;
  for (int p=0; p<nprtcl_thispack; ++p) {
    if (suppressed > UINT64_MAX-s(p).suppressed_reentries) ShockError("counter overflow");
    suppressed += s(p).suppressed_reentries;
    if (s(p).phase != shocks::outside) {
      shocks::Event e;
      shocks::FillEvent(s(p),e,shocks::right_censored);
      if (records.size() >= static_cast<size_t>(shock_capacity)) {
        ShockError("suspended event buffer overflow; increase event_capacity");
      }
      records.push_back(e);
    }
  }
  std::sort(records.begin(),records.end(),[](const auto &a,const auto &b) {
    return a.particle_id < b.particle_id;
  });
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<Real>::max_digits10)
      << "# AKSHOCK-SUSPENDED 1 " << shock_lineage << ' ' << shock_generation << ' '
      << records.size() << '\n' << "# " << shock_contract << '\n';
  if (!shock_binary) for (const auto &e : records) PutEvent(out,e);
  out << "# END " << records.size() << '\n';
  std::string data=out.str();
  if (shock_binary) {
    tracer_output::Table table("particle_suspended",tracer_output::EventColumns(false),
        records.size(),shock_generation,shock_time,
        "position=wrapped\nlineage="+shock_lineage+"\ncontract="+shock_contract+"\n");
    for (const auto &e : records) tracer_output::Event(table,e);
    data=table.Finish();
  }
  CommitFile(shock_directory + "/suspended" + std::to_string(shock_generation) +
               RankName() + tracer_output::Extension(shock_binary), data, true);
  std::cout << "Shock rank=" << global_variable::my_rank << " sampled=" << shock_sampled
            << " invalid=" << shock_invalid << " complete=" << shock_complete
            << " partial=" << shock_partial << " suspended=" << records.size()
            << " suppressed_reentries_surviving=" << suppressed
            << " event_highwater=" << shock_highwater
            << " diagnostic_seconds_this_process=" << shock_seconds << std::endl;
}

}  // namespace particles
