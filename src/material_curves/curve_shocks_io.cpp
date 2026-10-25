// Copyright (C) 2026 the AthenaK contributors
// Licensed under the 3-clause BSD License (the "LICENSE").
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "material_curves/material_curves.hpp"
#include "shocks/store.hpp"

namespace material_curves {
namespace {
constexpr uint64_t kOffset = UINT64_C(14695981039346656037);
constexpr uint64_t kContractCap = 16384;

struct Header {
  char magic[8] = {'A','K','C','S','H','C','K','1'};
  uint32_t version = 1, endian = 0x12345678, real_bytes = sizeof(Real);
  uint32_t state_bytes = sizeof(ShockNode), int_bytes = sizeof(int);
  int rank = 0, ranks = 0;
  uint64_t generation = 0, count = 0, contract_bytes = 0;
  uint64_t fluid_bytes = 0, curve_bytes = 0;
  uint64_t fluid_digest = 0, curve_digest = 0, history_digest = 0;
  Real time = 0;
};
static_assert(std::is_trivially_copyable<ShockNode>::value,
              "Curve checkpoint state must have no pointers");

uint64_t Size(const std::string &path) {
  std::error_code error;
  const auto size = std::filesystem::file_size(path,error);
  if (error) shocks::Fail("cannot size curve shock input: " + path);
  return size;
}
std::string Sidecar(const std::string &fluid) {
  return fluid + ".curves.shock" + shocks::RankName();
}
std::string Manifest(const MaterialCurves &m) {
  return "AKCURVE-SHOCK-EVENTS 1\n" + m.shock_lineage + "\n" +
         m.shock_contract + "\n";
}
void CheckManifest(const MaterialCurves &m) {
  const auto expected = Manifest(m);
  const auto path = m.shock_directory + "/manifest" + shocks::RankName();
  if (Size(path) != expected.size()) shocks::Fail("curve shock manifest mismatch");
  std::ifstream in(path,std::ios::binary);
  std::string data(expected.size(),'\0');
  in.read(data.data(),data.size());
  if (!in || data != expected) shocks::Fail("curve shock manifest mismatch");
}
uint64_t PrefixDigest(const MaterialCurves &m, uint64_t generation) {
  CheckManifest(m);
  uint64_t digest = kOffset;
  for (uint64_t g=0; g<=generation; ++g) {
    const auto path = m.shock_directory + "/cycle" + std::to_string(g) +
                      shocks::RankName() + ".csv";
    if (m.shock_capacity <= 0 || static_cast<uint64_t>(m.shock_capacity) >
        (UINT64_MAX-kContractCap)/8192 ||
        Size(path) > kContractCap+static_cast<uint64_t>(m.shock_capacity)*8192) {
      shocks::Fail("oversized curve shock history chunk");
    }
    std::ifstream in(path,std::ios::binary);
    if (!in) shocks::Fail("missing curve shock history chunk: " + path);
    char buffer[65536];
    while (in.read(buffer,sizeof(buffer)) || in.gcount()) {
      digest = shocks::Digest(digest,buffer,in.gcount());
    }
    if (!in.eof()) shocks::Fail("cannot validate curve shock history: " + path);
  }
  return digest;
}
std::unordered_map<NodeID, const Node *> Geometry(const MaterialCurves &m) {
  std::unordered_map<NodeID, const Node *> nodes;
  for (const auto &curve : m.curves) {
    if (curve.id % global_variable::nranks !=
        static_cast<uint64_t>(global_variable::my_rank)) {
      shocks::Fail("curve shock checkpoint home-rank mismatch");
    }
    for (const auto &node : curve.nodes) {
      if (nodes.size() >= m.max_nodes || !nodes.emplace(node.id,&node).second) {
        shocks::Fail("curve shock geometry count/identity mismatch");
      }
    }
  }
  return nodes;
}
void CheckObservation(const shocks::Observation &o, Real birth,
                      uint64_t generation, Real time) {
  // Check raw bool representations before reading a native-ABI checkpoint bool.
  const auto bytes = reinterpret_cast<const unsigned char *>(&o);
  if (bytes[offsetof(shocks::Observation,valid)] > 1 ||
      bytes[offsetof(shocks::Observation,inside)] > 1 ||
      !std::isfinite(o.time) || (o.inside && !o.valid)) {
    shocks::Fail("invalid curve shock observation metadata");
  }
  for (Real x : o.x) {
    if (!std::isfinite(x)) shocks::Fail("nonfinite curve shock position");
  }
  if (o.time != -1 && (o.time < birth || o.time > time || o.generation > generation)) {
    shocks::Fail("curve shock observation outside checkpoint/birth interval");
  }
  if (!o.valid) return;
  if (o.time < birth || o.rho <= 0 || o.pressure <= 0 || o.temperature <= 0 ||
      !std::isfinite(o.rho) || !std::isfinite(o.pressure) ||
      !std::isfinite(o.temperature) || std::isinf(o.mach_pressure) ||
      std::isinf(o.mach_temperature)) shocks::Fail("invalid curve shock fluid evidence");
  for (int d=0; d<3; ++d) {
    if (!std::isfinite(o.normal[d]) || !std::isfinite(o.upstream[d]) ||
        !std::isfinite(o.downstream[d]) || o.upstream[d] <= 0 ||
        o.downstream[d] <= 0) shocks::Fail("invalid curve shock endpoint evidence");
  }
}
void CheckState(const ShockNode &node, const Node &geometry,
                uint64_t generation, Real time, bool three_d) {
  const auto &s = node.encounter;
  const auto bytes = reinterpret_cast<const unsigned char *>(&s);
  if (s.particle_id != IDBits(geometry.id) || !std::isfinite(node.birth) ||
      node.birth != geometry.birth || node.birth < 0 || node.birth > time ||
      s.last_generation != generation || s.phase < shocks::outside ||
      s.phase > shocks::pending_exit || s.clean < 0 || s.clean > 1 ||
      (s.phase != shocks::pending_exit && s.clean != 0) ||
      (s.flags & ~1023u) ||
      bytes[offsetof(shocks::EncounterState,suppressing_reentry)] > 1) {
    shocks::Fail("invalid curve shock state/identity/birth");
  }
  for (const auto *o : {&s.previous,&s.first,&s.last,&s.post}) {
    CheckObservation(*o,node.birth,generation,time);
  }
  if (s.phase != shocks::outside &&
      (!s.first.valid || !s.first.inside || !s.last.valid || !s.last.inside ||
       s.last.time < s.first.time ||
       (s.previous.valid && s.previous.time > s.first.time) ||
       (s.phase == shocks::active && s.post.valid) ||
       (s.phase == shocks::pending_exit && (!s.post.valid || s.post.inside ||
                                           s.post.time < s.last.time)))) {
    shocks::Fail("invalid curve shock encounter brackets");
  }
  for (int which=0; which<2; ++which) {
    const int level = which ? node.last_level : node.first_level;
    const Real *spacing = which ? node.last_spacing : node.first_spacing;
    if (level < -1 || (s.phase != shocks::outside && level < 0)) {
      shocks::Fail("invalid curve shock refinement level");
    }
    for (int d=0; d<3; ++d) {
      if (!std::isfinite(spacing[d]) || spacing[d] < 0 ||
          (s.phase != shocks::outside && d<(three_d ? 3 : 2) && spacing[d] <= 0)) {
        shocks::Fail("invalid curve shock cell spacing");
      }
    }
  }
}
std::string Records(const MaterialCurves &m, std::vector<ShockEvent> events,
                    const char *kind) {
  if (m.shock_capacity <= 0 || events.size() >
      static_cast<size_t>(m.shock_capacity)) shocks::Fail("curve shock event overflow");
  std::sort(events.begin(),events.end(),[](const auto &a,const auto &b) {
    return std::pair(a.encounter.particle_id,a.encounter.sequence) <
           std::pair(b.encounter.particle_id,b.encounter.sequence);
  });
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<Real>::max_digits10)
      << "# AKCURVE-SHOCK-" << kind << " 1 " << m.shock_lineage << ' '
      << m.shock_generation << ' ' << events.size() << '\n'
      << "# " << m.shock_contract << '\n';
  uint64_t previous_id = 0, previous_sequence = 0;
  for (const auto &e : events) {
    if (e.encounter.particle_id == 0 || (e.encounter.particle_id == previous_id &&
        e.encounter.sequence == previous_sequence)) shocks::Fail("duplicate curve event");
    previous_id = e.encounter.particle_id;
    previous_sequence = e.encounter.sequence;
    shocks::PutEvent(out,e.encounter,false,
                     SignedID(static_cast<std::uint32_t>(e.encounter.particle_id)));
    out << ',' << e.birth << ',' << e.first_level;
    for (Real dx : e.first_spacing) out << ',' << dx;
    out << ',' << e.last_level;
    for (Real dx : e.last_spacing) out << ',' << dx;
    out << '\n';
  }
  out << "# END " << events.size() << '\n';
  return out.str();
}
}  // namespace

void MaterialCurves::PrepareShockStore(bool restart) {
  if (!shock_enabled) return;
  if (restart) {
    if (shock_generation == UINT64_MAX ||
        PrefixDigest(*this,shock_generation) != shock_history_digest) {
      shocks::Fail("curve shock committed prefix missing or corrupted");
    }
  } else {
    std::error_code error;
    std::filesystem::create_directories(shock_directory,error);
    if (error) shocks::Fail("cannot create curve shock directory: " + error.message());
    shocks::CommitFile(shock_directory+"/manifest"+shocks::RankName(),
                       Manifest(*this),false);
  }
}
void MaterialCurves::WriteShockChunk(const std::vector<ShockEvent> &events) {
  if (!shock_enabled) return;
  const auto data = Records(*this,events,"EVENTS");
  shocks::CommitFile(shock_directory+"/cycle"+std::to_string(shock_generation)+
                     shocks::RankName()+".csv",data,true);
  shock_history_digest = shocks::Digest(shock_history_digest,data.data(),data.size());
}
void MaterialCurves::WriteShockTrace(const std::vector<Node *> &nodes,
                                    const std::vector<ShockObservation> &samples,
                                    bool births_only) {
  if (!shock_enabled || !shock_trace) return;
  if (nodes.size() != samples.size() || nodes.size() > max_nodes) {
    shocks::Fail("curve shock trace extent mismatch");
  }
  const auto generation = static_cast<uint64_t>(mesh->ncycle);
  std::vector<size_t> order(nodes.size());
  for (size_t i=0; i<order.size(); ++i) order[i] = i;
  std::sort(order.begin(),order.end(),[&](size_t a,size_t b) {
    return nodes[a]->id < nodes[b]->id;
  });
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<Real>::max_digits10);
  NodeID previous = 0;
  for (size_t i : order) {
    const auto &node = *nodes[i];
    auto found = shock_states.find(node.id);
    if (node.id == previous || found == shock_states.end() ||
        found->second.birth != node.birth ||
        samples[i].observation.generation != generation) {
      shocks::Fail("curve shock trace identity/birth/generation mismatch");
    }
    previous = node.id;
    const auto &s = found->second.encounter;
    out << node.id << ',' << s.sequence << ',' << s.phase << ',' << s.flags << ',';
    shocks::PutObservation(out,samples[i].observation);
    out << ',' << node.birth << ',' << samples[i].level;
    for (Real dx : samples[i].spacing) out << ',' << dx;
    out << '\n';
  }
  shocks::CommitFile(shock_directory+(births_only ? "/birthtrace" : "/trace")+
                     std::to_string(generation)+shocks::RankName()+".csv",
                     out.str(),true);
}
void MaterialCurves::PreflightShockRestart(const std::string &fluid_path) const {
  if (!shock_enabled) return;
  if (global_variable::nranks>1 && mesh->pmb_pack->ppart) {
    shocks::Fail("MC/curve multi-rank checkpoints lack a compatible shared layout");
  }
  int exists = std::filesystem::exists(fluid_path) ||
               std::filesystem::exists(fluid_path+".curves") ||
               std::filesystem::exists(Sidecar(fluid_path));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE,&exists,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
  if (exists) shocks::Fail("existing paired curve shock checkpoint; use a new basename");
}
void MaterialCurves::SaveShockRestart(const std::string &fluid_path) const {
  if (!shock_enabled) return;
  if (shock_generation != static_cast<uint64_t>(mesh->ncycle) ||
      shock_time != mesh->time || shock_contract.size()+shock_lineage.size()+1 >
      kContractCap) shocks::Fail("curve shock checkpoint is not synchronized");
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);  // Existing shared geometry write must be complete.
#endif
  const auto geometry = Geometry(*this);
  if (geometry.size() != shock_states.size()) {
    shocks::Fail("missing curve node histories");
  }
  Header h;
  h.rank = global_variable::my_rank;
  h.ranks = global_variable::nranks;
  h.generation = shock_generation;
  h.time = shock_time;
  h.count = geometry.size();
  h.fluid_bytes = Size(fluid_path);
  h.curve_bytes = Size(fluid_path+".curves");
  for (const auto &path : {fluid_path,fluid_path+".curves"}) {
    const int fd = open(path.c_str(),O_RDONLY);
    if (fd < 0 || fsync(fd) || close(fd)) shocks::Fail("curve checkpoint flush failed");
  }
  h.fluid_digest = shocks::FileDigest(fluid_path);
  h.curve_digest = shocks::FileDigest(fluid_path+".curves");
  h.history_digest = shock_history_digest;
  const auto contract = shock_lineage+"\n"+shock_contract;
  h.contract_bytes = contract.size();
  std::vector<NodeID> order;
  for (const auto &entry : geometry) order.push_back(entry.first);
  std::sort(order.begin(),order.end());
  std::string data(reinterpret_cast<const char *>(&h),sizeof(h));
  data += contract;
  for (NodeID id : order) {
    const auto &state = shock_states.at(id);
    CheckState(state,*geometry.at(id),shock_generation,shock_time,mesh->three_d);
    data.append(reinterpret_cast<const char *>(&state),sizeof(state));
  }
  const auto checksum = shocks::Digest(kOffset,data.data(),data.size());
  data.append(reinterpret_cast<const char *>(&checksum),sizeof(checksum));
  shocks::CommitFile(Sidecar(fluid_path),data,false);
}
void MaterialCurves::LoadShockRestart(const std::string &fluid_path) {
  if (!shock_enabled) return;
  const auto path = Sidecar(fluid_path);
  std::ifstream in(path,std::ios::binary);
  Header h;
  in.read(reinterpret_cast<char *>(&h),sizeof(h));
  const auto geometry = Geometry(*this);
  if (!in || std::memcmp(h.magic,"AKCSHCK1",8) || h.version != 1 ||
      h.endian != 0x12345678 || h.real_bytes != sizeof(Real) ||
      h.state_bytes != sizeof(ShockNode) || h.int_bytes != sizeof(int) ||
      h.rank != global_variable::my_rank || h.ranks != global_variable::nranks ||
      h.generation != static_cast<uint64_t>(mesh->ncycle) || h.time != mesh->time ||
      h.count != geometry.size() || h.count > max_nodes ||
      h.contract_bytes > kContractCap) {
    shocks::Fail("missing/incompatible curve shock sidecar");
  }
  const uint64_t fixed = sizeof(h)+h.contract_bytes+sizeof(uint64_t);
  if (h.count > (UINT64_MAX-fixed)/sizeof(ShockNode) ||
      Size(path) != fixed+h.count*sizeof(ShockNode)) {
    shocks::Fail("truncated or malformed curve shock sidecar");
  }
  in.seekg(Size(path)-sizeof(uint64_t));
  uint64_t checksum;
  in.read(reinterpret_cast<char *>(&checksum),sizeof(checksum));
  if (!in || shocks::FileDigest(path,Size(path)-sizeof(uint64_t)) != checksum) {
    shocks::Fail("curve shock sidecar checksum mismatch");
  }
  in.seekg(sizeof(h));
  std::string contract(h.contract_bytes,'\0');
  in.read(contract.data(),contract.size());
  if (!in || contract != shock_lineage+"\n"+shock_contract) {
    shocks::Fail("curve shock lineage/source/config contract changed");
  }
  if (Size(fluid_path) != h.fluid_bytes || Size(fluid_path+".curves") != h.curve_bytes ||
      shocks::FileDigest(fluid_path) != h.fluid_digest ||
      shocks::FileDigest(fluid_path+".curves") != h.curve_digest) {
    shocks::Fail("fluid/geometry/curve shock checkpoint identity mismatch");
  }
  if (PrefixDigest(*this,h.generation) != h.history_digest) {
    shocks::Fail("curve shock committed prefix missing or corrupted");
  }
  std::unordered_map<NodeID, ShockNode> restored;
  for (uint64_t n=0; n<h.count; ++n) {
    ShockNode state;
    in.read(reinterpret_cast<char *>(&state),sizeof(state));
    if (state.encounter.particle_id > UINT32_MAX) {
      shocks::Fail("curve shock node ID exceeds 32 bits");
    }
    auto node = geometry.find(SignedID(static_cast<std::uint32_t>(state.encounter.particle_id)));
    if (!in || node == geometry.end()) {
      shocks::Fail("unknown node in curve shock sidecar");
    }
    CheckState(state,*node->second,h.generation,h.time,mesh->three_d);
    if (!restored.emplace(node->first,state).second) {
      shocks::Fail("duplicate node history");
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);  // No rank applies state until every rank validates.
#endif
  shock_states.swap(restored);
  shock_generation = h.generation;
  shock_time = h.time;
  shock_history_digest = h.history_digest;
}
void MaterialCurves::WriteShockSuspended() const {
  if (!shock_enabled) return;
  std::vector<ShockEvent> events;
  for (const auto &entry : shock_states) {
    const auto &s = entry.second;
    if (s.encounter.phase == shocks::outside) continue;
    if (events.size() >= static_cast<size_t>(shock_capacity)) {
      shocks::Fail("curve suspended event overflow; increase event_capacity");
    }
    ShockEvent event;
    shocks::FillEvent(s.encounter,event.encounter,shocks::right_censored);
    event.birth = s.birth;
    event.first_level = s.first_level;
    event.last_level = s.last_level;
    for (int d=0; d<3; ++d) {
      event.first_spacing[d] = s.first_spacing[d];
      event.last_spacing[d] = s.last_spacing[d];
    }
    events.push_back(event);
  }
  shocks::CommitFile(shock_directory+"/suspended"+std::to_string(shock_generation)+
                     shocks::RankName()+".csv",
                     Records(*this,std::move(events),"SUSPENDED"),
                     true);
}
}  // namespace material_curves
