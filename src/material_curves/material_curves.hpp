// Copyright (C) 2026 the AthenaK contributors
// Licensed under the 3-clause BSD License (the "LICENSE").
#ifndef MATERIAL_CURVES_MATERIAL_CURVES_HPP_
#define MATERIAL_CURVES_MATERIAL_CURVES_HPP_

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "particles/tracer_fields.hpp"
#include "shocks/detector.hpp"
#include "shocks/encounter.hpp"

namespace material_curves {
using Vec = std::array<Real, 3>;

// IDs never pass through Real. The (curve, sequence) namespace is bounded explicitly.
struct Node {
  std::uint64_t id = 0, prev = 0, next = 0;
  std::uint64_t parent_left = 0, parent_right = 0;
  Real a = 0, birth = 0;
  Vec x{}, old{};  // continuous lift; old is RK scratch, checkpoints are step-complete
  std::int32_t anchor = 1;
  std::int32_t placement = 0;  // 0 initial, 1 linear, 2 limited cubic, 3 analytic seed
  std::int32_t sample_owner = -1;  // cache only, invalid after mesh/rank changes
};

struct Curve {
  std::uint64_t id = 0, next_sequence = 1;
  bool closed = false;
  std::array<std::int64_t, 3> winding{};
  Real period = 0;  // initial arc length; labels are never renumbered
  Real unresolved_since = -1;
  std::int32_t reason = 0;  // 0 valid, 1 budget, 2 label, 3 degeneracy, 4 geometry loss
  std::uint64_t inserted = 0, migrations = 0, samples = 0;
  std::vector<Node> nodes;  // canonical order by immutable a, audited independently
};

// History is optional and owned by stable node ID, never copied with geometry.
struct ShockNode {
  shocks::EncounterState encounter;
  Real birth = 0;
  int first_level = -1, last_level = -1;
  Real first_spacing[3] = {}, last_spacing[3] = {};
};
struct ShockEvent {
  shocks::Event encounter;
  Real birth = 0;
  int first_level = -1, last_level = -1;
  Real first_spacing[3] = {}, last_spacing[3] = {};
};
struct ShockObservation {
  shocks::Observation observation;
  int level = -1;
  Real spacing[3] = {};
};

class MaterialCurves {
 public:
  MaterialCurves(Mesh *mesh, ParameterInput *pin);
  void Initialize();
  void ConfigureShocks(ParameterInput *pin);
  void ObserveShocks(bool births_only = false);
  void PrepareShockStore(bool restart);
  void WriteShockChunk(const std::vector<ShockEvent> &events);
  void WriteShockTrace(const std::vector<Node *> &nodes,
                       const std::vector<ShockObservation> &samples, bool births_only);
  void PreflightShockRestart(const std::string &fluid_path) const;
  void SaveShockRestart(const std::string &fluid_path) const;
  void LoadShockRestart(const std::string &fluid_path);
  void WriteShockSuspended() const;
  void Advance(int stage);
  void CompleteStep();
  void RebuildMesh();
  void Adapt(Real birth_time, bool initial, Real affected_time = -1);
  void Audit() const;
  void Snapshot(const std::string &path) const;
  std::string ThermoHistory(const std::vector<particles::TracerField> &fields);
  void Save(const std::string &fluid_path) const;
  void Load(const std::string &fluid_path);
  [[noreturn]] static void Fail(const std::string &message);

  // Stable persistence interface. Curves have home rank id % nranks; sampling
  // requests go to current leaf owners. No whole-curve replication across ranks.
  Mesh *mesh;
  std::vector<Curve> curves;
  std::uint64_t max_nodes = 1000000;
  Real max_spacing = 0, max_angle = 0.2, chord_error = 0;
  Real spacing_cells = 1;
  bool adaptive = true, binary_output = false;
  int output_dcycle = 0;
  std::string basename;

  bool shock_enabled = false, shock_trace = false;
  shocks::DetectorConfig shock_config;
  int shock_capacity = 4096;
  std::string shock_directory, shock_lineage, shock_contract;
  std::unordered_map<std::uint64_t, ShockNode> shock_states;
  std::uint64_t shock_generation = UINT64_MAX;
  Real shock_time = -1;
  std::uint64_t shock_history_digest = UINT64_C(14695981039346656037);

 private:
  ParameterInput *pin_;
  struct CircleSeed {Vec center, e1, e2; Real radius;};
  std::unordered_map<std::uint64_t, CircleSeed> circle_seeds_;
  std::unordered_map<LogicalLocation, int, LogicalLocationHash> leaves_;
  int finest_level_ = 0;
  Real mesh_spacing_ = 0;
  std::vector<Vec> Sample(std::vector<Node *> &nodes,
                          std::vector<ShockObservation> *observations = nullptr,
                          const std::vector<particles::TracerField> &fields = {},
                          std::vector<Real> *values = nullptr);
  int Locate(const Vec &lift, Vec &wrapped) const;
  Vec Edge(const Curve &curve, std::size_t i) const;
  void Links(Curve &curve);
  std::uint64_t NewID(Curve &curve);
};
}  // namespace material_curves
#endif  // MATERIAL_CURVES_MATERIAL_CURVES_HPP_
