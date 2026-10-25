// Isolated persistence check: real mesh/curve structs, synthetic completed state.
// Mesh/curve constructors, sampling and Audit are stubs here; production-path
// validation separately exercises those components. Link with material_curves_io.cpp.
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "material_curves/material_curves.hpp"
#include "globals.hpp"
// No physics objects are instantiated in this serialization-only harness.

namespace global_variable {
int my_rank = 0, nranks = 1;
}
Mesh::Mesh(ParameterInput *) {}
Mesh::~Mesh() {}
MeshBlockTree::~MeshBlockTree() {}
MeshRefinement::MeshRefinement(Mesh *, ParameterInput *) :
    ncyc_since_ref("test_refinement_clock", 2) {}
MeshRefinement::~MeshRefinement() {}
GeodesicGrid::~GeodesicGrid() {}
SphericalGrid::~SphericalGrid() {}
namespace material_curves {
MaterialCurves::MaterialCurves(Mesh *m, ParameterInput *p) : mesh(m), pin_(p) {}
void MaterialCurves::Audit() const {}
void MaterialCurves::RebuildMesh() {}
std::vector<Vec> MaterialCurves::Sample(std::vector<Node *> &,
    std::vector<ShockObservation> *, const std::vector<particles::TracerField> &,
    std::vector<Real> *) {
  Fail("unexpected sampling in persistence-only check");
}
int MaterialCurves::Locate(const Vec &, Vec &) const {
  Fail("unexpected cell lookup in persistence-only check");
}
[[noreturn]] void MaterialCurves::Fail(const std::string &text) {
  std::cerr << text << '\n';
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  std::exit(1);
}
}

int main(int argc, char **argv) {
  static_assert(sizeof(material_curves::Node::id)==4);
  static_assert(sizeof(material_curves::Curve::next_sequence)==8);
  static_assert(sizeof(material_curves::Node::prev)==4);
  static_assert(sizeof(material_curves::Node::next)==4);
  static_assert(sizeof(material_curves::Node::parent_left)==4);
  static_assert(sizeof(material_curves::Node::parent_right)==4);
#if MPI_PARALLEL_ENABLED
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &global_variable::my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &global_variable::nranks);
#endif
  assert(argc == 3);
  const std::string mode = argv[1], path = argv[2];
  Kokkos::initialize();
  Mesh mesh(nullptr);
  mesh.adaptive = true;
  mesh.nmb_total = 2;
  LogicalLocation leaves[2]{{0, 0, 0, 1}, {1, 0, 0, 1}};
  mesh.lloc_eachmb = leaves;
  mesh.pmr = new MeshRefinement(&mesh, nullptr);
  mesh.pmr->ncyc_check_amr = 3;
  mesh.pmr->refinement_interval = 5;
  mesh.pmr->prolong_prims = true;
  mesh.pmr->ncyc_since_ref(0) = mode == "save" ? 2 : 0;
  mesh.pmr->ncyc_since_ref(1) = mode == "save" ? 4 : 0;
  mesh.time = 2;
  mesh.ncycle = 17;
  mesh.mesh_size = {0, 0, 0, 1, 1, 1, 0.125, 0.125, 0.125};
  mesh.mesh_indcs.nx1 = mesh.mesh_indcs.nx2 = 8;
  mesh.mesh_indcs.nx3 = 1;
  for (auto &bc : mesh.mesh_bcs) bc = BoundaryFlag::periodic;
  material_curves::MaterialCurves curves(&mesh, nullptr);
  curves.max_nodes = 100;
  curves.max_spacing = 0.1;
  curves.max_angle = 0.2;
  curves.basename = "roundtrip";
  if (mode == "save") {
    if (global_variable::my_rank == 0) {
      std::ofstream fluid(path, std::ios::binary);
      fluid << "synthetic fluid checkpoint\0completed";
    }
#if MPI_PARALLEL_ENABLED
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    for (material_curves::NodeID id : {1u,255u}) {
      if (id % global_variable::nranks !=
          static_cast<std::uint64_t>(global_variable::my_rank)) continue;
      material_curves::Curve curve;
      curve.id = id;
      curve.closed = id != 1;
      const int count = curve.closed ? 4 : 3;
      curve.period = 3;
      curve.next_sequence = count+1;
      curve.inserted = 1;
      curve.samples = 9;
      curve.migrations = 7;
      curve.unresolved_since = 1.5;
      curve.reason = 1;
      for (int i = 0; i < count; ++i) {
        material_curves::Node node;
        const int open_sequence[3]{1, 3, 2};
        const int closed_sequence[4]{1, 2, 3, 4};
        const int *sequence = curve.closed ? closed_sequence : open_sequence;
        node.id = (id << material_curves::SequenceBits(curves.curve_id_bits)) | sequence[i];
        node.a = curve.closed ? (i == 3 ? 2.5 : i) : (i == 2 ? 3 : i);
        if (i == 0) node.a = 0;
        node.x = {2*node.a, 0, 0};
        node.old = node.x;
        if (i == (curve.closed ? 3 : 1)) {
          node.parent_left = (id << material_curves::SequenceBits(curves.curve_id_bits)) | (curve.closed ? 3 : 1);
          node.parent_right = (id << material_curves::SequenceBits(curves.curve_id_bits)) | (curve.closed ? 1 : 2);
          node.birth = 1;
          node.anchor = 0;
          node.placement = 2;
        }
        curve.nodes.push_back(node);
      }
      for (int i = 0; i < count; ++i) {
        curve.nodes[i].prev = i == 0 ? (curve.closed ? curve.nodes.back().id : 0)
                                    : curve.nodes[i-1].id;
        curve.nodes[i].next = i+1 == count
            ? (curve.closed ? curve.nodes.front().id : 0) : curve.nodes[i+1].id;
      }
      curves.curves.push_back(curve);
    }
    curves.Save(path);
    curves.Snapshot(path+".snapshot");
  } else if (mode == "load") {
    curves.Load(path);
    assert(mesh.pmr->ncyc_since_ref(0) == 2 && mesh.pmr->ncyc_since_ref(1) == 4);
    std::uint64_t count = curves.curves.size();
    for (const auto &curve : curves.curves) {
      assert(curve.id % global_variable::nranks ==
             static_cast<std::uint64_t>(global_variable::my_rank));
      assert(curve.next_sequence == (curve.closed ? 5 : 4));
      assert(curve.nodes.size() == (curve.closed ? 4 : 3));
      assert(curve.reason == 1 && curve.unresolved_since == 1.5);
      assert(curve.inserted == 1 && curve.samples == 9 && curve.migrations == 7);
      const auto &born = curve.nodes[curve.closed ? 3 : 1];
      assert(born.anchor == 0 && born.placement == 2);
      assert(born.birth == 1 && born.a == (curve.closed ? 2.5 : 1));
      for (const auto &node : curve.nodes) {
        assert(node.sample_owner == -1 && node.old == node.x);
      }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
    assert(count == 2);
    curves.Snapshot(path+".restored");
  } else {
    assert(false);
  }
  delete mesh.pmr;
  mesh.pmr = nullptr;
#if MPI_PARALLEL_ENABLED
  MPI_Finalize();
#endif
}
