// Native MC AMR remap/state check with minimal constructor/count stand-ins.
// Actual RemapAfterMeshRefinement and neighbor methods execute; not a hydro oracle.
#include <cassert>
#include <cstring>
#include <iostream>
#include <string>
#include <set>
#include <vector>
#include "athena.hpp"
#include "globals.hpp"
#include "particles/particles.hpp"
#include "outputs/outputs.hpp"
#include "pgen/pgen.hpp"

static int count_calls = 0;
static std::vector<uint64_t> exit_ids;
Mesh::Mesh(ParameterInput *) {}
Mesh::~Mesh() {}
MeshBlockTree::~MeshBlockTree() {}
SphericalGrid::~SphericalGrid() {}
GeodesicGrid::~GeodesicGrid() {}
MeshBlockPack::MeshBlockPack(Mesh *pm, int lo, int hi) :
    pmesh(pm), gids(lo), gide(hi), nmb_thispack(hi-lo+1) {}
MeshBlockPack::~MeshBlockPack() {}
MeshBlock::MeshBlock(MeshBlockPack *pack, int lo, int nmb) :
    mb_gid("fixture_gid",nmb), mb_lev("fixture_level",nmb),
    mb_size("fixture_size",nmb), nghbr("fixture_neighbors",nmb,56), pmy_pack(pack) {}
void Mesh::UpdateParticleCounts() {
  ++count_calls;
  nprtcl_thisrank = pmb_pack->ppart->nprtcl_thispack;
}
namespace particles {
Particles::Particles(MeshBlockPack *pack, ParameterInput *) : pmy_pack(pack) {
  shock_enabled = true;
  particle_type = ParticleType::lagrangian_mc;
  pusher = ParticlesPusher::lagrangian_mc;
  nrdata = LMC_NREAL;
  nidata = PSEEDID+1;
  nprtcl_thispack = 4;
  Kokkos::realloc(prtcl_rdata,nrdata,nprtcl_thispack);
  Kokkos::realloc(prtcl_idata,nidata,nprtcl_thispack);
  Kokkos::realloc(shock_state,nprtcl_thispack);
}
Particles::~Particles() {}
void Particles::EmitShockExits(const std::vector<int> &indices) {
  auto states = Kokkos::create_mirror_view_and_copy(HostMemSpace(),shock_state);
  for (int p : indices) exit_ids.push_back(states(p).particle_id);
}
ParticlesBoundaryValues::ParticlesBoundaryValues(Particles *part, ParameterInput *) :
    sendlist("fixture_sendlist",1), send_count("fixture_sendcount",1), pmy_part(part) {
  nprtcl_send = nprtcl_recv = nsends = nrecvs = 0;
#if MPI_PARALLEL_ENABLED
  nsends_eachrank.resize(global_variable::nranks);
  MPI_Comm_dup(MPI_COMM_WORLD,&mpi_comm_part);
#endif
}
ParticlesBoundaryValues::~ParticlesBoundaryValues() {}
}  // namespace particles

static shocks::EncounterState MakeState(int id) {
  shocks::EncounterState s{};
  s.particle_id = (uint64_t{1}<<60)+id;
  s.sequence = (uint64_t{1}<<54)+2*id;
  s.last_generation = 500+id;
  s.phase = id%3;
  s.clean = id%2;
  s.flags = 512+id;
  s.suppressing_reentry = true;
  s.suppressed_reentries = 900+id;
  shocks::Observation *obs[] = {&s.previous,&s.first,&s.last,&s.post};
  for (int n=0; n<4; ++n) {
    auto &o = *obs[n];
    o.time = id+0.01*n;
    o.generation = 80*id+n;
    o.valid = true;
    o.inside = n%2;
    o.rho = 10*id+n;
    o.pressure = 20*id+n;
    o.temperature = 30*id+n;
    o.mach_pressure = 40*id+n;
    o.mach_temperature = 50*id+n;
    for (int d=0; d<3; ++d) {
      o.x[d] = id+0.1*n+0.01*d;
      o.image[d] = 10*id+n+d;
      o.normal[d] = id+2*n+d;
      o.upstream[d] = id+3*n+d;
      o.downstream[d] = id+4*n+d;
      s.image[d] = 200+id+d;
    }
  }
  return s;
}

static void CheckState(const shocks::EncounterState &a,
                       const shocks::EncounterState &b) {
  assert(a.particle_id == b.particle_id && a.sequence == b.sequence);
  assert(a.last_generation == b.last_generation && a.phase == b.phase);
  assert(a.clean == b.clean && a.flags == b.flags);
  assert(a.suppressing_reentry == b.suppressing_reentry);
  assert(a.suppressed_reentries == b.suppressed_reentries);
  for (int d=0; d<3; ++d) assert(a.image[d] == b.image[d]);
  const shocks::Observation *aa[] = {&a.previous,&a.first,&a.last,&a.post};
  const shocks::Observation *bb[] = {&b.previous,&b.first,&b.last,&b.post};
  for (int n=0; n<4; ++n) {
    const auto &x = *aa[n], &y = *bb[n];
    assert(x.time == y.time && x.generation == y.generation);
    assert(x.valid == y.valid && x.inside == y.inside);
    assert(x.rho == y.rho && x.pressure == y.pressure && x.temperature == y.temperature);
    assert(x.mach_pressure == y.mach_pressure && x.mach_temperature == y.mach_temperature);
    for (int d=0; d<3; ++d) {
      assert(x.x[d] == y.x[d] && x.image[d] == y.image[d]);
      assert(x.normal[d] == y.normal[d]);
      assert(x.upstream[d] == y.upstream[d] && x.downstream[d] == y.downstream[d]);
    }
  }
}

static void CheckNeighbors() {
  Mesh mesh(nullptr);
  mesh.multi_d = true;
  mesh.three_d = false;
  mesh.mesh_size.x1min = mesh.mesh_size.x2min = mesh.mesh_size.x3min = 0;
  mesh.mesh_size.x1max = mesh.mesh_size.x2max = mesh.mesh_size.x3max = 1;
  for (int d=0; d<6; ++d) mesh.mesh_bcs[d] = BoundaryFlag::periodic;
  MeshBlockPack pack(&mesh,0,3);
  MeshBlock block(&pack,0,4);
  pack.pmb = &block;
  mesh.pmb_pack = &pack;
  particles::Particles part(&pack,nullptr);
  pack.ppart = &part;
  particles::ParticlesBoundaryValues bvals(&part,nullptr);
  auto reals = Kokkos::create_mirror(part.prtcl_rdata);
  auto ints = Kokkos::create_mirror(part.prtcl_idata);
  auto states = Kokkos::create_mirror(part.shock_state);
  for (int pass=0; pass<2; ++pass) {
    auto &size = block.mb_size.h_view(0);
    size = mesh.mesh_size;
    size.x1max = size.x2max = pass == 0 ? 0.5 : 0.25;
    block.mb_lev.h_view(0) = pass == 0 ? 0 : 1;
    for (int n=0; n<56; ++n) {
      block.nghbr.h_view(0,n).gid = -1;
      block.nghbr.h_view(0,n).lev = -1;
    }
    for (int n=4; n<6; ++n) {
      auto &neighbor = block.nghbr.h_view(0,n);
      neighbor.gid = pass == 0 ? n-3 : (n == 4 ? -1 : 2);
      neighbor.rank = global_variable::my_rank;
      neighbor.lev = pass == 0 ? 1 : 0;
    }
    for (int p=0; p<4; ++p) {
      reals(IPX,p) = pass == 0 ? 0.5625 : 0.3125;
      reals(IPY,p) = pass == 0 && p%2 ? 0.375 : 0.125;
      reals(IPZ,p) = 0.5;
      ints(PGID,p) = 0;
      states(p) = MakeState(p+1);
    }
    block.mb_size.modify<HostMemSpace>(); block.mb_size.sync<DevExeSpace>();
    block.mb_lev.modify<HostMemSpace>(); block.mb_lev.sync<DevExeSpace>();
    block.nghbr.modify<HostMemSpace>(); block.nghbr.sync<DevExeSpace>();
    Kokkos::deep_copy(part.prtcl_rdata,reals);
    Kokkos::deep_copy(part.prtcl_idata,ints);
    Kokkos::deep_copy(part.shock_state,states);
    assert(bvals.SetNewPrtclGID() == TaskStatus::complete);
    auto after_int = Kokkos::create_mirror_view_and_copy(HostMemSpace(),part.prtcl_idata);
    auto after_state = Kokkos::create_mirror_view_and_copy(HostMemSpace(),part.shock_state);
    for (int p=0; p<4; ++p) {
      assert(after_int(PGID,p) == (pass == 0 ? (p%2 ? 2 : 1) : 2));
      CheckState(after_state(p),MakeState(p+1));
    }
  }
  std::cout << "PASS native AMR fine-subface and coarse-fallback selection\n";
}

static void CheckRefineOccupancy() {
  Mesh mesh(nullptr);
  mesh.multi_d = true;
  mesh.three_d = false;
  mesh.nmb_packs_thisrank = 1;
  mesh.root_level = 0;
  mesh.nmb_rootx1 = mesh.nmb_rootx2 = mesh.nmb_rootx3 = 1;
  mesh.mesh_size.x1min = mesh.mesh_size.x2min = mesh.mesh_size.x3min = 0;
  mesh.mesh_size.x1max = mesh.mesh_size.x2max = mesh.mesh_size.x3max = 1;
  mesh.mb_indcs.nx1 = mesh.mb_indcs.nx2 = 4;
  mesh.mb_indcs.nx3 = 1;
  for (int f=0; f<6; ++f) mesh.mesh_bcs[f] = BoundaryFlag::periodic;
  LogicalLocation locations[4] = {};
  int ranks[4] = {};
  mesh.lloc_eachmb = locations;
  mesh.rank_eachmb = ranks;
  mesh.nmb_total = 4;
  for (int g=0; g<4; ++g) {
    locations[g].level = 1;
    locations[g].lx1 = g%2;
    locations[g].lx2 = g/2;
  }
  MeshBlockPack pack(&mesh,0,3);
  mesh.pmb_pack = &pack;
  particles::Particles part(&pack,nullptr);
  pack.ppart = &part;
  const int count = 4096;
  part.nprtcl_thispack = count;
  Kokkos::realloc(part.prtcl_rdata,part.nrdata,count);
  Kokkos::realloc(part.prtcl_idata,part.nidata,count);
  Kokkos::realloc(part.shock_state,count);
  auto reals = Kokkos::create_mirror(part.prtcl_rdata);
  auto ints = Kokkos::create_mirror(part.prtcl_idata);
  auto states = Kokkos::create_mirror(part.shock_state);
  for (int p=0; p<count; ++p) {
    reals(IPX,p) = reals(IPY,p) = 0.125;  // a coarse dx=1/4 cell center
    reals(IPZ,p) = 0.5;
    ints(PGID,p) = 0;
    ints(PTAG,p) = p;
    states(p) = MakeState(p+1);
  }
  Kokkos::deep_copy(part.prtcl_rdata,reals);
  Kokkos::deep_copy(part.prtcl_idata,ints);
  Kokkos::deep_copy(part.shock_state,states);
  part.RemapAfterMeshRefinement();
  auto final = Kokkos::create_mirror_view_and_copy(HostMemSpace(),part.prtcl_rdata);
  int children[4] = {};
  for (int p=0; p<part.nprtcl_thispack; ++p) {
    const int ix = static_cast<int>(final(IPX,p)/0.125)%2;
    const int iy = static_cast<int>(final(IPY,p)/0.125)%2;
    ++children[ix+2*iy];
  }
  std::cout << "AMR_OCCUPANCY {\"count\":" << count << ",\"children\":["
            << children[0] << ',' << children[1] << ',' << children[2] << ','
            << children[3] << "],\"uniform_child_mass_fraction\":0.25,"
            << "\"first_position\":[" << final(IPX,0) << ',' << final(IPY,0)
            << "],\"scope\":\"native remap only, uniform-child-mass oracle\"}\n";
}

int main(int argc, char **argv) {
#if MPI_PARALLEL_ENABLED
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&global_variable::my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&global_variable::nranks);
#else
  global_variable::my_rank = 0;
  global_variable::nranks = 1;
#endif
  Kokkos::initialize(argc,argv);
  if (argc > 1 && std::string(argv[1]) == "occupancy") {
    CheckRefineOccupancy();
  } else {
  CheckNeighbors();
  {
    const std::string mode = argc > 1 ? argv[1] : "refine_coarsen";
    const int rank = global_variable::my_rank;
    Mesh mesh(nullptr);
    mesh.multi_d = true;
    mesh.three_d = false;
    mesh.nmb_packs_thisrank = 1;
    mesh.root_level = 0;
    mesh.nmb_rootx1 = mesh.nmb_rootx2 = mesh.nmb_rootx3 = 1;
    mesh.mesh_size.x1min = mesh.mesh_size.x2min = mesh.mesh_size.x3min = 0;
    mesh.mesh_size.x1max = mesh.mesh_size.x2max = mesh.mesh_size.x3max = 1;
    mesh.mb_indcs.nx1 = mesh.mb_indcs.nx2 = 4;
    mesh.mb_indcs.nx3 = 1;
    for (int d=0; d<6; ++d) mesh.mesh_bcs[d] = BoundaryFlag::periodic;
    LogicalLocation locations[4] = {};
    int ranks[4] = {};
    mesh.lloc_eachmb = locations;
    mesh.rank_eachmb = ranks;
    mesh.nmb_total = 4;
    for (int g=0; g<4; ++g) {
      locations[g].level = 1;
      locations[g].lx1 = g%2;
      locations[g].lx2 = g/2;
      ranks[g] = global_variable::nranks > 1 ? (g+1)%global_variable::nranks : 0;
    }
    MeshBlockPack pack(&mesh,0,3);
    mesh.pmb_pack = &pack;
    particles::Particles part(&pack,nullptr);
    pack.ppart = &part;
    auto reals = Kokkos::create_mirror(part.prtcl_rdata);
    auto ints = Kokkos::create_mirror(part.prtcl_idata);
    auto states = Kokkos::create_mirror(part.shock_state);
    for (int p=0; p<4; ++p) {
      for (int i=0; i<part.nrdata; ++i) reals(i,p) = rank*10000+100*p+i;
      for (int i=0; i<part.nidata; ++i) ints(i,p) = rank*10000+200*p+i;
      reals(IPX,p) = p%2 ? 0.75 : 0.25;
      reals(IPY,p) = p/2 ? 0.75 : 0.25;
      reals(IPZ,p) = 0.5;
      ints(PGID,p) = 0;
      ints(PTAG,p) = 100*rank+p;
      states(p) = MakeState(100*rank+p+1);
    }
    if (mode == "outside") reals(IPX,0) = 1;
    Kokkos::deep_copy(part.prtcl_rdata,reals);
    Kokkos::deep_copy(part.prtcl_idata,ints);
    Kokkos::deep_copy(part.shock_state,states);
    if (mode == "extent") Kokkos::resize(part.shock_state,3);
    part.RemapAfterMeshRefinement();
    auto check = [&](int expected_level) {
      const auto state = Kokkos::create_mirror_view_and_copy(HostMemSpace(),part.shock_state);
      const auto ri = Kokkos::create_mirror_view_and_copy(HostMemSpace(),part.prtcl_rdata);
      const auto ii = Kokkos::create_mirror_view_and_copy(HostMemSpace(),part.prtcl_idata);
      assert(state.extent(0) == static_cast<size_t>(part.nprtcl_thispack));
      std::set<int> seen;
      for (int p=0; p<part.nprtcl_thispack; ++p) {
        const int tag = ii(PTAG,p), source_rank = tag/100, source_p = tag%100;
        assert(seen.insert(tag).second);
        auto expected = MakeState(100*source_rank+source_p+1);
        CheckState(state(p),expected);
        assert(ii(PLASTLEVEL,p) == expected_level && ii(PLASTMOVE,p) == 0);
        assert(ri(LMC_CREATE_TIME,p) == source_rank*10000+100*source_p+LMC_CREATE_TIME);
        assert(ri(IPX,p) >= 0 && ri(IPX,p) < 1 && ri(IPY,p) >= 0 && ri(IPY,p) < 1);
        assert(ranks[ii(PGID,p)] == rank);
      }
      for (int source_rank=0; source_rank<global_variable::nranks; ++source_rank) {
        for (int source_p=0; source_p<4; ++source_p) {
          const int destination = ranks[expected_level == 1 ? source_p : 0];
          assert(seen.count(100*source_rank+source_p) == (destination == rank ? 1u : 0u));
        }
      }
    };
    check(1);
    const int after_refine = part.nprtcl_thispack;
    // The native mesh transaction has replaced the four children by their parent.
    mesh.nmb_total = 1;
    locations[0] = LogicalLocation{};
    ranks[0] = global_variable::nranks-1;
    pack.gids = 0;
    pack.gide = 0;
    pack.nmb_thispack = 1;
    part.RemapAfterMeshRefinement();
    check(0);
    assert(part.nprtcl_thispack == (rank == global_variable::nranks-1 ?
                                   4*global_variable::nranks : 0));
    assert(count_calls == 2);
    std::cout << "PASS native remap refine/coarsen rank=" << rank
              << " after_refine=" << after_refine
              << " after_coarsen=" << part.nprtcl_thispack
              << " full_state_preserved=1 count_calls=2\n";
  }
  }
  Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
  MPI_Finalize();
#endif
}
