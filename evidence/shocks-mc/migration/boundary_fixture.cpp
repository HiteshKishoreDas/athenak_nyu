// Bounded test of production bvals_part.cpp with constructor/count/exit stand-ins.
// This does not instantiate the hydro solver or qualify inter-rank migration.
#include <cassert>
#include <cstring>
#include <iostream>
#include <string>
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
  shocks::Observation *obs[] = {&s.previous,&s.first,&s.last,&s.post,&s.representative};
  for (int n=0; n<5; ++n) {
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

int main(int argc, char **argv) {
#if MPI_PARALLEL_ENABLED
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&global_variable::my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&global_variable::nranks);
  assert(global_variable::nranks == 1);
#else
  global_variable::my_rank = 0;
  global_variable::nranks = 1;
#endif
  Kokkos::initialize(argc,argv);
  {
    const std::string mode = argc > 1 ? argv[1] : "periodic2d";
    Mesh mesh(nullptr);
    mesh.multi_d = true;
    mesh.three_d = mode == "periodic3d";
    mesh.nmb_packs_thisrank = 1;
    mesh.mesh_size.x1min = mesh.mesh_size.x2min = mesh.mesh_size.x3min = 0;
    mesh.mesh_size.x1max = mesh.mesh_size.x2max = mesh.mesh_size.x3max = 1;
    for (int d=0; d<6; ++d) mesh.mesh_bcs[d] = BoundaryFlag::periodic;
    if (mode == "outflow" || mode == "mixed" || mode == "duplicate") mesh.mesh_bcs[outer_x1] = BoundaryFlag::outflow;
    MeshBlockPack pack(&mesh,0,0);
    MeshBlock block(&pack,0,1);
    pack.pmb = &block;
    mesh.pmb_pack = &pack;
    block.mb_size.h_view(0) = mesh.mesh_size;
    block.mb_lev.h_view(0) = 0;
    for (int n=0; n<56; ++n) {
      auto &neighbor = block.nghbr.h_view(0,n);
      neighbor.gid = (mode == "outflow" || mode == "mixed" || mode == "duplicate") ? -1 : 0;
      neighbor.rank = 0;
      neighbor.lev = 0;
    }
    block.mb_size.modify<HostMemSpace>();
    block.mb_size.sync<DevExeSpace>();
    block.mb_lev.modify<HostMemSpace>();
    block.mb_lev.sync<DevExeSpace>();
    block.nghbr.modify<HostMemSpace>();
    block.nghbr.sync<DevExeSpace>();
    particles::Particles part(&pack,nullptr);
    pack.ppart = &part;
    particles::ParticlesBoundaryValues bvals(&part,nullptr);
    auto reals = Kokkos::create_mirror_view(part.prtcl_rdata);
    auto ints = Kokkos::create_mirror_view(part.prtcl_idata);
    auto states = Kokkos::create_mirror(part.shock_state);
    for (int p=0; p<4; ++p) {
      for (int i=0; i<part.nrdata; ++i) reals(i,p) = 100*p+i;
      for (int i=0; i<part.nidata; ++i) ints(i,p) = 200*p+i;
      reals(IPX,p) = 0.25;
      reals(IPY,p) = 0.5;
      reals(IPZ,p) = 0.75;
      ints(PGID,p) = 0;
      states(p) = MakeState(p+1);
    }
    reals(IPX,1) = 1;  // exact upper half-open seam or outflow exit
    if (mode == "periodic2d") reals(IPY,2) = -0.125;
    if (mode == "periodic3d") reals(IPZ,2) = 1.125;
    if (mode == "image_max") states(1).image[0] = INT64_MAX;
    if (mode == "far") reals(IPX,1) = 2;
    if (mode == "extent") Kokkos::resize(part.shock_state,3);
    if (mode == "outflow" || mode == "mixed" || mode == "duplicate") reals(IPX,3) = 1.125;
    Kokkos::deep_copy(part.prtcl_rdata,reals);
    Kokkos::deep_copy(part.prtcl_idata,ints);
    if (mode != "extent") Kokkos::deep_copy(part.shock_state,states);
    if (mode == "disabled") {
      part.shock_enabled = false;
      Kokkos::realloc(part.shock_state,0);
    }
    assert(bvals.SetNewPrtclGID() == TaskStatus::complete);
    if (mode == "duplicate") {
      bvals.nprtcl_send = 1;
      Kokkos::resize(bvals.sendlist,1);
      bvals.sendlist.h_view(0) = {1,0,0};
      bvals.sendlist.modify<HostMemSpace>();
      bvals.sendlist.sync<DevExeSpace>();
    }
#if MPI_PARALLEL_ENABLED
    if (mode == "mixed") {
      bvals.nprtcl_send = 1;
      Kokkos::resize(bvals.sendlist,1);
      bvals.sendlist.h_view(0) = {0,0,0};
      bvals.sendlist.modify<HostMemSpace>();
      bvals.sendlist.sync<DevExeSpace>();
    }
    assert(bvals.CountSendsAndRecvs() == TaskStatus::complete);
    if (mode == "byte_limit") {
      bvals.sends_allranks = {ParticleMessageData(0,0,
          INT_MAX/static_cast<int>(sizeof(shocks::EncounterState))+1)};
    }
    assert(bvals.InitPrtclRecv() == TaskStatus::complete);
    assert(bvals.PackAndSendPrtcls() == TaskStatus::complete);
#endif
    int attempts = 0;
    while (bvals.RecvAndUnpackPrtcls() == TaskStatus::incomplete) {
      assert(count_calls == 0 && ++attempts < 10000);
    }
    assert(count_calls == 1);
    auto after = Kokkos::create_mirror_view_and_copy(HostMemSpace(),part.shock_state);
    auto ri = Kokkos::create_mirror_view_and_copy(HostMemSpace(),part.prtcl_rdata);
    auto ii = Kokkos::create_mirror_view_and_copy(HostMemSpace(),part.prtcl_idata);
    const int survivors[] = {0,2};
    const int count = (mode == "outflow" || mode == "mixed" || mode == "duplicate") ? 2 : 4;
    assert(part.nprtcl_thispack == count);
    for (int p=0; p<count; ++p) {
      const int old = mode == "mixed" ? (p == 0 ? 2 : 0) :
                      (mode == "outflow" ? survivors[p] : p);
      auto expected = states(old);
      if (mode != "outflow" && mode != "mixed" && old == 1) ++expected.image[0];
      if (mode == "periodic2d" && old == 2) --expected.image[1];
      if (mode == "periodic3d" && old == 2) ++expected.image[2];
      if (mode != "disabled") {
        assert(std::memcmp(&after(p),&expected,sizeof(expected)) == 0);
      }
      for (int i=0; i<part.nidata; ++i) assert(ii(i,p) == ints(i,old));
      for (int i=0; i<part.nrdata; ++i) {
        if (i != IPX && i != IPY && i != IPZ) assert(ri(i,p) == reals(i,old));
      }
    }
    if (mode == "outflow" || mode == "mixed" || mode == "duplicate") {
      assert(exit_ids.size() == 2);
      assert(exit_ids[0] == states(1).particle_id && exit_ids[1] == states(3).particle_id);
    } else {
      assert(ri(IPX,1) == (mode == "disabled" ? 1 : 0));
      if (mode == "periodic2d") assert(ri(IPY,2) == 0.875);
      if (mode == "periodic3d") assert(ri(IPZ,2) == 0.125);
    }
    if (mode == "disabled") {
      assert(part.shock_state.extent(0) == 0);
      assert(bvals.shock_destroylist.extent(0) == 0);
      assert(bvals.shock_destroy_count.extent(0) == 0);
#if MPI_PARALLEL_ENABLED
      assert(bvals.shock_sendbuf.extent(0) == 0 && bvals.shock_recvbuf.extent(0) == 0);
      assert(bvals.shock_send_req.empty() && bvals.shock_recv_req.empty());
#endif
    }
    assert(bvals.ClearPrtclSend() == TaskStatus::complete);
    assert(bvals.ClearPrtclRecv() == TaskStatus::complete);
    std::cout << "PASS " << mode << " full-state identity and single count call\n";
  }
  Kokkos::finalize();
#if MPI_PARALLEL_ENABLED
  MPI_Finalize();
#endif
}
