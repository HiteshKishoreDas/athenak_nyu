// Copyright (C) 2026 the AthenaK contributors
// Licensed under the 3-clause BSD License (the "LICENSE").
// Passive material curves: separate identity/geometry from fluid leaf ownership.
#include "material_curves.hpp"
#include "particles/particles.hpp"
#include "particles/fluid_velocity.hpp"
#include "shocks/observation.hpp"

#include <algorithm>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <unordered_set>

#include "globals.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "srcterms/srcterms.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace material_curves {
namespace {
Real Norm(const Vec &v) {
  return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
Real Dot(const Vec &a, const Vec &b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
Vec Difference(const Vec &a, const Vec &b) {
  return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
Vec Cross(const Vec &a, const Vec &b) {
  return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
std::uint64_t Integer(const std::string &text) {
  if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos) {
    MaterialCurves::Fail("expected an unsigned decimal integer: " + text);
  }
  try {
    return std::stoull(text);
  } catch (...) {
    MaterialCurves::Fail("integer out of range: " + text);
  }
}
std::int64_t SignedInteger(const std::string &text) {
  std::size_t start = (!text.empty() && text[0] == '-') ? 1 : 0;
  if (text.size() <= start ||
      text.find_first_not_of("0123456789", start) != std::string::npos) {
    MaterialCurves::Fail("expected signed decimal winding: " + text);
  }
  try {
    return std::stoll(text);
  } catch (...) {
    MaterialCurves::Fail("winding out of range: " + text);
  }
}
}  // namespace

[[noreturn]] void MaterialCurves::Fail(const std::string &message) {
  std::cerr << "### FATAL ERROR: material curves: " << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

MaterialCurves::MaterialCurves(Mesh *pm, ParameterInput *pin) : mesh(pm), pin_(pin) {
  auto pp = mesh->pmb_pack;
  if (pp->phydro == nullptr || pp->pmhd || pp->prad || pp->pionn || pp->pz4c ||
      pp->pdyngr || pp->padm || (pp->ppart && !pp->ppart->IsFluidTracer()) ||
      pp->pcoord->is_special_relativistic || pp->pcoord->is_general_relativistic ||
      pp->pcoord->is_dynamical_relativistic) {
    Fail("requires passive Newtonian single-fluid HD; only fluid tracers may coexist");
  }
  if (pin->GetOrAddString("time", "integrator", "rk2") != "rk2" ||
      pin->GetOrAddString("time", "evolution", "dynamic") != "dynamic") {
    Fail("only dynamic HD with RK2 is supported");
  }
  if (!mesh->multi_d) Fail("only 2D and 3D Cartesian curves are supported");
  for (int d=0; d<(mesh->three_d ? 3 : 2); ++d) {
    if (mesh->mesh_bcs[2*d] != BoundaryFlag::periodic ||
        mesh->mesh_bcs[2*d+1] != BoundaryFlag::periodic) {
      Fail("active dimensions must have periodic boundaries; exits are unsupported");
    }
  }
  auto src = pp->phydro->psrc;
  if (src && (src->CoolingEnabled() || src->rel_cooling)) {
    Fail("cooling source stages are unqualified");
  }
  if (pp->phydro->psbox_u || pp->phydro->porb_u) {
    Fail("shearing/orbital advection stages are unqualified");
  }
  if (pp->pframe_tracker) {
    Fail("frame tracking is unqualified for material curve coordinates");
  }
#if SINGLE_PRECISION_ENABLED
  Fail("single precision is unqualified");
#endif
#if defined(KOKKOS_ENABLE_SYCL)
  Fail("SYCL material curves are unsupported");
#endif
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  if (!pin->GetOrAddBoolean("material_curves", "allow_unqualified_accelerator", false)) {
    Fail("CUDA/HIP requires allow_unqualified_accelerator=true until device gates pass");
  }
#endif
  if (global_variable::my_rank==0) {
    std::cout << "Material curve sampling/push execution space: "
              << DevExeSpace::name() << "; topology/adaptation: host" << std::endl;
  }
  const auto bits=Integer(pin->GetOrAddString("material_curves", "curve_id_bits", "8"));
  if (bits>30) Fail("curve_id_bits must lie in [0,30]");
  curve_id_bits=static_cast<unsigned>(bits);
  max_nodes = Integer(pin->GetOrAddString("material_curves", "max_nodes", "1000000"));
  if (max_nodes == 0 ||
      max_nodes > static_cast<std::uint64_t>(std::numeric_limits<int>::max()/3)) {
    Fail("max_nodes exceeds bounded MPI/device batch count");
  }
  spacing_cells = pin->GetOrAddReal("material_curves", "spacing_cells", 1.0);
  max_spacing = pin->GetOrAddReal("material_curves", "max_spacing", 0.0);
  max_angle = pin->GetOrAddReal("material_curves", "max_angle", 0.2);
  chord_error = pin->GetOrAddReal("material_curves", "chord_error", 0.0);
  adaptive = pin->GetOrAddBoolean("material_curves", "adaptive", true);
  output_dcycle = pin->GetOrAddInteger("material_curves", "output_dcycle", 0);
  basename = pin->GetString("job", "basename");
  if (!(spacing_cells > 0) || !std::isfinite(spacing_cells) ||
      max_spacing < 0 || !std::isfinite(max_spacing) ||
      !(max_angle > 0 && max_angle < M_PI) || !std::isfinite(max_angle) ||
      chord_error < 0 || !std::isfinite(chord_error) || output_dcycle < 0) {
    Fail("invalid spacing, angle, error tolerance, or output cadence");
  }
  RebuildMesh();
  ConfigureShocks(pin);
}

NodeID MaterialCurves::NewID(Curve &c) {
  if (c.id == 0 || c.id > MaxCurveID(curve_id_bits) || c.next_sequence == 0 ||
      c.next_sequence > MaxSequence(curve_id_bits)) {
    Fail("curve/sequence ID namespace exhausted");
  }
  const NodeID prefix=curve_id_bits ? c.id << SequenceBits(curve_id_bits) : 0;
  return prefix | static_cast<NodeID>(c.next_sequence++);
}

void MaterialCurves::Links(Curve &c) {
  for (std::size_t i=0; i<c.nodes.size(); ++i) {
    c.nodes[i].prev = i ? c.nodes[i-1].id : (c.closed ? c.nodes.back().id : 0);
    c.nodes[i].next = i+1<c.nodes.size() ? c.nodes[i+1].id :
                     (c.closed ? c.nodes.front().id : 0);
  }
}

void MaterialCurves::Initialize() {
  int nc = pin_->GetOrAddInteger("material_curves", "ncurves", 1);
  if (nc < 0 || static_cast<std::uint64_t>(nc) > MaxCurveID(curve_id_bits) || nc>100000) {
    Fail("ncurves exceeds curve_id_bits namespace or 100000-curve IO cap");
  }
  std::uint64_t total = 0;
  std::unordered_set<NodeID> ids;
  for (int k=1; k<=nc; ++k) {
    const std::string block = "curve" + std::to_string(k);
    Curve c;
    const auto requested_id = Integer(
        pin_->GetOrAddString(block, "id", std::to_string(k)));
    if (requested_id == 0 || requested_id > MaxCurveID(curve_id_bits)) {
      Fail("curve IDs must be unique within curve_id_bits namespace");
    }
    c.id = static_cast<NodeID>(requested_id);
    if (c.id == 0 || c.id > MaxCurveID(curve_id_bits) || !ids.insert(c.id).second) {
      Fail("curve IDs must be unique within curve_id_bits namespace");
    }
    auto type = pin_->GetOrAddString(block, "type", "line");
    c.closed = pin_->GetOrAddBoolean(block, "closed", type == "circle");
    for (int d=0; d<3; ++d) {
      c.winding[d] = SignedInteger(pin_->GetOrAddString(block,
                                  "winding"+std::to_string(d+1), "0"));
    }
    if ((!c.closed && (c.winding[0] || c.winding[1] || c.winding[2])) ||
        (!mesh->three_d && c.winding[2])) Fail("invalid winding for this curve");
    std::uint64_t requested=0;
    if (type=="polyline") {
      const auto text=pin_->GetString(block,"points");
      if (text.empty() || text.back()==';') Fail("empty/trailing polyline point");
      requested=1+std::count(text.begin(),text.end(),';');
    } else if (type=="line" || type=="circle") {
      int count=pin_->GetOrAddInteger(block,"count",type=="circle" ? 32 : 16);
      if (count<(c.closed ? 3 : 2)) Fail("invalid seed sample count");
      requested=count;
    } else {
      Fail("seed type must be line, circle, or polyline");
    }
    if (requested>max_nodes-total) Fail("global seed count exceeds max_nodes");
    total+=requested;
    // Only the home rank constructs samples; invalid input aborts all ranks.
    if (c.id%global_variable::nranks!=
        static_cast<unsigned>(global_variable::my_rank)) continue;
    std::vector<Vec> positions;
    if (type == "polyline") {
      if (pin_->GetString(block, "coordinates") != "lifted") {
        Fail("polyline requires coordinates=lifted (unambiguous periodic images)");
      }
      std::istringstream stream(pin_->GetString(block, "points"));
      std::string item;
      while (std::getline(stream, item, ';')) {
        Vec x;
        std::istringstream point(item);
        std::string extra;
        if (!(point >> x[0] >> x[1] >> x[2]) || (point >> extra)) {
          Fail("polyline points must be semicolon-separated x y z triples");
        }
        positions.push_back(x);
        if (positions.size() > max_nodes) Fail("seed exceeds max_nodes");
      }
    } else if (type == "line" || type == "circle") {
      int count = static_cast<int>(requested);
      Vec origin;
      for (int d=0; d<3; ++d) {
        origin[d] = pin_->GetOrAddReal(block, "x"+std::to_string(d+1), 0.0);
      }
      if (type == "line") {
        Vec end;
        for (int d=0; d<3; ++d) {
          end[d] = pin_->GetOrAddReal(block, "end"+std::to_string(d+1), origin[d]);
        }
        if (c.closed) {
          Vec closure = origin;
          Real length[3] = {mesh->mesh_size.x1max-mesh->mesh_size.x1min,
                            mesh->mesh_size.x2max-mesh->mesh_size.x2min,
                            mesh->mesh_size.x3max-mesh->mesh_size.x3min};
          for (int d=0; d<3; ++d) closure[d] += c.winding[d]*length[d];
          if (Norm(Difference(end, closure)) >
              32*std::numeric_limits<Real>::epsilon()*Norm(Difference(end, origin))) {
            Fail("closed line end must equal origin plus winding times domain");
          }
        }
        for (int i=0; i<count; ++i) {
          Real a = static_cast<Real>(i)/(c.closed ? count : count-1);
          Vec x;
          for (int d=0; d<3; ++d) x[d] = origin[d]+a*(end[d]-origin[d]);
          positions.push_back(x);
        }
      } else {
        if (!c.closed || c.winding[0] || c.winding[1] || c.winding[2]) {
          Fail("circle seed must be a closed zero-winding loop");
        }
        Real radius = pin_->GetReal(block, "radius");
        Vec normal;
        for (int d=0; d<3; ++d) {
          normal[d] = pin_->GetOrAddReal(block, "normal"+std::to_string(d+1),
                                        d == 2 ? 1.0 : 0.0);
        }
        Real n = Norm(normal);
        if (!(radius > 0) || !std::isfinite(radius) || !(n > 0) || !std::isfinite(n)) {
          Fail("circle radius and normal must be finite and nonzero");
        }
        for (auto &x : normal) x /= n;
        int axis = std::min_element(normal.begin(), normal.end(),
                   [](Real a, Real b) {return std::abs(a)<std::abs(b);})-normal.begin();
        Vec basis{};
        basis[axis] = 1;
        Vec e1 = Cross(normal, basis);
        Real norm = Norm(e1);
        for (auto &x : e1) x /= norm;
        Vec e2 = Cross(normal, e1);
        circle_seeds_[c.id] = {origin, e1, e2, radius};
        for (int i=0; i<count; ++i) {
          Real a = 2*M_PI*i/count;
          Vec x;
          for (int d=0; d<3; ++d) {
            x[d] = origin[d]+radius*(std::cos(a)*e1[d]+std::sin(a)*e2[d]);
          }
          positions.push_back(x);
        }
        c.period = 2*M_PI*radius;
      }
    } else {
      Fail("seed type must be line, circle, or polyline");
    }
    if (positions.size() < (c.closed ? 3u : 2u)) Fail("too few seed samples");
    Real a = 0;
    for (std::size_t i=0; i<positions.size(); ++i) {
      Vec wrapped;
      Locate(positions[i], wrapped);  // finite/lift validation, including at seams
      if (!mesh->three_d && positions[i][2] != positions[0][2]) {
        Fail("2D curves must be planar with constant x3");
      }
      if (i) {
        Real ds = Norm(Difference(positions[i], positions[i-1]));
        if (!(ds > 0) || !std::isfinite(ds)) Fail("degenerate seed edge");
        a += ds;
      }
      Node node;
      node.id = NewID(c);
      node.x = node.old = positions[i];
      node.a = type == "circle" ? c.period*i/positions.size() : a;
      node.birth = mesh->time;
      c.nodes.push_back(node);
    }
    if (type != "circle") {
      if (c.closed) {
        Real closure = Norm(Edge(c, c.nodes.size()-1));
        if (!(closure > 0) || !std::isfinite(closure)) Fail("degenerate closure edge");
        a += closure;
      }
      c.period = a;
    }
    Links(c);
    curves.push_back(std::move(c));
  }
  Audit();
  if (adaptive) Adapt(mesh->time, true);
  circle_seeds_.clear();
}

void MaterialCurves::RebuildMesh() {
  leaves_.clear();
  finest_level_ = mesh->root_level;
  for (int g=0; g<mesh->nmb_total; ++g) {
    leaves_.emplace(mesh->lloc_eachmb[g], g);
    finest_level_ = std::max(finest_level_, static_cast<int>(mesh->lloc_eachmb[g].level));
  }
  Real h = std::min(mesh->mesh_size.dx1, mesh->mesh_size.dx2);
  if (mesh->three_d) h = std::min(h, mesh->mesh_size.dx3);
  // Conservative: a long edge cannot miss an interior fine patch.
  mesh_spacing_ = spacing_cells*std::ldexp(h, mesh->root_level-finest_level_);
}

int MaterialCurves::Locate(const Vec &x, Vec &wrapped) const {
  const Real low[3] = {mesh->mesh_size.x1min, mesh->mesh_size.x2min,
                       mesh->mesh_size.x3min};
  const Real length[3] = {mesh->mesh_size.x1max-low[0], mesh->mesh_size.x2max-low[1],
                          mesh->mesh_size.x3max-low[2]};
  Real fraction[3]{};
  for (int d=0; d<3; ++d) {
    if (!std::isfinite(x[d])) Fail("nonfinite lifted position");
    long double f = (static_cast<long double>(x[d])-low[d])/length[d];
    if (std::abs(f) > (1LL<<40)) Fail("lift exceeds reliable floating image range");
    fraction[d] = static_cast<Real>(f-std::floor(f));
    wrapped[d] = low[d] + fraction[d]*length[d];
    if (fraction[d]>=1 || wrapped[d]>=low[d]+length[d]) {
      fraction[d]=0;
      wrapped[d]=low[d];
    }
  }
  for (int level=finest_level_; level>=mesh->root_level; --level) {
    int scale = level-mesh->root_level;
    if (scale > 30 ||
        std::ldexp(mesh->nmb_rootx1,scale)>INT32_MAX ||
        std::ldexp(mesh->nmb_rootx2,scale)>INT32_MAX ||
        std::ldexp(mesh->nmb_rootx3,scale)>INT32_MAX) {
      Fail("logical mesh index exceeds supported integer range");
    }
    LogicalLocation loc;
    loc.level = level;
    loc.lx1 = std::floor(fraction[0]*std::ldexp(mesh->nmb_rootx1, scale));
    loc.lx2 = std::floor(fraction[1]*std::ldexp(mesh->nmb_rootx2, scale));
    loc.lx3 = mesh->three_d ?
              std::floor(fraction[2]*std::ldexp(mesh->nmb_rootx3, scale)) : 0;
    auto leaf = leaves_.find(loc);
    if (leaf != leaves_.end()) return leaf->second;
  }
  Fail("position has no unique active leaf");
}

Vec MaterialCurves::Edge(const Curve &c, std::size_t i) const {
  std::size_t j = (i+1)%c.nodes.size();
  Vec v = Difference(c.nodes[j].x, c.nodes[i].x);
  if (j == 0) {
    Real length[3] = {mesh->mesh_size.x1max-mesh->mesh_size.x1min,
                      mesh->mesh_size.x2max-mesh->mesh_size.x2min,
                      mesh->mesh_size.x3max-mesh->mesh_size.x3min};
    for (int d=0; d<3; ++d) v[d] += c.winding[d]*length[d];
  }
  return v;
}

std::vector<Vec> MaterialCurves::Sample(std::vector<Node *> &nodes,
                                        std::vector<ShockObservation> *observations,
                                        const std::vector<particles::TracerField> &fields,
                                        std::vector<Real> *values) {
  const int ranks = global_variable::nranks;
  std::vector<std::vector<int>> buckets(ranks);
  std::vector<Vec> wrapped(nodes.size());
  std::vector<int> gids(nodes.size());
  for (std::size_t p=0; p<nodes.size(); ++p) {
    gids[p] = Locate(nodes[p]->x, wrapped[p]);
    int rank = mesh->rank_eachmb[gids[p]];
    buckets[rank].push_back(p);
  }
  std::vector<int> sc(ranks), rc(ranks), sd(ranks), rd(ranks);
  for (int r=0; r<ranks; ++r) sc[r] = buckets[r].size();
#if MPI_PARALLEL_ENABLED
  MPI_Alltoall(sc.data(), 1, MPI_INT, rc.data(), 1, MPI_INT, MPI_COMM_WORLD);
#else
  rc = sc;
#endif
  std::uint64_t incoming = 0;
  for (int r=0; r<ranks; ++r) incoming += rc[r];
  if (incoming > max_nodes || incoming > static_cast<std::uint64_t>(INT_MAX/3)) {
    Fail("sampling rank batch exceeds bounded allocation/count limit");
  }
  for (int r=1; r<ranks; ++r) {
    sd[r] = sd[r-1]+sc[r-1];
    rd[r] = rd[r-1]+rc[r-1];
  }
  int ns = nodes.size(), nr = incoming;
  std::vector<Real> send(3*ns), recv(3*nr), velocity(3*nr), returned(3*ns);
  std::vector<int> sendg(ns), recvg(nr), order(ns);
  std::vector<NodeID> sendid(ns), recvid(nr), backid(ns);
  for (int r=0; r<ranks; ++r) {
    for (int j=0; j<sc[r]; ++j) {
      int p = buckets[r][j], q = sd[r]+j;
      for (int d=0; d<3; ++d) send[3*q+d] = wrapped[p][d];
      sendg[q] = gids[p];
      sendid[q] = nodes[p]->id;
      order[q] = p;
    }
  }
  std::vector<int> sc3(ranks), rc3(ranks), sd3(ranks), rd3(ranks);
  for (int r=0; r<ranks; ++r) {
    sc3[r]=3*sc[r]; rc3[r]=3*rc[r]; sd3[r]=3*sd[r]; rd3[r]=3*rd[r];
  }
#if MPI_PARALLEL_ENABLED
  MPI_Alltoallv(send.data(), sc3.data(), sd3.data(), MPI_DOUBLE,
                recv.data(), rc3.data(), rd3.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(sendg.data(), sc.data(), sd.data(), MPI_INT,
                recvg.data(), rc.data(), rd.data(), MPI_INT, MPI_COMM_WORLD);
  MPI_Alltoallv(sendid.data(), sc.data(), sd.data(), MPI_UINT32_T,
                recvid.data(), rc.data(), rd.data(), MPI_UINT32_T, MPI_COMM_WORLD);
#else
  recv=send; recvg=sendg; recvid=sendid;
#endif
  DvceArray2D<Real> q("curve_sample_positions", nr, 3);
  DvceArray2D<Real> v("curve_sample_velocity", nr, 3);
  DvceArray1D<int> block("curve_sample_block", nr), status("curve_sample_status", nr);
  auto qh=Kokkos::create_mirror_view(q);
  auto bh=Kokkos::create_mirror_view(block);
  for (int p=0; p<nr; ++p) {
    bh(p)=recvg[p]-mesh->pmb_pack->gids;
    if (bh(p)<0 || bh(p)>=mesh->pmb_pack->nmb_thispack) Fail("sample leaf/rank mismatch");
    for (int d=0; d<3; ++d) qh(p,d)=recv[3*p+d];
  }
  Kokkos::deep_copy(q,qh);
  Kokkos::deep_copy(block,bh);
  auto w=mesh->pmb_pack->phydro->w0;
  auto sizes=mesh->pmb_pack->pmb->mb_size.d_view;
  auto ind=mesh->mb_indcs;
  bool three=mesh->three_d;
  if (values) {
    const std::size_t nf=fields.size();
    if (!nf || nf>INT_MAX || ns>INT_MAX/nf || nr>INT_MAX/nf) {
      Fail("thermo field batch exceeds MPI count limit");
    }
    // Reuse the particle diagnostic evaluator on the leaf that owns each request.
    // ponytail: host mirror matches particle output; device evaluation if profiling warrants.
    HostArray5D<Real> hw=Kokkos::create_mirror_view_and_copy(HostMemSpace(),w);
    HostArray5D<Real> hb;
    auto hs=Kokkos::create_mirror_view_and_copy(HostMemSpace(),sizes);
    std::vector<Real> sampled(nf*nr), returned_values(nf*ns);
    auto hydro=mesh->pmb_pack->phydro;
    Real gamma=hydro->peos->eos_data.gamma;
    for (int p=0; p<nr; ++p) {
      int m=bh(p);
      auto size=hs(m);
      int i=std::clamp(static_cast<int>((qh(p,0)-size.x1min)/size.dx1)+ind.is,
                       ind.is,ind.ie);
      int j=std::clamp(static_cast<int>((qh(p,1)-size.x2min)/size.dx2)+ind.js,
                       ind.js,ind.je);
      int k=three ? std::clamp(static_cast<int>((qh(p,2)-size.x3min)/size.dx3)+ind.ks,
                               ind.ks,ind.ke) : ind.ks;
      for (std::size_t n=0; n<nf; ++n) {
        Real value=particles::EvaluateTracerFieldHost(fields[n],hw,hb,false,gamma,
                                                       hydro->nhydro,m,k,j,i);
        if (!std::isfinite(value)) Fail("nonfinite curve thermo field: "+fields[n].name);
        sampled[nf*p+n]=value;
      }
    }
#if MPI_PARALLEL_ENABLED
    for (int r=0; r<ranks; ++r) {
      sc3[r]=nf*sc[r]; rc3[r]=nf*rc[r]; sd3[r]=nf*sd[r]; rd3[r]=nf*rd[r];
    }
    MPI_Alltoallv(sampled.data(),rc3.data(),rd3.data(),MPI_ATHENA_REAL,
                  returned_values.data(),sc3.data(),sd3.data(),MPI_ATHENA_REAL,
                  MPI_COMM_WORLD);
    MPI_Alltoallv(recvid.data(),rc.data(),rd.data(),MPI_UINT32_T,
                  backid.data(),sc.data(),sd.data(),MPI_UINT32_T,MPI_COMM_WORLD);
#else
    returned_values=std::move(sampled); backid=recvid;
#endif
    values->resize(nf*ns);
    for (int q=0; q<ns; ++q) {
      int p=order[q];
      if (backid[q]!=nodes[p]->id) Fail("thermo response identity mismatch");
      for (std::size_t n=0; n<nf; ++n) (*values)[nf*p+n]=returned_values[nf*q+n];
    }
    return {};
  }
  if (observations) {
    if (ns>INT_MAX/static_cast<int>(sizeof(ShockObservation)) ||
        nr>INT_MAX/static_cast<int>(sizeof(ShockObservation))) {
      Fail("shock observation batch exceeds MPI byte-count limit");
    }
    DvceArray1D<ShockObservation> sample("curve_shock_observations",nr);
    auto config=shock_config;
    auto levels=mesh->pmb_pack->pmb->mb_lev.d_view;
    Real time=mesh->time;
    std::uint64_t generation=mesh->ncycle;
    par_for("material_curve_shock", DevExeSpace(), 0, nr-1,
    KOKKOS_LAMBDA(int p) {
      int m=block(p);
      auto size=sizes(m);
      int i=static_cast<int>(floor((q(p,0)-size.x1min)/size.dx1))+ind.is;
      int j=static_cast<int>(floor((q(p,1)-size.x2min)/size.dx2))+ind.js;
      int k=three ? static_cast<int>(floor((q(p,2)-size.x3min)/size.dx3))+ind.ks : 0;
      Real position[3]={q(p,0),q(p,1),q(p,2)};
      ShockObservation out;
      out.level=levels(m);
      out.spacing[0]=size.dx1; out.spacing[1]=size.dx2; out.spacing[2]=size.dx3;
      auto signal=shocks::DetectCell(w,m,k,j,i,size.dx1,size.dx2,size.dx3,three,config);
      out.observation=shocks::ObserveSignal(signal,position,time,generation);
      sample(p)=out;
    });
    auto host=Kokkos::create_mirror_view_and_copy(HostMemSpace(),sample);
    std::vector<ShockObservation> returned_samples(ns);
#if MPI_PARALLEL_ENABLED
    for (int r=0; r<ranks; ++r) {
      sc3[r]=sc[r]*sizeof(ShockObservation); rc3[r]=rc[r]*sizeof(ShockObservation);
      sd3[r]=sd[r]*sizeof(ShockObservation); rd3[r]=rd[r]*sizeof(ShockObservation);
    }
    MPI_Alltoallv(host.data(),rc3.data(),rd3.data(),MPI_BYTE,
                  returned_samples.data(),sc3.data(),sd3.data(),MPI_BYTE,MPI_COMM_WORLD);
    MPI_Alltoallv(recvid.data(),rc.data(),rd.data(),MPI_UINT32_T,
                  backid.data(),sc.data(),sd.data(),MPI_UINT32_T,MPI_COMM_WORLD);
#else
    for (int p=0; p<nr; ++p) returned_samples[p]=host(p);
    backid=recvid;
#endif
    observations->resize(ns);
    for (int q=0; q<ns; ++q) {
      int p=order[q];
      if (backid[q]!=nodes[p]->id) Fail("shock response identity mismatch");
      (*observations)[p]=returned_samples[q];
    }
    return {};
  }
  par_for("material_curve_velocity", DevExeSpace(), 0, nr-1,
  KOKKOS_LAMBDA(int p) {
    int m=block(p);
    Real position[3]={q(p,0),q(p,1),q(p,2)}, velocity[3];
    status(p)=particles::FluidVelocity(w,m,position,sizes(m),ind,three,velocity);
    for (int d=0; d<3; ++d) v(p,d)=velocity[d];
  });
  auto vh=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),v);
  auto sh=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),status);
  for (int p=0; p<nr; ++p) {
    if (sh(p)) Fail(sh(p)==1 ? "interpolation stencil out of bounds" :
                             "invalid density/velocity in interpolation stencil");
    for (int d=0; d<3; ++d) velocity[3*p+d]=vh(p,d);
  }
#if MPI_PARALLEL_ENABLED
  MPI_Alltoallv(velocity.data(), rc3.data(), rd3.data(), MPI_DOUBLE,
                returned.data(), sc3.data(), sd3.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Alltoallv(recvid.data(), rc.data(), rd.data(), MPI_UINT32_T,
                backid.data(), sc.data(), sd.data(), MPI_UINT32_T, MPI_COMM_WORLD);
#else
  returned=velocity; backid=recvid;
#endif
  std::vector<Vec> result(ns);
  for (int q=0; q<ns; ++q) {
    int p=order[q];
    if (backid[q] != nodes[p]->id) Fail("sample response identity mismatch");
    for (int d=0; d<3; ++d) result[p][d]=returned[3*q+d];
  }
  for (auto &c : curves) {
    for (auto &n : c.nodes) {
      Vec x;
      int owner=mesh->rank_eachmb[Locate(n.x,x)];
      if (n.sample_owner>=0 && owner!=n.sample_owner) ++c.migrations;
      n.sample_owner=owner;
      ++c.samples;
    }
  }
  return result;
}

void MaterialCurves::Advance(int stage) {
  if (stage != 1 && stage != 2) Fail("invalid RK2 curve stage");
  if (stage==1) {
    RebuildMesh();
    if (adaptive) Adapt(mesh->time, false);
    if (shock_enabled) ObserveShocks(true);
  }
  std::vector<Node *> nodes;
  for (auto &c : curves) for (auto &n : c.nodes) nodes.push_back(&n);
  auto velocity=Sample(nodes);
  int count=nodes.size(), dimensions=mesh->three_d ? 3 : 2;
  DvceArray2D<Real> x("curve_push_lift", count, 3);
  DvceArray2D<Real> old("curve_push_old", count, 3);
  DvceArray2D<Real> v("curve_push_velocity", count, 3);
  auto xh=Kokkos::create_mirror_view(x);
  auto oh=Kokkos::create_mirror_view(old);
  auto vh=Kokkos::create_mirror_view(v);
  for (int i=0; i<count; ++i) {
    for (int d=0; d<3; ++d) {
      xh(i,d)=nodes[i]->x[d];
      oh(i,d)=nodes[i]->old[d];
      vh(i,d)=velocity[i][d];
    }
  }
  Kokkos::deep_copy(x,xh);
  Kokkos::deep_copy(old,oh);
  Kokkos::deep_copy(v,vh);
  Real dt=mesh->dt;
  par_for("material_curve_push", DevExeSpace(), 0, count-1,
  KOKKOS_LAMBDA(int i) {
    for (int d=0; d<3; ++d) {
      if (stage==1) old(i,d)=x(i,d);
      if (d<dimensions) {
        Real euler=x(i,d)+dt*v(i,d);
        x(i,d)=stage==1 ? euler : 0.5*old(i,d)+0.5*euler;
      }
    }
  });
  Kokkos::deep_copy(xh,x);
  if (stage==1) Kokkos::deep_copy(oh,old);
  for (int i=0; i<count; ++i) {
    for (int d=0; d<3; ++d) {
      nodes[i]->x[d]=xh(i,d);
      if (stage==1) nodes[i]->old[d]=oh(i,d);
    }
  }
}

void MaterialCurves::Audit() const {
  std::uint64_t count=0;
  std::unordered_set<NodeID> curve_ids;
  for (const auto &c : curves) {
    if (c.id==0 || c.id>MaxCurveID(curve_id_bits) || !curve_ids.insert(c.id).second ||
        c.id%global_variable::nranks!=static_cast<unsigned>(global_variable::my_rank) ||
        c.nodes.size()<(c.closed ? 3u : 2u) || !(c.period>0) ||
        !std::isfinite(c.period) || c.next_sequence==0 ||
        c.next_sequence>MaxSequence(curve_id_bits)+1) {
      Fail("invalid curve identity, owner, period, count, or allocator");
    }
    if ((c.reason==0)!=(c.unresolved_since==-1) || c.reason<0 || c.reason>4 ||
        !std::isfinite(c.unresolved_since)) Fail("invalid persistent quality state");
    if ((!c.closed && (c.winding[0] || c.winding[1] || c.winding[2])) ||
        (!mesh->three_d && c.winding[2]) || c.nodes.front().a!=0 ||
        (!c.closed && c.nodes.back().a!=c.period) ||
        (c.reason && (c.unresolved_since<0 || c.unresolved_since>mesh->time))) {
      Fail("invalid winding, material endpoints, or quality time");
    }
    for (auto w : c.winding) {
      if (std::abs(static_cast<long double>(w))>(1LL<<40)) {
        Fail("winding exceeds reliable lifted geometry range");
      }
    }
    std::unordered_map<NodeID, const Node *> by_id;
    for (std::size_t i=0; i<c.nodes.size(); ++i) {
      const auto &n=c.nodes[i];
      if (CurveID(n.id,curve_id_bits)!=c.id || (n.id&MaxSequence(curve_id_bits))==0 ||
          (n.id&MaxSequence(curve_id_bits))>=c.next_sequence ||
          !by_id.emplace(n.id, &n).second ||
          !std::isfinite(n.a) || n.a<0 || n.a>c.period ||
          (c.closed && n.a>=c.period) || (i && n.a<=c.nodes[i-1].a) ||
          !std::isfinite(n.birth) || n.birth<0 || n.birth>mesh->time ||
          (n.anchor!=0 && n.anchor!=1) || n.placement<0 || n.placement>3) {
        Fail("invalid sample ID, label, birth, or provenance");
      }
      auto prev=i ? c.nodes[i-1].id : (c.closed ? c.nodes.back().id : 0);
      auto next=i+1<c.nodes.size() ? c.nodes[i+1].id :
                (c.closed ? c.nodes.front().id : 0);
      if (n.prev!=prev || n.next!=next) Fail("nonreciprocal or unordered connectivity");
      Vec wrapped;
      Locate(n.x, wrapped);
      if (!mesh->three_d && n.x[2]!=c.nodes.front().x[2]) {
        Fail("nonplanar curve in a 2D mesh");
      }
    }
    for (const auto &n : c.nodes) {
      if (n.anchor) {
        if (n.parent_left || n.parent_right || n.placement) {
          Fail("invalid anchor provenance");
        }
      } else if (!by_id.count(n.parent_left) || !by_id.count(n.parent_right) ||
                 n.parent_left==n.parent_right || n.parent_left>=n.id ||
                 n.parent_right>=n.id || !n.placement) {
        Fail("invalid insertion parents/method");
      } else {
        const auto &left=*by_id.at(n.parent_left), &right=*by_id.at(n.parent_right);
        Real upper=c.closed && right.a==0 ? c.period : right.a;
        if (!(left.a<n.a && n.a<upper) ||
            std::max(left.birth,right.birth)>n.birth) {
          Fail("insertion does not lie in its recorded parent material interval");
        }
      }
    }
    count+=c.nodes.size();
  }
#if MPI_PARALLEL_ENABLED
  std::uint64_t total;
  MPI_Allreduce(&count,&total,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
  count=total;
#endif
  if (count>max_nodes) Fail("global node count exceeds budget");
}

void MaterialCurves::Adapt(Real birth_time, bool initial, Real affected_time) {
  if (affected_time<0) affected_time=birth_time;
  Real spacing=max_spacing>0 ? std::min(max_spacing,mesh_spacing_) : mesh_spacing_;
  Real tolerance=chord_error>0 ? chord_error : 0.02*spacing;
  // Transactions iterate collectively; cap rejection leaves every curve unchanged.
  for (int iteration=0; iteration<64; ++iteration) {
    struct Split {std::size_t curve, edge; Vec x; Real a; int method;};
    std::vector<Split> splits;
    std::uint64_t local_count=0;
    for (std::size_t ci=0; ci<curves.size(); ++ci) {
      auto &c=curves[ci];
      local_count+=c.nodes.size();
      std::size_t n=c.nodes.size(), edges=c.closed ? n : n-1;
      for (std::size_t i=0; i<edges; ++i) {
        std::size_t j=(i+1)%n;
        Vec edge=Edge(c,i);
        Real length=Norm(edge);
        if (!(length>0) || !std::isfinite(length)) {
          if (!c.reason) {c.reason=3; c.unresolved_since=birth_time;}
          continue;
        }
        if (!adaptive) continue;  // fixed sampling still records degeneracy
        Real da=j ? c.nodes[j].a-c.nodes[i].a : c.period-c.nodes[i].a;
        Vec left{}, right{};
        for (int d=0; d<3; ++d) left[d]=right[d]=edge[d]/da;
        Real angle=0;
        if (i || c.closed) {
          std::size_t p=(i+n-1)%n;
          Vec v=Edge(c,p);
          Real dp=i ? c.nodes[i].a-c.nodes[p].a : c.period-c.nodes[p].a;
          Real vp=Norm(v);
          if (vp>0) {
            angle=std::max(angle,std::acos(std::clamp(Dot(v,edge)/(vp*length),
                                                     Real(-1),Real(1))));
            for (int d=0; d<3; ++d) left[d]=(da*v[d]/dp+dp*edge[d]/da)/(dp+da);
          }
        }
        if (j+1<n || c.closed) {
          Vec v=Edge(c,j);
          Real dn=j+1<n ? c.nodes[j+1].a-c.nodes[j].a : c.period-c.nodes[j].a;
          Real vn=Norm(v);
          if (vn>0) {
            angle=std::max(angle,std::acos(std::clamp(Dot(edge,v)/(length*vn),
                                                     Real(-1),Real(1))));
            for (int d=0; d<3; ++d) right[d]=(dn*edge[d]/da+da*v[d]/dn)/(da+dn);
          }
        }
        Vec correction;
        for (int d=0; d<3; ++d) correction[d]=da*(left[d]-right[d])/8;
        Real error=Norm(correction);
        bool known_polygon=initial && !circle_seeds_.count(c.id);
        if (known_polygon && angle>max_angle && !c.reason) {
          c.reason=4; c.unresolved_since=birth_time;  // true nonsmooth seed corner
        }
        if (length<=spacing && (known_polygon || angle<=max_angle) &&
            (known_polygon || error<=tolerance)) continue;
        if (!initial && iteration==0 &&
            (length>2*spacing || angle>2*max_angle || error>4*tolerance) && !c.reason) {
          c.reason=4;
          c.unresolved_since=std::max(Real(0),affected_time);
        }
        Real a=c.nodes[i].a+0.5*da;
        if (!(a>c.nodes[i].a && a<(j ? c.nodes[j].a : c.period)) ||
            c.next_sequence>MaxSequence(curve_id_bits)) {
          if (!c.reason) {c.reason=2; c.unresolved_since=birth_time;}
          continue;
        }
        int method=2;
        if (known_polygon || error>0.25*length ||
            std::abs(Dot(correction,edge))>0.25*length*length) {
          correction={}; method=1;
        }
        Vec x;
        for (int d=0; d<3; ++d) x[d]=c.nodes[i].x[d]+0.5*edge[d]+correction[d];
        if (initial && circle_seeds_.count(c.id)) {
          const auto &seed=circle_seeds_.at(c.id);
          Real phase=a/seed.radius;
          for (int d=0; d<3; ++d) {
            x[d]=seed.center[d]+seed.radius*(std::cos(phase)*seed.e1[d]+
                                           std::sin(phase)*seed.e2[d]);
          }
          method=3;
        }
        splits.push_back({ci,i,x,a,method});
      }
    }
    std::uint64_t requested=splits.size(), global_count=local_count;
#if MPI_PARALLEL_ENABLED
    std::uint64_t global_requested;
    MPI_Allreduce(&requested,&global_requested,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
    requested=global_requested;
    MPI_Allreduce(&local_count,&global_count,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
#endif
    if (!requested) break;
    if (requested>max_nodes-global_count) {
      for (const auto &s : splits) {
        auto &c=curves[s.curve];
        if (!c.reason) {c.reason=1; c.unresolved_since=birth_time;}
      }
      break;
    }
    // Merge each accepted split batch in material order: O(nodes + splits).
    std::size_t cursor=0;
    for (std::size_t ci=0; ci<curves.size(); ++ci) {
      auto &c=curves[ci];
      std::size_t begin=cursor;
      while (cursor<splits.size() && splits[cursor].curve==ci) ++cursor;
      if (begin==cursor) continue;
      if (cursor-begin > MaxSequence(curve_id_bits)+1-c.next_sequence) {
        if (!c.reason) {c.reason=2; c.unresolved_since=birth_time;}
        continue;
      }
      std::vector<Node> updated;
      updated.reserve(c.nodes.size()+cursor-begin);
      std::size_t split=begin;
      for (std::size_t i=0; i<c.nodes.size(); ++i) {
        updated.push_back(c.nodes[i]);
        if (split<cursor && splits[split].edge==i) {
          const auto &s=splits[split++];
          Node node;
          node.id=NewID(c);
          node.a=s.a; node.birth=birth_time; node.anchor=0; node.placement=s.method;
          node.x=node.old=s.x;
          node.parent_left=c.nodes[i].id;
          node.parent_right=c.nodes[(i+1)%c.nodes.size()].id;
          updated.push_back(node);
          ++c.inserted;
        }
      }
      c.nodes.swap(updated);
      Links(c);
    }
    if (iteration==63) {
      for (auto &c : curves) {
        if (!c.reason) {c.reason=2; c.unresolved_since=birth_time;}
      }
    }
  }
  Audit();
}

void MaterialCurves::CompleteStep() {
  Adapt(mesh->time, false, mesh->time-mesh->dt);
  if (shock_enabled) ObserveShocks();
  if (output_dcycle>0 && mesh->ncycle%output_dcycle==0) {
    std::ostringstream path;
    path << basename << ".curves." << std::setw(8) << std::setfill('0')
         << mesh->ncycle << ".txt";
    Snapshot(path.str());
  }
}
}  // namespace material_curves
