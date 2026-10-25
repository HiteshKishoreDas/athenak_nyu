// Copyright(C) 2026 the AthenaK contributors
// Licensed under the 3-clause BSD License (the "LICENSE").
// Versioned material-curve snapshots and shared-fluid checkpoint sidecars.
#include "material_curves.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <string>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "globals.hpp"

namespace material_curves {
void MaterialCurves::ThermoHistory(
    const std::vector<particles::TracerField> &fields,
    std::vector<int> &integers, std::vector<Real> &reals) {
  static_assert(sizeof(int)==sizeof(NodeID), "THP requires 32-bit integer slots");
  std::vector<Node *> nodes;
  for (auto &c : curves) for (auto &n : c.nodes) nodes.push_back(&n);
  std::sort(nodes.begin(),nodes.end(),[](const Node *a,const Node *b) {
    return a->id<b->id;
  });
  std::vector<Real> values;
  Sample(nodes,nullptr,fields,&values);
  bool shock_fields=false;
  for (const auto &field : fields) shock_fields |= particles::IsShockField(field.kind);
  std::vector<ShockObservation> observations;
  if (shock_fields) {
    if (!shock_enabled || shock_generation!=static_cast<std::uint64_t>(mesh->ncycle) ||
        shock_time!=mesh->time) Fail("stale curve shock output generation");
    Sample(nodes,&observations);
  }
  integers.resize(4*nodes.size());
  const std::size_t width=7+fields.size();  // time, x[3], birth, a, level, fields
  reals.resize(width*nodes.size());
  for (std::size_t p=0; p<nodes.size(); ++p) {
    const auto &node=*nodes[p];
    Vec wrapped;
    int gid=Locate(node.x,wrapped);
    integers[4*p]=mesh->ncycle;
    integers[4*p+1]=node.id;
    integers[4*p+2]=CurveID(node.id,curve_id_bits);
    integers[4*p+3]=gid;
    reals[width*p]=mesh->time;
    for (int d=0; d<3; ++d) reals[width*p+1+d]=node.x[d];
    reals[width*p+4]=node.birth;
    reals[width*p+5]=node.a;
    reals[width*p+6]=mesh->lloc_eachmb[gid].level;
    for (std::size_t n=0; n<fields.size(); ++n) {
      reals[width*p+7+n]=particles::IsShockField(fields[n].kind) ?
          particles::EvaluateShockField(fields[n].kind,observations[p].observation) :
          values[fields.size()*p+n];
    }
  }
}

namespace {
constexpr std::uint64_t kOffset = 14695981039346656037ULL;
constexpr std::uint64_t kPrime = 1099511628211ULL;
constexpr std::uint64_t kCurveCap = 100000;
constexpr std::uint64_t kMeshLeafCap = 10000000;
const char *kCurveColumns = "curve_columns id:i32 closed:i32 winding:i64[3] "
    "period:real unresolved_since:real reason:i32 next_sequence:u64 inserted:u64 "
    "migrations:u64 samples:u64 count:u64 length:real";
const char *kNodeColumns = "node_columns id:i32 prev:i32 next:i32 parent_left:i32 "
    "parent_right:i32 a:real birth:real anchor:i32 placement:i32 x:real[3] image:i64[3] "
    "edge_length:real stretch:real curvature:real curvature_mask:i32";

void Require(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error("material curves: " + message);
}
void Hash(std::uint64_t &value, const char *bytes, std::size_t size) {
  for (std::size_t i = 0; i < size; ++i) {
    value ^= static_cast<unsigned char>(bytes[i]);
    value *= kPrime;
  }
}
std::pair<std::uint64_t, std::uint64_t> FluidIdentity(const std::string &path) {
  std::ifstream input(path, std::ios::binary);
  Require(input.is_open(), "cannot open fluid checkpoint " + path);
  std::uint64_t size = 0, hash = kOffset;
  char buffer[65536];
  while (input.read(buffer, sizeof(buffer)) || input.gcount() != 0) {
    auto count = static_cast<std::size_t>(input.gcount());
    Require(size <= UINT64_MAX - count, "fluid file size overflow");
    size += count;
    Hash(hash, buffer, count);
  }
  Require(input.eof(), "cannot read fluid checkpoint " + path);
  return {size, hash};
}
std::array<Real, 6> Domain(const Mesh *mesh) {
  const auto &s = mesh->mesh_size;
  return {s.x1min, s.x1max, s.x2min, s.x2max, s.x3min, s.x3max};
}
std::ostringstream Stream() {
  std::ostringstream stream;
  stream.imbue(std::locale::classic());
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10);
  return stream;
}
Vec Difference(const Curve &curve, std::size_t i, const std::array<Real, 6> &d) {
  Vec result{};
  if (!curve.closed && i + 1 == curve.nodes.size()) return result;
  std::size_t j = (i + 1) % curve.nodes.size();
  for (int axis = 0; axis < 3; ++axis) {
    result[axis] = curve.nodes[j].x[axis] - curve.nodes[i].x[axis];
    if (j == 0) result[axis] += curve.winding[axis] * (d[2*axis+1]-d[2*axis]);
  }
  return result;
}
Real Norm(const Vec &x) { return std::hypot(x[0], x[1], x[2]); }
std::string Encode(const Curve &curve, const std::array<Real, 6> &domain) {
  auto out = Stream();
  std::vector<Vec> edges(curve.nodes.size());
  Real length = 0;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    edges[i] = Difference(curve, i, domain);
    length += Norm(edges[i]);
  }
  out << "curve " << curve.id << ' ' << curve.closed;
  for (auto winding : curve.winding) out << ' ' << winding;
  out << ' ' << curve.period << ' ' << curve.unresolved_since << ' ' << curve.reason
      << ' ' << curve.next_sequence << ' ' << curve.inserted << ' ' << curve.migrations
      << ' ' << curve.samples << ' ' << curve.nodes.size() << ' ' << length << '\n';
  for (std::size_t i = 0; i < curve.nodes.size(); ++i) {
    const auto &node = curve.nodes[i];
    out << "node " << node.id << ' ' << node.prev << ' ' << node.next << ' '
        << node.parent_left << ' ' << node.parent_right << ' ' << node.a << ' '
        << node.birth << ' ' << node.anchor << ' ' << node.placement;
    for (auto x : node.x) out << ' ' << x;
    for (int axis = 0; axis < 3; ++axis) {
      long double image = std::floor((static_cast<long double>(node.x[axis])
          - domain[2*axis]) / (domain[2*axis+1] - domain[2*axis]));
      Require(std::isfinite(image) && image >= INT64_MIN && image <= INT64_MAX,
              "periodic image diagnostic overflow");
      out << ' ' << static_cast<std::int64_t>(image);
    }
    Real edge_length = Norm(edges[i]), stretch = 0;
    if (curve.closed || i + 1 < curve.nodes.size()) {
      Real da = (i+1 == curve.nodes.size() ? curve.period : curve.nodes[i+1].a)
          - node.a;
      stretch = edge_length / da;
    }
    int mask = 0;
    Real curvature = std::numeric_limits<Real>::quiet_NaN();
    if (!curve.closed && (i == 0 || i+1 == curve.nodes.size())) {
      mask = 1;  // endpoint: no two-sided estimator
    } else {
      const auto &left = edges[(i + edges.size()-1) % edges.size()];
      const auto &right = edges[i];
      Vec cross{left[1]*right[2]-left[2]*right[1],
                left[2]*right[0]-left[0]*right[2],
                left[0]*right[1]-left[1]*right[0]};
      Vec chord{left[0]+right[0], left[1]+right[1], left[2]+right[2]};
      Real denominator = Norm(left)*Norm(right)*Norm(chord);
      if (denominator > 0 && std::isfinite(denominator)) {
        curvature = 2*Norm(cross)/denominator;
        if (!std::isfinite(curvature)) mask = 2;
      } else {
        mask = 2;  // collapsed edge or reversing/degenerate three-point stencil
      }
    }
    if (mask != 0) curvature = std::numeric_limits<Real>::quiet_NaN();
    Require(std::isfinite(length) && std::isfinite(edge_length)
            && std::isfinite(stretch), "nonfinite geometry diagnostics");
    out << ' ' << edge_length << ' ' << stretch << ' ' << curvature << ' ' << mask
        << '\n';
  }
  return out.str();
}
std::vector<std::string> Tokens(const std::string &line) {
  std::istringstream stream(line);
  stream.imbue(std::locale::classic());
  std::vector<std::string> fields;
  std::string token;
  while (stream >> token) fields.push_back(token);
  return fields;
}
template<typename T> T Integer(const std::string &text) {
  T value{};
  auto result = std::from_chars(text.data(), text.data()+text.size(), value);
  Require(result.ec == std::errc() && result.ptr == text.data()+text.size(),
          "invalid exact-width integer: " + text);
  return value;
}
Real Number(const std::string &text, bool diagnostic = false) {
  std::istringstream stream(text);
  stream.imbue(std::locale::classic());
  Real value{};
  if (diagnostic && text == "nan") return std::numeric_limits<Real>::quiet_NaN();
  Require(static_cast<bool>(stream >> value) && stream.peek() == EOF
          && std::isfinite(value), "invalid real: " + text);
  return value;
}
Curve DecodeCurve(const std::vector<std::string> &f) {
  Curve curve;
  curve.id = Integer<NodeID>(f[1]);
  int closed = Integer<int>(f[2]);
  Require(closed == 0 || closed == 1, "invalid closed flag");
  curve.closed = closed != 0;
  for (int i = 0; i < 3; ++i) curve.winding[i] = Integer<std::int64_t>(f[3+i]);
  curve.period = Number(f[6]);
  curve.unresolved_since = Number(f[7]);
  curve.reason = Integer<std::int32_t>(f[8]);
  curve.next_sequence = Integer<std::uint64_t>(f[9]);
  curve.inserted = Integer<std::uint64_t>(f[10]);
  curve.migrations = Integer<std::uint64_t>(f[11]);
  curve.samples = Integer<std::uint64_t>(f[12]);
  return curve;
}

Node DecodeNode(const std::vector<std::string> &f, bool unsigned_ids = false) {
  auto id = [unsigned_ids](const std::string &text) {
    return unsigned_ids ? SignedID(Integer<std::uint32_t>(text)) : Integer<NodeID>(text);
  };
  Node node;
  node.id = id(f[1]);
  node.prev = id(f[2]);
  node.next = id(f[3]);
  node.parent_left = id(f[4]);
  node.parent_right = id(f[5]);
  node.a = Number(f[6]);
  node.birth = Number(f[7]);
  node.anchor = Integer<std::int32_t>(f[8]);
  node.placement = Integer<std::int32_t>(f[9]);
  for (int j = 0; j < 3; ++j) node.x[j] = Number(f[10+j]);
  node.old = node.x;
  return node;
}

struct Header {
  std::string kind, basename;
  int precision = 0, cycle = 0, adaptive = 0, dcycle = 0, version = 0;
  unsigned curve_id_bits = 8;
  Real time = 0, spacing = 0, spacing_cells = 1, angle = 0, error = 0;
  std::uint64_t count = 0, nodes = 0, cap = 0, fluid_size = 0, fluid_hash = 0;
  int amr = 0, check_interval = 0, refinement_interval = 0, prolong = 0;
  std::vector<std::int32_t> amr_clocks;
  std::array<Real, 6> domain{};
  std::array<std::int32_t, 6> bcs{};
  std::array<std::int32_t, 3> dimensions{};
};
class Reader {
 public:
  explicit Reader(const std::string &path, const Mesh *mesh = nullptr) :
      input_(path, std::ios::binary) {
    Require(input_.is_open(), "cannot open curve file " + path);
    auto fields = Fields(7, "MATERIAL_CURVES");
    Require((fields[1] == "3" || fields[1] == "4" || fields[1] == "5") && fields[2] == "ASCII"
            && fields[4] == (fields[1] == "5" ? "I32" : "U32")
            && fields[5] == "I64" && fields[6] == "I32",
            "unsupported curve schema; restart requires v3/v4/v5 compact IDs");
    header.precision = Integer<int>(fields[3]);
    Require(header.precision == 32 || header.precision == 64,
            "unsupported Real precision");
    header.version = Integer<int>(fields[1]);
    if (header.version>=4) {
      auto bits=Fields(2,"id_bits");
      header.curve_id_bits=Integer<unsigned>(bits[1]);
      Require(header.curve_id_bits<=30,"invalid curve_id_bits metadata");
    }
    fields = Fields(6, "state");
    header.kind = fields[1];
    Require(header.kind == "snapshot" || header.kind == "checkpoint", "invalid kind");
    header.time = Number(fields[2]);
    header.cycle = Integer<int>(fields[3]);
    header.count = Integer<std::uint64_t>(fields[4]);
    header.nodes = Integer<std::uint64_t>(fields[5]);
    Require(header.time >= 0 && header.cycle >= 0 &&
            header.count <= std::min<std::uint64_t>(
                kCurveCap,MaxCurveID(header.curve_id_bits)),
            "invalid time, cycle or curve count");
    fields = Fields(16, "domain");
    for (int i = 0; i < 6; ++i) header.domain[i] = Number(fields[i+1]);
    for (int i = 0; i < 6; ++i) header.bcs[i] = Integer<std::int32_t>(fields[i+7]);
    for (int i = 0; i < 3; ++i) {
      header.dimensions[i] = Integer<std::int32_t>(fields[i+13]);
      Require(header.domain[2*i+1] > header.domain[2*i]
              && header.dimensions[i] > 0, "invalid domain");
    }
    fields = Fields(8, "config");
    header.cap = Integer<std::uint64_t>(fields[1]);
    header.spacing = Number(fields[2]);
    header.spacing_cells = Number(fields[3]);
    header.angle = Number(fields[4]);
    header.error = Number(fields[5]);
    header.adaptive = Integer<int>(fields[6]);
    header.dcycle = Integer<int>(fields[7]);
    Require(header.cap > 0 && header.nodes <= header.cap && header.spacing >= 0
            && header.spacing_cells > 0 && header.angle > 0 && header.error >= 0
            && header.dcycle >= 0
            && (header.adaptive == 0 || header.adaptive == 1), "invalid config/caps");
    std::istringstream name(Line());
    std::string key;
    Require(static_cast<bool>(name >> key >> std::quoted(header.basename))
            && key == "basename" && (name >> std::ws).eof(), "invalid basename");
    fields = Fields(3, "fluid");
    header.fluid_size = Integer<std::uint64_t>(fields[1]);
    header.fluid_hash = Integer<std::uint64_t>(fields[2]);
    fields = Fields(6, "amr");
    header.amr = Integer<int>(fields[1]);
    header.check_interval = Integer<std::int32_t>(fields[2]);
    header.refinement_interval = Integer<std::int32_t>(fields[3]);
    header.prolong = Integer<int>(fields[4]);
    auto leaves = Integer<std::uint64_t>(fields[5]);
    Require((header.amr == 0 || header.amr == 1) && leaves <= kMeshLeafCap
            && (header.prolong == 0 || header.prolong == 1), "invalid AMR metadata");
    if (header.amr) {
      Require(header.check_interval > 0 && header.refinement_interval > 0 && leaves > 0,
              "invalid AMR intervals/count");
    } else {
      Require(leaves == 0 && header.check_interval == 0
              && header.refinement_interval == 0 && header.prolong == 0,
              "inactive AMR has state");
    }
    if (mesh) {
      Require(header.amr == static_cast<int>(mesh->adaptive), "AMR mode mismatch");
      Require(leaves == (mesh->adaptive ?
              static_cast<std::uint64_t>(mesh->nmb_total) : 0),
              "AMR leaf count mismatch");
      if (mesh->adaptive) {
        Require(mesh->pmr && header.check_interval == mesh->pmr->ncyc_check_amr
                && header.refinement_interval == mesh->pmr->refinement_interval
                && header.prolong == static_cast<int>(mesh->pmr->prolong_prims),
                "AMR configuration mismatch");
      }
    }
    header.amr_clocks.reserve(static_cast<std::size_t>(leaves));
    for (std::uint64_t i = 0; i < leaves; ++i) {
      fields = Fields(6, "leaf");
      LogicalLocation loc;
      loc.level = Integer<std::int32_t>(fields[1]);
      loc.lx1 = Integer<std::int32_t>(fields[2]);
      loc.lx2 = Integer<std::int32_t>(fields[3]);
      loc.lx3 = Integer<std::int32_t>(fields[4]);
      auto clock = Integer<std::int32_t>(fields[5]);
      Require(loc.level >= 0 && loc.lx1 >= 0 && loc.lx2 >= 0 && loc.lx3 >= 0
              && clock >= 0 && clock <= header.cycle, "invalid AMR leaf/clock");
      if (mesh) Require(loc == mesh->lloc_eachmb[i], "AMR fluid leaf identity mismatch");
      header.amr_clocks.push_back(clock);
    }
    std::string curve_columns=kCurveColumns, node_columns=kNodeColumns;
    if (header.version<5) {
      curve_columns.replace(curve_columns.find("id:i32"),6,"id:u32");
      for (const auto *name : {"id", "prev", "next", "parent_left", "parent_right"}) {
        std::string column=std::string(name)+":i32";
        node_columns.replace(node_columns.find(column),column.size(),
                             std::string(name)+":u32");
      }
    }
    if (header.version==3) {
      const auto pos=curve_columns.find("next_sequence:u64");
      curve_columns.replace(pos,17,"next_sequence:u32");
    }
    Require(Line() == curve_columns && Line() == node_columns,
            "unsupported column schema");
  }
  Curve Next() {
    Require(curves_ < header.count, "too many curves");
    auto f = Fields(15, "curve");
    Curve curve = DecodeCurve(f);
    auto count = Integer<std::uint64_t>(f[13]);
    Require(Number(f[14]) >= 0, "negative diagnostic length");
    Require(curve.id > 0 && curve.id <= MaxCurveID(header.curve_id_bits)
            && ids_.insert(curve.id).second,
            "duplicate or invalid curve ID");
    Require(count >= (curve.closed ? 3u : 2u) && count <= header.nodes-nodes_,
            "invalid node count");
    Require(curve.period > 0 && curve.next_sequence > 0
            && curve.next_sequence <= MaxSequence(header.curve_id_bits)+1
            && curve.reason >= 0
            && curve.reason <= 4 && ((curve.reason == 0 && curve.unresolved_since == -1)
                || (curve.reason != 0 && curve.unresolved_since >= 0
                    && curve.unresolved_since <= header.time)), "invalid curve metadata");
    for (auto winding : curve.winding) {
      Require(winding >= -(INT64_C(1) << 40) && winding <= (INT64_C(1) << 40),
              "winding exceeds reliable lift range");
    }
    if (header.dimensions[2] == 1) {
      Require(curve.winding[2] == 0, "2D curve has transverse winding");
    }
    if (!curve.closed) Require(curve.winding == std::array<std::int64_t, 3>{},
                               "open string has winding");
    curve.nodes.reserve(static_cast<std::size_t>(count));
    for (std::uint64_t i = 0; i < count; ++i) {
      f = Fields(20, "node");
      Node node = DecodeNode(f,header.version<5);
      for (int j = 0; j < 3; ++j) {
        long double lift = (static_cast<long double>(node.x[j])-header.domain[2*j])
            / (static_cast<long double>(header.domain[2*j+1])-header.domain[2*j]);
        Require(std::isfinite(lift) && std::abs(lift) <= (INT64_C(1) << 40),
                "node coordinate exceeds reliable lift range");
      }
      if (header.dimensions[2] == 1 && !curve.nodes.empty()) {
        Require(node.x[2] == curve.nodes.front().x[2], "2D curve is not planar");
      }
      for (int j = 0; j < 3; ++j) Integer<std::int64_t>(f[13+j]);
      Require(Number(f[16]) >= 0 && Number(f[17]) >= 0,
              "negative edge/stretch diagnostic");
      Real curvature = Number(f[18], true);
      int mask = Integer<int>(f[19]);
      Require(mask >= 0 && mask <= 2 && ((mask == 0 && std::isfinite(curvature)
              && curvature >= 0) || (mask != 0 && std::isnan(curvature))),
              "invalid curvature diagnostic/mask");
      Require(CurveID(node.id,header.curve_id_bits) == curve.id
              && (IDBits(node.id) & MaxSequence(header.curve_id_bits)) != 0
              && (IDBits(node.id) & MaxSequence(header.curve_id_bits)) < curve.next_sequence,
              "node ID namespace or allocator mismatch");
      Require(node.a >= 0 && node.a <= curve.period
              && (!curve.closed || node.a < curve.period)
              && node.birth >= 0 && node.birth <= header.time
              && (node.anchor == 0 || node.anchor == 1), "invalid node metadata");
      curve.nodes.push_back(node);
    }
    std::sort(curve.nodes.begin(), curve.nodes.end(),
              [](const Node &a, const Node &b) { return a.a < b.a; });
    Require(curve.nodes.front().a == 0
            && (curve.closed || curve.nodes.back().a == curve.period),
            "material interval endpoints missing");
    std::unordered_map<NodeID, const Node *> nodes;
    for (const auto &n : curve.nodes) {
      Require(nodes.emplace(n.id, &n).second, "duplicate node ID");
    }
    std::uint64_t inserted = 0;
    for (std::size_t i = 0; i < curve.nodes.size(); ++i) {
      const auto &n = curve.nodes[i];
      auto prev = i == 0 ? (curve.closed ? curve.nodes.back().id : 0)
                         : curve.nodes[i-1].id;
      auto next = i+1 == curve.nodes.size() ? (curve.closed ? curve.nodes.front().id : 0)
                                          : curve.nodes[i+1].id;
      Require(n.prev == prev && n.next == next
              && (i == 0 || curve.nodes[i-1].a < n.a),
              "nonreciprocal links or material ordering");
      if (n.anchor) {
        Require(n.parent_left == 0 && n.parent_right == 0 && n.birth == 0
                && n.placement == 0,
                "invalid anchor provenance");
      } else {
        ++inserted;
        auto left = nodes.find(n.parent_left), right = nodes.find(n.parent_right);
        Require((n.placement >= 1 && n.placement <= 3)
                && left != nodes.end() && right != nodes.end()
                && n.parent_left != n.id && n.parent_right != n.id
                && n.parent_left != n.parent_right
                && IDBits(n.parent_left) < IDBits(n.id)
                && IDBits(n.parent_right) < IDBits(n.id)
                && left->second->birth <= n.birth
                && right->second->birth <= n.birth, "invalid insertion provenance");
        Real lower = left->second->a, upper = right->second->a;
        if (curve.closed && upper == 0) upper = curve.period;
        Require(lower < n.a && n.a < upper,
                "insertion parents do not bracket material label");
      }
    }
    Require(curve.inserted == inserted, "insertion counter mismatch");
    nodes_ += count;
    ++curves_;
    return curve;
  }
  void Finish() {
    Require(curves_ == header.count && nodes_ == header.nodes, "count mismatch");
    auto fields = Tokens(Line(false));
    Require(fields.size() == 3 && fields[0] == "commit"
            && Integer<std::uint64_t>(fields[1]) == hash_
            && Integer<std::uint64_t>(fields[2]) == bytes_, "invalid commit/checksum");
    Require(input_.peek() == EOF && !input_.bad(), "trailing data or read error");
  }
  Header header;

 private:
  std::string Line(bool hashed = true) {
    // Fixed records are tiny; do not let a malformed file allocate unbounded strings.
    std::string line;
    for (int i = 0; i < 8192; ++i) {
      int c = input_.get();
      Require(c != EOF, "truncated curve record");
      if (c == '\n') {
        if (hashed) {
          Hash(hash_, line.data(), line.size());
          Hash(hash_, "\n", 1);
          bytes_ += line.size()+1;
        }
        return line;
      }
      line.push_back(static_cast<char>(c));
    }
    throw std::runtime_error("material curves: oversized curve record");
  }
  std::vector<std::string> Fields(std::size_t count, const std::string &key) {
    auto fields = Tokens(Line());
    Require(fields.size() == count && fields[0] == key, "invalid " + key + " record");
    return fields;
  }
  std::ifstream input_;
  std::uint64_t hash_ = kOffset, bytes_ = 0, curves_ = 0, nodes_ = 0;
  std::set<NodeID> ids_;
};

#if MPI_PARALLEL_ENABLED
void Send(const std::string &data, int rank) {
  std::uint64_t size = data.size();
  MPI_Send(&size, 1, MPI_UINT64_T, rank, 27141, MPI_COMM_WORLD);
  for (std::size_t offset = 0; offset < data.size(); offset += 65536) {
    int count = static_cast<int>(std::min<std::size_t>(65536, data.size()-offset));
    MPI_Send(data.data()+offset, count, MPI_CHAR, rank, 27142, MPI_COMM_WORLD);
  }
}
std::string Receive(int rank, std::uint64_t cap) {
  std::uint64_t size;
  MPI_Recv(&size, 1, MPI_UINT64_T, rank, 27141, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // Encoded records have fewer than 1024 bytes per node plus a short curve record.
  Require(cap <= (UINT64_MAX-8192)/1024 && size <= cap*1024+8192,
          "oversized MPI curve record");
  std::string data(static_cast<std::size_t>(size), '\0');
  for (std::size_t offset = 0; offset < data.size(); offset += 65536) {
    int count = static_cast<int>(std::min<std::size_t>(65536, data.size()-offset));
    MPI_Recv(&data[offset], count, MPI_CHAR, rank, 27142, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  }
  return data;
}
#endif
void ValidateHeader(const Header &h, const MaterialCurves &m,
                    const std::pair<std::uint64_t, std::uint64_t> &identity) {
  Require(h.curve_id_bits == m.curve_id_bits,"restart curve_id_bits mismatch");
  Require(h.kind == "checkpoint" && h.precision == 8*sizeof(Real),
          "checkpoint kind or precision mismatch");
  Require(h.time == m.mesh->time && h.cycle == m.mesh->ncycle,
          "fluid/curve checkpoint time or cycle mismatch");
  Require(h.fluid_size == identity.first && h.fluid_hash == identity.second,
          "fluid/curve checkpoint content identity mismatch");
  Require(h.domain == Domain(m.mesh), "checkpoint domain mismatch");
  for (int i = 0; i < 6; ++i) {
    Require(h.bcs[i] == static_cast<int>(m.mesh->mesh_bcs[i]), "boundary mismatch");
  }
  Require(h.dimensions == std::array<std::int32_t, 3>{m.mesh->mesh_indcs.nx1,
          m.mesh->mesh_indcs.nx2, m.mesh->mesh_indcs.nx3}, "mesh dimension mismatch");
  Require(h.nodes <= m.max_nodes && h.spacing == m.max_spacing
          && h.spacing_cells == m.spacing_cells && h.angle == m.max_angle
          && h.error == m.chord_error
          && h.adaptive == static_cast<int>(m.adaptive), "restart config/cap mismatch");
}
void WriteFile(const MaterialCurves &m, const std::string &path,
               const std::string &fluid) {
  const int rank = global_variable::my_rank, ranks = global_variable::nranks;
  std::uint64_t counts[2]{static_cast<std::uint64_t>(m.curves.size()), 0};
  for (const auto &curve : m.curves) counts[1] += curve.nodes.size();
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, counts, 2, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  Real times[2]{m.mesh->time, m.mesh->time};
  int cycles[2]{m.mesh->ncycle, m.mesh->ncycle};
  MPI_Allreduce(MPI_IN_PLACE, &times[0], 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &times[1], 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &cycles[0], 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &cycles[1], 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  Require(times[0] == times[1] && cycles[0] == cycles[1], "unsynchronized checkpoint");
#endif
  Require(counts[0] <= kCurveCap && counts[1] <= m.max_nodes, "output count cap");
  auto amr = Stream();
  std::uint64_t local_amr_hash = kOffset;
  if (m.mesh->adaptive) {
    Require(m.mesh->pmr && m.mesh->nmb_total > 0
            && static_cast<std::uint64_t>(m.mesh->nmb_total) <= kMeshLeafCap,
            "invalid AMR mesh leaf count");
    const auto &ref = *m.mesh->pmr;
    Require(ref.ncyc_check_amr > 0 && ref.refinement_interval > 0
            && ref.ncyc_since_ref.extent(0) ==
                static_cast<std::size_t>(m.mesh->nmb_total),
            "invalid AMR refinement state");
    amr << "amr 1 " << ref.ncyc_check_amr << ' ' << ref.refinement_interval << ' '
        << ref.prolong_prims << ' ' << m.mesh->nmb_total << '\n';
    auto amr_header = amr.str();
    Hash(local_amr_hash, amr_header.data(), amr_header.size());
    for (int i = 0; i < m.mesh->nmb_total; ++i) {
      const auto &loc = m.mesh->lloc_eachmb[i];
      int clock = ref.ncyc_since_ref(i);
      Require(clock >= 0 && clock <= m.mesh->ncycle, "invalid AMR refinement clock");
      auto row = Stream();
      row << "leaf " << loc.level << ' ' << loc.lx1 << ' ' << loc.lx2 << ' '
          << loc.lx3 << ' ' << clock << '\n';
      auto record = row.str();
      Hash(local_amr_hash, record.data(), record.size());
    }
  } else {
    amr << "amr 0 0 0 0 0\n";
    auto record = amr.str();
    Hash(local_amr_hash, record.data(), record.size());
  }
#if MPI_PARALLEL_ENABLED
  std::uint64_t amr_hash[2]{local_amr_hash, local_amr_hash};
  MPI_Allreduce(MPI_IN_PLACE, &amr_hash[0], 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &amr_hash[1], 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  Require(amr_hash[0] == amr_hash[1], "inconsistent replicated AMR clocks/tree");
#endif
  std::ofstream file;
  std::uint64_t hash = kOffset, bytes = 0;
  auto write = [&](const std::string &text) {
    file.write(text.data(), static_cast<std::streamsize>(text.size()));
    Require(static_cast<bool>(file), "curve output write failure: " + path);
    Hash(hash, text.data(), text.size());
    bytes += text.size();
  };
  if (rank == 0) {
    auto identity = fluid.empty() ? std::pair<std::uint64_t, std::uint64_t>{0, 0}
                                  : FluidIdentity(fluid);
    file.open(path+".tmp", std::ios::binary | std::ios::trunc);
    Require(file.is_open(), "cannot create curve output " + path);
    auto out = Stream();
    out << "MATERIAL_CURVES 5 ASCII " << 8*sizeof(Real) << " I32 I64 I32\nid_bits "
        << m.curve_id_bits << "\nstate "
        << (fluid.empty() ? "snapshot" : "checkpoint") << ' ' << m.mesh->time << ' '
        << m.mesh->ncycle << ' ' << counts[0] << ' ' << counts[1] << "\ndomain";
    for (auto value : Domain(m.mesh)) out << ' ' << value;
    for (auto bc : m.mesh->mesh_bcs) out << ' ' << static_cast<int>(bc);
    out << ' ' << m.mesh->mesh_indcs.nx1 << ' ' << m.mesh->mesh_indcs.nx2 << ' '
        << m.mesh->mesh_indcs.nx3 << "\nconfig " << m.max_nodes << ' ' << m.max_spacing
        << ' ' << m.spacing_cells << ' ' << m.max_angle << ' ' << m.chord_error
        << ' ' << m.adaptive << ' '
        << m.output_dcycle << "\nbasename " << std::quoted(m.basename) << "\nfluid "
        << identity.first << ' ' << identity.second << '\n';
    write(out.str());
    write(amr.str());
    if (m.mesh->adaptive) {
      for (int i = 0; i < m.mesh->nmb_total; ++i) {
        const auto &loc = m.mesh->lloc_eachmb[i];
        auto row = Stream();
        row << "leaf " << loc.level << ' ' << loc.lx1 << ' ' << loc.lx2 << ' '
            << loc.lx3 << ' ' << m.mesh->pmr->ncyc_since_ref(i) << '\n';
        write(row.str());
      }
    }
    write(std::string(kCurveColumns)+"\n"+kNodeColumns+"\n");
  }
  for (int source = 0; source < ranks; ++source) {
    std::uint64_t count = rank == source ? m.curves.size() : 0;
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(&count, 1, MPI_UINT64_T, source, MPI_COMM_WORLD);
#endif
    for (std::uint64_t i = 0; i < count; ++i) {
      if (rank == source) {
        Require(m.curves[i].id % ranks == static_cast<std::uint64_t>(rank),
                "curve home ownership mismatch");
        auto encoded = Encode(m.curves[i], Domain(m.mesh));
        if (rank == 0) write(encoded);
#if MPI_PARALLEL_ENABLED
        else
          Send(encoded, 0);
#endif
      } else if (rank == 0) {
#if MPI_PARALLEL_ENABLED
        write(Receive(source, m.max_nodes));
#endif
      }
    }
  }
  if (rank == 0) {
    file << "commit " << hash << ' ' << bytes << '\n';
    file.flush();
    Require(static_cast<bool>(file), "curve commit write failure: " + path);
    file.close();
    Require(!file.fail(), "curve close failure: " + path);
    Require(std::rename((path+".tmp").c_str(), path.c_str()) == 0,
            "cannot atomically publish curve output " + path);
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
}  // namespace

void MaterialCurves::Snapshot(const std::string &path) const {
  try {
    Audit();
    WriteFile(*this, path, "");
  } catch (const std::exception &error) {
    Fail(error.what());
  }
}
void MaterialCurves::Save(const std::string &fluid_path) const {
  try {
    Audit();
    WriteFile(*this, fluid_path+".curves", fluid_path);
  } catch (const std::exception &error) {
    Fail(error.what());
  }
}
void MaterialCurves::Load(const std::string &fluid_path) {
  try {
    const int rank = global_variable::my_rank, ranks = global_variable::nranks;
    const std::string path = fluid_path+".curves";
    std::vector<Curve> restored;
    std::vector<std::int32_t> clocks;
    if (rank == 0) {
      // Full-file preflight holds only one curve and a small set of curve IDs.
      Reader check(path, mesh);
      ValidateHeader(check.header, *this, FluidIdentity(fluid_path));
      for (std::uint64_t i = 0; i < check.header.count; ++i) check.Next();
      check.Finish();
      clocks = std::move(check.header.amr_clocks);
      Reader input(path, mesh);
      ValidateHeader(input.header, *this, FluidIdentity(fluid_path));
      for (std::uint64_t i = 0; i < input.header.count; ++i) {
        auto curve = input.Next();
        int home = static_cast<int>(curve.id % ranks);
        if (home == 0) restored.push_back(std::move(curve));
#if MPI_PARALLEL_ENABLED
        else
          Send(Encode(curve, input.header.domain), home);
#endif
      }
      input.Finish();
      Require(input.header.amr_clocks == clocks, "AMR clocks changed during read");
#if MPI_PARALLEL_ENABLED
      for (int home = 1; home < ranks; ++home) Send("", home);
#endif
    } else {
#if MPI_PARALLEL_ENABLED
      while (true) {
        auto record = Receive(0, max_nodes);
        if (record.empty()) break;
        // Root has independently audited these bounded, canonical transfer records.
        // Only rank0 reads global data; transfers contain one complete curve.
        std::istringstream lines(record);
        std::string line;
        std::getline(lines, line);
        auto f = Tokens(line);
        Curve c = DecodeCurve(f);
        auto count = Integer<std::uint64_t>(f[13]);
        c.nodes.reserve(static_cast<std::size_t>(count));
        for (std::uint64_t i = 0; i < count; ++i) {
          std::getline(lines, line);
          f = Tokens(line);
          Node n = DecodeNode(f);
          c.nodes.push_back(n);
        }
        Require(c.id % ranks == static_cast<std::uint64_t>(rank), "load ownership");
        restored.push_back(std::move(c));
      }
#endif
    }
#if MPI_PARALLEL_ENABLED
    std::uint64_t leaf_count = rank == 0 ? clocks.size() : 0;
    MPI_Bcast(&leaf_count, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    Require(leaf_count == (mesh->adaptive ?
            static_cast<std::uint64_t>(mesh->nmb_total) : 0),
            "distributed AMR leaf count mismatch");
    if (rank != 0) clocks.resize(static_cast<std::size_t>(leaf_count));
    MPI_Bcast(clocks.data(), static_cast<int>(leaf_count), MPI_INT32_T, 0,
              MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (mesh->adaptive) {
      for (int i = 0; i < mesh->nmb_total; ++i) mesh->pmr->ncyc_since_ref(i) = clocks[i];
    }
    curves.swap(restored);  // no live state changes before collective preflight success
    RebuildMesh();
    Audit();
  } catch (const std::exception &error) {
    Fail(error.what());
  }
}
}  // namespace material_curves
