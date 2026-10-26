"""Read AthenaK material-curve v2/v3/v4/v5 ASCII snapshots/checkpoint sidecars.

Current v5 IDs/links are signed32 with a declared curve-bit split; negative IDs
are valid and zero is missing. Previous v4 IDs are unsigned32; v3 uses unsigned 8/24.
Historical v2 IDs
retain their exact unsigned 64-bit values. IDs never pass through floats; winding/images
are signed 64-bit and flags signed 32-bit. Reals declare IEEE width 32 or 64,
written with max_digits10; ASCII has no endian dependence. Continuous lifts and
integer winding are authoritative. Labels a measure initial arc length (exact
circle circumference, initial
polygon arc length for line/polyline seeds),
so edge stretch is current chord length / initial material interval da. Curvature
is the three-point circumcircle estimator; mask 1 denotes open endpoints, mask 2
collapsed/reversing/nonfinite stencils, and masked curvature is nan. Curve reason
and unresolved_since carry unresolved quality independently of the stencil mask.

The FNV-1a 64-bit commit checksum covers every byte before the final commit row.
Checkpoints also fingerprint the complete shared fluid file (byte count and FNV).
Mandatory v2 AMR metadata preserves global logical leaf refinement clocks and
refinement scheduling options; v1 is explicitly rejected. Analysis reading uses
a separate max_mesh_leaves cap, independent of the material-node budget.
This detects accidental corruption; it is not a cryptographic authenticity check.
Snapshots are atomically replaced. Cycle-derived filenames govern continuation;
there is no append/history stream. Birth/provenance prevents invented pre-birth
history. Placement 0 is an anchor, 1 linear fallback, 2 limited cubic,
and 3 evaluation of a known analytic initial seed. At a completed step RK scratch
old is restored from x, and sample_owner
is discarded. Restart allows a larger node cap and changed output preferences,
but geometry thresholds/adaptive mode, precision/domain/boundaries must match.
Thresholds are inactive when adaptive is false; reason 0 then records structural
and degeneracy validity, without certifying geometric tolerances. Scientific
convergence and error budgets are reported by external validation sweeps.

read() accepts shuffled node rows, sorts by immutable a and independently audits
namespace, allocator, links, endpoints, ancestry, births and counters. geometry()
reconstructs closure using winding rather than a nearest-image shortcut. Root-only
C++ checkpoint preflight streams one curve at a time, then redistributes by ID;
this convenience analysis reader holds up to max_nodes nodes on its caller.
"""

import argparse
import json
import math
import re
import shlex

OFFSET = 14695981039346656037
PRIME = 1099511628211
MASK64 = (1 << 64) - 1
CURVE_COLUMNS = ("curve_columns id:u64 closed:i32 winding:i64[3] "
                 "period:real unresolved_since:real reason:i32 next_sequence:u64 "
                 "inserted:u64 migrations:u64 samples:u64 count:u64 length:real")
NODE_COLUMNS = ("node_columns id:u64 prev:u64 next:u64 parent_left:u64 "
                "parent_right:u64 a:real birth:real anchor:i32 placement:i32 x:real[3] "
                "image:i64[3] edge_length:real stretch:real curvature:real "
                "curvature_mask:i32")


def node_curve_id(node_id, sequence_bits=None, *, curve_id_bits=None, signed=False):
    """Decode declared current splits, or historical fixed-width identities."""
    if curve_id_bits is not None:
        if type(curve_id_bits) is not int or not 0 <= curve_id_bits <= 30:
            raise ValueError('invalid curve_id_bits')
        sequence_bits,width = 32-curve_id_bits,32
    else:
        if sequence_bits is None:
            sequence_bits = 24 if type(node_id) is int and node_id < 1 << 32 else 32
        if sequence_bits not in (24,32):
            raise ValueError('invalid historical ID encoding')
        width = 32 if sequence_bits == 24 else 64
    if signed:
        if type(node_id) is not int or not -(1 << 31) <= node_id < 1 << 31:
            raise ValueError('invalid signed32 curve node ID')
        node_id &= (1 << 32)-1
    if (type(node_id) is not int or not 0 < node_id < 1 << width or
            not node_id & ((1 << sequence_bits)-1)):
        raise ValueError('invalid curve node ID namespace/width')
    curve = 1 if curve_id_bits == 0 else node_id >> sequence_bits
    if not curve:
        raise ValueError('invalid curve node ID namespace/width')
    return curve


def _require(condition, message):
    if not condition:
        raise ValueError("material curves: " + message)


def _integer(value, bits=64, signed=False):
    pattern = r"-?[0-9]+" if signed else r"[0-9]+"
    _require(re.fullmatch(pattern, value) is not None, "invalid integer " + value)
    result = int(value)
    lower, upper = (-(1 << (bits-1)), (1 << (bits-1))-1) if signed else (0, (1 << bits)-1)
    _require(lower <= result <= upper, "integer width overflow " + value)
    return result


def _real(value, diagnostic=False):
    if diagnostic and value == "nan":
        return math.nan
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError("material curves: invalid real " + value) from error
    _require(math.isfinite(number), "nonfinite real " + value)
    return number


def fingerprint(path):
    """Stream the same byte-count/FNV identity used by the C++ writer."""
    size, checksum = 0, OFFSET
    with open(path, "rb") as stream:
        while True:
            data = stream.read(65536)
            if not data:
                break
            size += len(data)
            for byte in data:
                checksum = ((checksum ^ byte) * PRIME) & MASK64
    return size, checksum


class _Reader:
    def __init__(self, stream):
        self.stream = stream
        self.hash = OFFSET
        self.size = 0

    def line(self, hashed=True):
        data = self.stream.readline(8193)
        _require(data.endswith(b"\n") and len(data) <= 8192, "truncated/oversized record")
        if hashed:
            for byte in data:
                self.hash = ((self.hash ^ byte) * PRIME) & MASK64
            self.size += len(data)
        try:
            return data[:-1].decode("ascii")
        except UnicodeDecodeError as error:
            raise ValueError("material curves: non-ASCII record") from error

    def fields(self, count, key, hashed=True):
        fields = self.line(hashed).split()
        _require(len(fields) == count and fields[0] == key, "invalid " + key + " row")
        return fields


def _audit(curve, time, domain, dimensions, sequence_bits=32, curve_id_bits=None, signed=False):
    nodes = curve["nodes"]
    nodes.sort(key=lambda node: node["a"])
    _require(nodes[0]["a"] == 0 and
             (curve["closed"] or nodes[-1]["a"] == curve["period"]),
             "missing material endpoints")
    by_id = {node["id"]: node for node in nodes}
    _require(len(by_id) == len(nodes), "duplicate node ID")
    _require(all(abs(w) <= 2**40 for w in curve["winding"]),
             "winding exceeds reliable lift range")
    _require(dimensions[2] != 1 or curve["winding"][2] == 0,
             "2D curve has transverse winding")
    id_mask = (1 << 32)-1 if signed else MASK64
    inserted = 0
    for i, node in enumerate(nodes):
        _require(all(abs((node["x"][axis]-domain[2*axis]) /
                         (domain[2*axis+1]-domain[2*axis])) <= 2**40
                     for axis in range(3)), "node coordinate exceeds reliable lift range")
        _require(dimensions[2] != 1 or node["x"][2] == nodes[0]["x"][2],
                 "2D curve is not planar")
        _require(node_curve_id(node["id"],sequence_bits,curve_id_bits=curve_id_bits,signed=signed) == curve["id"] and
                 0 < (node["id"] & ((1 << sequence_bits)-1)) < curve["next_sequence"],
                 "node namespace/allocator mismatch")
        _require(0 <= node["a"] <= curve["period"] and
                 (not curve["closed"] or node["a"] < curve["period"]) and
                 (i == 0 or nodes[i-1]["a"] < node["a"]), "invalid material ordering")
        previous = nodes[i-1]["id"] if i or curve["closed"] else 0
        following = (nodes[(i+1) % len(nodes)]["id"]
                     if i+1 < len(nodes) or curve["closed"] else 0)
        _require(node["prev"] == previous and node["next"] == following, "invalid links")
        _require(node["anchor"] in (0, 1) and 0 <= node["birth"] <= time,
                 "invalid birth/anchor")
        if node["anchor"]:
            _require(node["parent_left"] == node["parent_right"] == 0 and
                     node["birth"] == 0 and node["placement"] == 0,
                     "invalid anchor provenance")
        else:
            inserted += 1
            left, right = node["parent_left"], node["parent_right"]
            _require(node["placement"] in (1, 2, 3) and left in by_id and right in by_id
                     and left != right and
                     (left & id_mask) < (node["id"] & id_mask) and
                     (right & id_mask) < (node["id"] & id_mask) and
                     max(by_id[left]["birth"], by_id[right]["birth"]) <= node["birth"],
                     "invalid insertion provenance")
            lower, upper = by_id[left]["a"], by_id[right]["a"]
            if curve["closed"] and upper == 0:
                upper = curve["period"]
            _require(lower < node["a"] < upper,
                     "insertion parents do not bracket material label")
    _require(inserted == curve["inserted"], "insertion count mismatch")


def geometry(curve, domain):
    """Reconstruct lifted edges and diagnostics independently of stored columns."""
    nodes = curve["nodes"]
    edges, lengths, stretch = [], [], []
    for i, node in enumerate(nodes):
        if not curve["closed"] and i+1 == len(nodes):
            edge, da = [0.0, 0.0, 0.0], 0.0
        else:
            next_node = nodes[(i+1) % len(nodes)]
            edge = [next_node["x"][j] - node["x"][j] for j in range(3)]
            if i+1 == len(nodes):
                edge = [edge[j] + curve["winding"][j] *
                        (domain[2*j+1]-domain[2*j]) for j in range(3)]
                da = curve["period"] - node["a"]
            else:
                da = next_node["a"] - node["a"]
        edges.append(edge)
        lengths.append(math.hypot(*edge))
        stretch.append(lengths[-1]/da if da else 0.0)
    curvature, masks = [], []
    for i in range(len(nodes)):
        if not curve["closed"] and i in (0, len(nodes)-1):
            curvature.append(math.nan)
            masks.append(1)
            continue
        left, right = edges[i-1], edges[i]
        chord = [left[j]+right[j] for j in range(3)]
        cross = [left[1]*right[2]-left[2]*right[1],
                 left[2]*right[0]-left[0]*right[2],
                 left[0]*right[1]-left[1]*right[0]]
        denominator = math.hypot(*left)*math.hypot(*right)*math.hypot(*chord)
        value = 2*math.hypot(*cross)/denominator if denominator else math.nan
        mask = 0 if math.isfinite(value) else 2
        curvature.append(value if mask == 0 else math.nan)
        masks.append(mask)
    return {"edges": edges, "length": math.fsum(lengths), "edge_length": lengths,
            "stretch": stretch, "curvature": curvature, "curvature_mask": masks}


def read(path, max_nodes=1000000, fluid_path=None, max_mesh_leaves=1000000):
    """Read and audit a bounded file; optionally verify shared-fluid identity."""
    with open(path,"rb") as probe:
        binary=probe.read(16)==b'ATHKCRV'+b'\0'*9
    if binary:
        _require(fluid_path is None,'CRV snapshots are not restart checkpoints')
        from curve_binary import read as read_binary
        return read_binary(path,max_nodes)

    _require(isinstance(max_mesh_leaves, int) and max_mesh_leaves >= 0,
             "invalid mesh reader cap")
    _require(isinstance(max_nodes, int) and max_nodes >= 0, "invalid reader cap")
    with open(path, "rb") as stream:
        reader = _Reader(stream)
        f = reader.fields(7, "MATERIAL_CURVES")
        version = _integer(f[1],32)
        bits = 32 if version >= 3 else 64
        curve_bits = 8 if version == 3 else None
        sequence_bits = 24 if version == 3 else 32
        _require(version in (2,3,4,5) and f[2] == "ASCII" and
                 f[4:] == ["I32" if version == 5 else "U32" if version >= 3 else "U64", "I64", "I32"],
                 "unsupported version/encoding/widths")
        precision = _integer(f[3], 32)
        _require(precision in (32, 64), "unsupported Real precision")
        if version >= 4:
            curve_bits = _integer(reader.fields(2,"id_bits")[1],32)
            _require(curve_bits <= 30,'invalid curve_id_bits metadata')
            sequence_bits = 32-curve_bits
        max_curve = ((1 << curve_bits)-1 if curve_bits else 1) if version >= 3 else (1 << 32)-1
        f = reader.fields(6, "state")
        state = {"kind": f[1], "time": _real(f[2]), "cycle": _integer(f[3], 32, True),
                 "curve_count": _integer(f[4]), "node_count": _integer(f[5]),
                 "precision": precision,"version":version,
                 "id_layout":(bits,curve_bits,True) if version == 5 else
                 (bits,curve_bits if version>=3 else 32)}
        _require(state["kind"] in ("snapshot", "checkpoint") and state["time"] >= 0
                 and state["cycle"] >= 0 and state["curve_count"] <= min(max_curve,100000)
                 and state["node_count"] <= max_nodes, "invalid state/count cap")
        f = reader.fields(16, "domain")
        state["domain"] = [_real(value) for value in f[1:7]]
        state["bcs"] = [_integer(value, 32, True) for value in f[7:13]]
        state["dimensions"] = [_integer(value, 32, True) for value in f[13:16]]
        _require(all(state["domain"][2*i+1] > state["domain"][2*i] and
                     state["dimensions"][i] > 0 for i in range(3)), "invalid domain")
        f = reader.fields(8, "config")
        state["config"] = {"max_nodes": _integer(f[1]), "max_spacing": _real(f[2]),
                           "spacing_cells": _real(f[3]), "max_angle": _real(f[4]),
                           "chord_error": _real(f[5]),
                           "adaptive": _integer(f[6], 32, True),
                           "output_dcycle": _integer(f[7], 32, True)}
        c = state["config"]
        _require(c["max_nodes"] > 0 and state["node_count"] <= c["max_nodes"] and
                 c["max_spacing"] >= 0 and c["spacing_cells"] > 0 and c["max_angle"] > 0
                 and c["chord_error"] >= 0
                 and c["adaptive"] in (0, 1) and c["output_dcycle"] >= 0,
                 "invalid config")
        name = shlex.split(reader.line())
        _require(len(name) == 2 and name[0] == "basename", "invalid basename")
        state["basename"] = name[1]
        f = reader.fields(3, "fluid")
        state["fluid_identity"] = (_integer(f[1]), _integer(f[2]))
        f = reader.fields(6, "amr")
        amr = {"active": _integer(f[1], 32, True),
               "ncycle_check": _integer(f[2], 32, True),
               "refinement_interval": _integer(f[3], 32, True),
               "prolong_primitives": _integer(f[4], 32, True), "leaves": []}
        leaves = _integer(f[5])
        _require(amr["active"] in (0, 1) and amr["prolong_primitives"] in (0, 1)
                 and 0 <= leaves <= min(max_mesh_leaves, 10000000), "invalid AMR metadata/cap")
        if amr["active"]:
            _require(amr["ncycle_check"] > 0 and amr["refinement_interval"] > 0
                     and leaves > 0, "invalid AMR interval/count")
        else:
            _require(leaves == 0 and amr["ncycle_check"] == 0 and
                     amr["refinement_interval"] == 0 and amr["prolong_primitives"] == 0,
                     "inactive AMR has state")
        leaf_ids = set()
        for _ in range(leaves):
            f = reader.fields(6, "leaf")
            values = [_integer(value, 32, True) for value in f[1:]]
            identity = tuple(values[:4])
            _require(all(value >= 0 for value in values) and values[4] <= state["cycle"]
                     and identity not in leaf_ids, "invalid/duplicate AMR leaf or clock")
            leaf_ids.add(identity)
            amr["leaves"].append(dict(zip(("level", "lx1", "lx2", "lx3",
                                           "cycles_since_ref"), values)))
        state["amr"] = amr
        curve_columns = CURVE_COLUMNS
        node_columns = NODE_COLUMNS
        if version >= 3:
            curve_columns = curve_columns.replace('id:u64','id:u32')
            node_columns = node_columns.replace(':u64',':u32')
            if version == 5:
                curve_columns = curve_columns.replace('id:u32','id:i32')
                node_columns = node_columns.replace(':u32',':i32')
            if version == 3:
                curve_columns = curve_columns.replace('next_sequence:u64','next_sequence:u32')
        _require(reader.line() == curve_columns and reader.line() == node_columns,
                 "unsupported columns")
        curves, ids, total = [], set(), 0
        for _ in range(state["curve_count"]):
            f = reader.fields(15, "curve")
            curve = {"id": _integer(f[1],bits,version == 5), "closed": _integer(f[2], 32, True),
                     "winding": [_integer(value, 64, True) for value in f[3:6]],
                     "period": _real(f[6]), "unresolved_since": _real(f[7]),
                     "reason": _integer(f[8], 32, True), "next_sequence": _integer(f[9],64 if version >= 4 else bits),
                     "inserted": _integer(f[10]), "migrations": _integer(f[11]),
                     "samples": _integer(f[12]), "nodes": [], "length": _real(f[14])}
            count = _integer(f[13])
            _require(0 < curve["id"] <= max_curve and curve["id"] not in ids and
                     curve["closed"] in (0, 1) and curve["period"] > 0 and
                     0 < curve["next_sequence"] <= (1 << sequence_bits) and
                     0 <= curve["reason"] <= 4 and curve["length"] >= 0,
                     "invalid curve metadata/ID")
            status = curve["unresolved_since"]
            _require((curve["reason"] == 0 and status == -1) or
                     (curve["reason"] != 0 and 0 <= status <= state["time"]),
                     "invalid unresolved status")
            _require(curve["closed"] or curve["winding"] == [0, 0, 0], "open winding")
            _require((3 if curve["closed"] else 2) <= count <= state["node_count"]-total,
                     "invalid curve node count")
            ids.add(curve["id"])
            for _ in range(count):
                f = reader.fields(20, "node")
                node = {key: _integer(f[i+1],bits,version == 5) for i, key in enumerate(
                    ("id", "prev", "next", "parent_left", "parent_right"))}
                node.update(a=_real(f[6]), birth=_real(f[7]),
                            anchor=_integer(f[8], 32, True),
                            placement=_integer(f[9], 32, True),
                            x=[_real(value) for value in f[10:13]],
                            image=[_integer(value, 64, True) for value in f[13:16]],
                            edge_length=_real(f[16]), stretch=_real(f[17]),
                            curvature=_real(f[18], True),
                            curvature_mask=_integer(f[19], 32, True))
                mask, value = node["curvature_mask"], node["curvature"]
                _require(node["edge_length"] >= 0 and node["stretch"] >= 0 and
                         ((mask == 0 and math.isfinite(value) and value >= 0) or
                          (mask in (1, 2) and math.isnan(value))), "invalid diagnostics")
                curve["nodes"].append(node)
            _audit(curve, state["time"], state["domain"], state["dimensions"],sequence_bits,curve_bits,version == 5)
            curve["reconstructed"] = geometry(curve, state["domain"])
            curves.append(curve)
            total += count
        _require(total == state["node_count"], "total count mismatch")
        f = reader.fields(3, "commit", hashed=False)
        _require((_integer(f[1]), _integer(f[2])) == (reader.hash, reader.size),
                 "checksum/commit mismatch")
        _require(stream.read(1) == b"", "trailing data")
    if fluid_path is not None:
        _require(state["kind"] == "checkpoint" and
                 fingerprint(fluid_path) == state["fluid_identity"],
                 "fluid checkpoint identity mismatch")
    state["curves"] = curves
    return state


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path")
    parser.add_argument("--max-nodes", type=int, default=1000000)
    parser.add_argument("--fluid")
    parser.add_argument("--max-mesh-leaves", type=int, default=1000000)
    args = parser.parse_args()
    data = read(args.path, args.max_nodes, args.fluid, args.max_mesh_leaves)
    print(json.dumps({"cycle": data["cycle"], "time": data["time"],
                      "curves": data["curve_count"], "nodes": data["node_count"],
                      "lengths": [c["reconstructed"]["length"] for c in data["curves"]]}))
