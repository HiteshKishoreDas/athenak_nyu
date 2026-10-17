"""Independent schema/geometry/corruption checks; no AthenaK runtime required."""
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest

MODULE = Path(__file__).resolve().parents[3] / "vis/python/material_curves.py"
spec = importlib.util.spec_from_file_location("material_curves", MODULE)
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)


def complete(body):
    checksum = reader.OFFSET
    data = body.encode("ascii")
    for byte in data:
        checksum = ((checksum ^ byte) * reader.PRIME) & reader.MASK64
    return data + f"commit {checksum} {len(data)}\n".encode()


def fixture(closed=False, shuffled=False, empty=False):
    count = 0 if empty else 1
    cid = (1 << 32)-1  # IDs near UINT64_MAX, beyond binary64's integer range.
    n = 4 if closed else 3
    total = 0 if empty else n
    body = (f"MATERIAL_CURVES 2 ASCII 64 U64 I64 I32\n"
            f"state snapshot 2 17 {count} {total}\n"
            "domain 0 1 0 1 0 1 1 1 1 1 1 1 8 8 8\n"
            "config 100 0.1 1 0.2 0 1 1\nbasename \"curves test\"\nfluid 0 0\n"
            "amr 0 0 0 0 0\n"
            f"{reader.CURVE_COLUMNS}\n{reader.NODE_COLUMNS}\n")
    if empty:
        return body
    ids = [(cid << 32) | (i+1) for i in range(n)]
    if not closed:
        ids[1], ids[2] = ids[2], ids[1]
    if closed:
        points = [(0, 0), (1, 0), (1, 1), (0, 1)]
        # Material interval is 4, one positively wound seam edge ends at (1,0).
        winding, period, inserted = "1 0 0", 4, 0
    else:
        points = [(0, 0), (2, 0), (6, 0)]
        winding, period, inserted = "0 0 0", 3, 1
    body += (f"curve {cid} {int(closed)} {winding} {period} -1 0 {n+1} "
             f"{inserted} 9 20 {n} 6\n")
    rows = []
    for i, (x, y) in enumerate(points):
        previous = ids[i-1] if i or closed else 0
        following = ids[(i+1) % n] if i+1 < n or closed else 0
        born = not closed and i == 1
        parents = f"{ids[0]} {ids[-1]}" if born else "0 0"
        a = i if closed or i < 2 else 3
        curvature, mask = ("0", "0") if i == 1 and not closed else ("nan", "1")
        if closed:
            curvature, mask = "1.4142135623730951", "0"
        rows.append(f"node {ids[i]} {previous} {following} {parents} {a} "
                    f"{1 if born else 0} {0 if born else 1} {2 if born else 0} "
                    f"{x} {y} 0 {x} {y} 0 0 0 {curvature} {mask}\n")
    if shuffled:
        rows.reverse()
    return body + "".join(rows)


class ReaderTests(unittest.TestCase):
    def parse(self, data, **kwargs):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.curves"
            path.write_bytes(data)
            return reader.read(path, **kwargs)

    def test_open_birth_large_ids_and_geometry(self):
        data = self.parse(complete(fixture()))
        curve = data["curves"][0]
        self.assertEqual(curve["id"], (1 << 32)-1)
        self.assertEqual(curve["nodes"][0]["id"], ((1 << 32)-1) << 32 | 1)
        self.assertEqual(curve["nodes"][1]["birth"], 1)
        self.assertEqual(curve["nodes"][1]["placement"], 2)
        self.assertEqual(curve["reconstructed"]["length"], 6)
        self.assertEqual(curve["reconstructed"]["stretch"], [2, 2, 0])
        self.assertEqual(curve["reconstructed"]["curvature"][1], 0)
        self.assertEqual(curve["reconstructed"]["curvature_mask"], [1, 0, 1])

    def test_winding_and_shuffled_rows(self):
        original = self.parse(complete(fixture(closed=True)))
        shuffled = self.parse(complete(fixture(closed=True, shuffled=True)))
        a, b = original["curves"][0], shuffled["curves"][0]
        self.assertEqual(a["nodes"], b["nodes"])
        self.assertEqual(a["winding"], [1, 0, 0])
        self.assertEqual(a["reconstructed"]["edges"][-1], [1, -1, 0])
        self.assertAlmostEqual(a["reconstructed"]["length"], 3 + math.sqrt(2))

    def test_analytic_seed_provenance(self):
        body = fixture().replace(" 1 0 2 2 0 0", " 1 0 3 2 0 0")
        self.assertNotEqual(body, fixture())
        data = self.parse(complete(body))
        self.assertEqual(data["curves"][0]["nodes"][1]["placement"], 3)

    def test_closure_parent_interval(self):
        rows = fixture(closed=True).splitlines()
        cid = (1 << 32)-1
        for i, row in enumerate(rows):
            fields = row.split()
            if fields[0] == "curve":
                fields[10] = "1"  # one born sample, three original anchors
            elif fields[0] == "node" and int(fields[1]) == (cid << 32) | 4:
                fields[4:6] = [str((cid << 32) | 3), str((cid << 32) | 1)]
                fields[7:10] = ["1", "0", "2"]
            rows[i] = " ".join(fields) if fields[0] in ("curve", "node") else row
        data = self.parse(complete("\n".join(rows) + "\n"))
        born = data["curves"][0]["nodes"][-1]
        self.assertEqual(born["parent_right"], (cid << 32) | 1)
        self.assertEqual(born["a"], 3)

    def test_preflight_geometry_and_parent_rejections(self):
        body = fixture().replace("8 8 8\n", "8 8 1\n")
        self.parse(complete(body))

        def alter(key, index, value, predicate=lambda fields: True, source=body):
            rows = source.splitlines()
            for i, row in enumerate(rows):
                fields = row.split()
                if fields[0] == key and predicate(fields):
                    edits = (zip(index, value) if isinstance(index, tuple)
                             else [(index, value)])
                    for column, replacement in edits:
                        fields[column] = str(replacement)
                    rows[i] = " ".join(fields)
                    break
            changed = "\n".join(rows) + "\n"
            self.assertNotEqual(changed, source)
            return complete(changed)

        cid = (1 << 32)-1
        closed = fixture(closed=True).replace("8 8 8\n", "8 8 1\n")
        malformed = [alter("node", 12, 0.25, lambda f: f[8] == "0"),
                     alter("node", 10, 2**40+1),
                     alter("curve", 3, -(2**40+1), source=closed),
                     alter("curve", 5, 1, source=closed),
                     alter("node", (4, 5), ((cid << 32) | 2, (cid << 32) | 1),
                           lambda f: f[8] == "0")]
        for data in malformed:
            with self.assertRaises(ValueError):
                self.parse(data)

    def test_amr_ledger_and_rejection(self):
        inactive = "amr 0 0 0 0 0\n"
        active = "amr 1 3 5 1 2\nleaf 1 0 0 0 2\nleaf 1 1 0 0 4\n"
        body = fixture().replace(inactive, active)
        amr = self.parse(complete(body))["amr"]
        self.assertEqual(amr["refinement_interval"], 5)
        self.assertEqual([leaf["cycles_since_ref"] for leaf in amr["leaves"]], [2, 4])
        malformed = [body.replace("amr 1 3 5 1 2\n", ""),
                     body.replace("amr 1 3 5 1 2", "amr 1 3 5 1 3"),
                     body.replace("amr 1 3 5 1 2", "amr 1 0 5 1 2"),
                     body.replace("leaf 1 1 0 0 4", "leaf 1 1 0 0 -1"),
                     body.replace("leaf 1 1 0 0 4", "leaf 1 1 0 0 18"),
                     body.replace("leaf 1 1 0 0 4", "leaf 1 0 0 0 4"),
                     fixture().replace(inactive, "")]
        for changed in malformed:
            with self.assertRaises(ValueError):
                self.parse(complete(changed))
        with self.assertRaises(ValueError):
            self.parse(complete(body), max_mesh_leaves=1)

    def test_empty(self):
        self.assertEqual(self.parse(complete(fixture(empty=True)))["curves"], [])

    def test_circle_estimator(self):
        nodes = [{"a": i*2*math.pi/32, "x": [math.cos(i*2*math.pi/32),
                  math.sin(i*2*math.pi/32), 0]} for i in range(32)]
        curve = {"nodes": nodes, "closed": True, "winding": [0, 0, 0],
                 "period": 2*math.pi}
        geometry = reader.geometry(curve, [0, 10, 0, 10, 0, 10])
        self.assertAlmostEqual(geometry["length"], 64*math.sin(math.pi/32))
        self.assertTrue(all(abs(k-1) < 1e-13 for k in geometry["curvature"]))

    def test_rejections(self):
        body = fixture()
        cid = (1 << 32)-1
        first = (cid << 32) | 1
        malformed = [
            complete(body)[:-1],  # missing commit newline
            complete(body)[:-30],
            complete(body) + b"extra\n",
            complete(body.replace("ASCII 64", "BINARY 64")),
            complete(body.replace("ASCII 64", "ASCII 128")),
            complete(body.replace("MATERIAL_CURVES 2", "MATERIAL_CURVES 1")),
            complete(body.replace("state snapshot 2 17 1 3", "state snapshot 2 17 1 4")),
            complete(body.replace(f"node {first} 0", f"node {first} 1")),
            complete(body.replace(" -1 0 4 ", " -1 4 4 ")),
            complete(body.replace(" 1 0 2 2 0 0", " 3 0 2 2 0 0")),  # future birth
            complete(body.replace(f"node {first} ", f"node {1 << 64} ")),
            complete(body.replace(" 2 0 0 2 0 0", " nan 0 0 2 0 0")),
            complete(body).replace(b"curves test", b"curves best"),
            complete(body.replace("nan 1", "inf 1")),
            complete(body.replace("nan 1", "nan 0")),
            b"x"*8193,
        ]
        for data in malformed:
            with self.subTest(data=data[:40]):
                with self.assertRaises(ValueError):
                    self.parse(data)
        with self.assertRaises(ValueError):
            self.parse(complete(body), max_nodes=2)

    def test_fluid_fingerprint(self):
        with tempfile.TemporaryDirectory() as directory:
            fluid = Path(directory) / "fluid.rst"
            fluid.write_bytes(b"a fluid state\x00")
            size, checksum = reader.fingerprint(fluid)
            body = fixture().replace("state snapshot", "state checkpoint")
            body = body.replace("fluid 0 0", f"fluid {size} {checksum}")
            data = complete(body)
            self.parse(data, fluid_path=fluid)
            fluid.write_bytes(b"a fluid state\x01")
            with self.assertRaises(ValueError):
                self.parse(data, fluid_path=fluid)


if __name__ == "__main__":
    unittest.main()
