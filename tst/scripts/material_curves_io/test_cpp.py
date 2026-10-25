"""Build isolated C++ persistence checks using configured AthenaK compiler flags.

Example (after configuring with CMAKE_EXPORT_COMPILE_COMMANDS=ON):
  python3 test_cpp.py --build BUILD --kokkos-library BUILD/kokkos/core/src/libkokkoscore.a
MPI configurations additionally use --launcher mpiexec (or installed prterun).
This harness stubs mesh constructors, sampling and structural Audit; it verifies
production serialization/parser/distribution, while production.py tests evolution.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("schema_tests", Path(__file__).with_name(
    "test_reader.py"))
SCHEMA = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SCHEMA)
READER = SCHEMA.reader


def run(command, expected_failure=False):
    result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    if (result.returncode == 0) == expected_failure:
        raise AssertionError(f"unexpected exit {result.returncode}: {command}\n"
                             + result.stdout + result.stderr)
    if expected_failure and "material curves:" not in result.stderr:
        raise AssertionError("rejection did not identify material curves: "
                             + result.stderr)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--kokkos-library", type=Path, required=True)
    parser.add_argument("--launcher", default="mpiexec")
    args = parser.parse_args()
    build = args.build.resolve()
    item = json.loads((build / "compile_commands.json").read_text())[0]
    command = shlex.split(item["command"])
    index = command.index("-o")
    del command[index:index+2]
    command.remove("-c")
    command.pop()
    mpi = "#define MPI_PARALLEL_ENABLED 1" in (build / "config.hpp").read_text()
    if mpi:
        command[0] = shutil.which("mpicxx") or "mpicxx"
    with tempfile.TemporaryDirectory(prefix="material-io-check-") as directory:
        directory = Path(directory)
        executable = directory / "roundtrip"
        command += [str(ROOT / "tst/scripts/material_curves_io/roundtrip.cpp"),
                    str(ROOT / "src/material_curves/material_curves_io.cpp"),
                    str(args.kokkos_library.resolve()), "-o", str(executable)]
        result = subprocess.run(command, cwd=item["directory"], capture_output=True,
                                text=True, timeout=60)
        if result.returncode:
            raise AssertionError(result.stdout + result.stderr)
        fluid = directory / "fluid.rst"

        def launch(mode, ranks=1, failure=False):
            prefix = [args.launcher, "-n", str(ranks)] if mpi else []
            return run(prefix + [str(executable), mode, str(fluid)], failure)

        launch("save", 2 if mpi else 1)
        checkpoint = fluid.with_suffix(".rst.curves")
        original = checkpoint.read_bytes()
        saved = READER.read(checkpoint, fluid_path=fluid)
        for ranks in ((4, 1, 2) if mpi else (1,)):
            launch("load", ranks)
            restored = READER.read(str(fluid) + ".restored")
            assert restored["amr"] == saved["amr"], "AMR clocks/config changed"
            saved_by_id = {c["id"]: c for c in saved["curves"]}
            restored_by_id = {c["id"]: c for c in restored["curves"]}
            # Dict equality excludes NaN-bearing diagnostic columns.
            for cid, curve in saved_by_id.items():
                other = restored_by_id[cid]
                for key in ("id", "closed", "winding", "period", "reason",
                            "unresolved_since", "next_sequence", "inserted",
                            "migrations", "samples"):
                    assert curve[key] == other[key], (cid, key)
                for left, right in zip(curve["nodes"], other["nodes"]):
                    for key in ("id", "prev", "next", "parent_left", "parent_right",
                                "a", "birth", "anchor", "placement", "x"):
                        assert left[key] == right[key], (cid, key)
        body = original[:original.rfind(b"commit ")].decode("ascii")
        def alter_record(key, index, value, predicate=lambda fields: True):
            rows = body.splitlines()
            for i, row in enumerate(rows):
                fields = row.split()
                if fields[0] == key and predicate(fields):
                    edits = (zip(index, value) if isinstance(index, tuple)
                             else [(index, value)])
                    for column, replacement in edits:
                        fields[column] = str(replacement)
                    rows[i] = " ".join(fields)
                    break
            return SCHEMA.complete("\n".join(rows) + "\n")

        malformed = [original[:-1], original[:-30], original + b"extra\n",
                     SCHEMA.complete(body.replace("MATERIAL_CURVES 4",
                                                  "MATERIAL_CURVES 1")),
                     SCHEMA.complete(body.replace("ASCII 64", "ASCII 32")),
                     alter_record("config", 3, 2),
                     alter_record("state", 2, 3),
                     SCHEMA.complete(body.replace(" 1.5 1 4 ", " -1 1 4 ")),
                     alter_record("curve", 5, 1, lambda f: f[2] == "1"),
                     alter_record("curve", 3, 2**40+1, lambda f: f[2] == "1"),
                     alter_record("node", 12, 0.25, lambda f: f[8] == "0"),
                     alter_record("node", 10, 2**40+1),
                     alter_record("node", (4, 5), ((1 << 24) | 2, (1 << 24) | 1),
                                  lambda f: f[8] == "0" and int(f[1]) >> 24 == 1),
                     alter_record("node", 5, (255 << 24) | 2,
                                  lambda f: f[8] == "0" and int(f[1]) >> 24 != 1),
                     alter_record("amr", 5, 1), alter_record("amr", 3, 6),
                     alter_record("leaf", 2, 10), alter_record("leaf", 5, -1),
                     alter_record("leaf", 5, 18),
                     SCHEMA.complete(body.replace("amr 1 3 5 1 2\n", ""))]
        for corrupted in malformed:
            assert corrupted != original, "ineffective fault injection"
            checkpoint.write_bytes(corrupted)
            launch("load", 2 if mpi else 1, failure=True)
        checkpoint.write_bytes(original)
        fluid.write_bytes(fluid.read_bytes() + b"mixed fluid state")
        launch("load", 2 if mpi else 1, failure=True)
        print("PASS C++ save/load, exact IDs and metadata, malformed/pair rejection"
              + (", MPI 2 -> 4/1/2 ranks with empty ranks" if mpi else ", serial"))


if __name__ == "__main__":
    main()
