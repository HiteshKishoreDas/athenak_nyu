# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this repo is

A fork of **AthenaK** — a block-based AMR astrophysical (M)HD code built on Kokkos for
performance portability (CPU / CUDA / HIP). Upstream is `dfielding14/athenak`
(itself a fork of `IAS-Astrophysics/athenak`); `origin` is `HiteshKishoreDas/athenak_nyu`.

The fork's own work lives mainly on `multiphase_turb` / `development` and targets
**multiphase turbulence**: optically-thin cooling, stochastic turbulence driving,
Lagrangian Monte-Carlo thermodynamic tracers, frame tracking, and new diagnostics.
Fork-added modules are documented in `docs/source/modules/` — read the relevant page
before changing that module.

## Build

Out-of-source only (in-source builds are a hard CMake error). The `kokkos/` submodule
must be initialized (`git submodule update --init --recursive`).

```bash
cmake -S . -B build -D PROBLEM=turb_cond -D Athena_ENABLE_MPI=On
cmake --build build -j 16          # executable: build/src/athena
```

Key options: `PROBLEM` (name of a file in `src/pgen/` without extension; default
`built_in_pgens`), `Athena_ENABLE_MPI`, `Athena_ENABLE_OPENMP`,
`Athena_SINGLE_PRECISION`, `CMAKE_BUILD_TYPE=Debug` (also turns on Kokkos bounds
checking), plus Kokkos flags (`Kokkos_ENABLE_CUDA`, `Kokkos_ARCH_HOPPER90`, …).
CUDA builds need `-D CMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper`.

`PROBLEM` is baked into the build, so **one build directory per problem generator**.
Existing trees here: `build/` (CPU, built-in pgens), `build_gpu/` (CUDA+MPI,
`turb_cond`), `build_cpu_debug/`. `build*/` is gitignored.

### On Rusty (Flatiron)

`rusty_build_script.sh` is the working GPU recipe (modules + cmake + make); adapt it
rather than reinventing module combinations. `rusty_get_gpu.sh` shows the `salloc`
line for an interactive H100 node, and `turb_drive/job_script.sh` is a working
4-GPU sbatch script. `configure_build_turb.sh` is the equivalent for Frontier (HIP/CCE).
Never mix Lmod modules with conda — ABI mismatches at runtime.

## Run

```bash
build/src/athena -i inputs/hydro/blast.athinput -d run_dir
```

`-i` input file, `-r` restart file, `-d` run directory, `-h` for all options.
Any `block/param=value` can be overridden on the command line. Input files are
INI-like `<block>` files; see `inputs/` for canonical examples per physics area and
`tst/inputs/` for regression inputs.

Run output directories in the repo root (`turb_drive*/`) are gitignored scratch, not
source. Large data belongs on `/mnt/ceph`, not in the home filesystem.

## Test and lint

From `tst/`:

```bash
python run_test_suite.py --style     # cpplint + flake8 + custom style checks
python run_test_suite.py --cpu       # pytest suites named *_cpu.py
python run_test_suite.py --mpicpu    # *_mpicpu.py
python run_test_suite.py --gpu "-DKokkos_ARCH_AMPERE80=On -DCMAKE_CXX_COMPILER=.../nvcc_wrapper"
```

The suite builds its own tree at `tst/build/` and writes `tst/test_log.txt`. Tests are
plain pytest files under `tst/test_suite/<area>/`; the name suffix (`_cpu`, `_mpicpu`,
`_gpu`) is what selects them — no registration needed. `tst/run_tests.py` is the older
regression harness (`tst/scripts/`) still referenced by `.gitlab-ci.yml`.

New functionality is expected to come with a regression test (see CONTRIBUTING.md).

## Style

- C++: Google C++ Style Guide, **90-column** limit, `.hpp`/`.cpp`, enforced by
  `tst/scripts/style/check_athena_cpp_style.sh` with the filters in `CPPLINT.cfg`.
  No tab characters — the style check greps for them and fails.
- Python: flake8, 90 columns (`setup.cfg`).
- Files carry the AthenaXXX copyright header and a `//! \file` / `//! \brief` block.

## Architecture

- `Mesh` (`src/mesh/`) owns the AMR tree and the `MeshBlockPack` — physics operates on
  *packs* of MeshBlocks, not single blocks, so kernels stay GPU-friendly.
- `MeshBlockPack` holds the optional physics modules (`phydro`, `pmhd`, `prad`,
  `ppart`, `pz4c`, `pdyngr`, `pturb_driver`, …), each constructed only when its
  `<block>` exists in the input file.
- `Driver` (`src/driver/`) runs named task lists per stage; each module contributes
  tasks via its `*_tasks.cpp`.
- Loops are Kokkos `par_for` over pack/block/k/j/i; anything called inside must be
  `KOKKOS_INLINE_FUNCTION` and device-safe.
- `ParameterInput` (`src/parameter_input.*`) parses the athinput blocks; validate
  parameters at construction time and fail loudly on unknown/legacy keys — the cooling
  module does this deliberately.

### Adding a problem generator

Add `src/pgen/<name>.cpp` defining `ProblemGenerator::UserProblem(ParameterInput*, bool restart)`,
then build with `-D PROBLEM=<name>`. Optional hooks are set by assigning the
corresponding member in `UserProblem`: `user_bcs_func`, `user_srcs_func`,
`user_hist_func`, `user_ref_func`, `user_cooling_func`, … (see `src/pgen/pgen.hpp`;
`src/pgen/turb_cond.cpp` is a good template). Built-in pgens are dispatched from
`src/pgen/pgen.cpp` and must handle `restart == true` by re-registering hooks without
re-initializing the fluid state.

### Fork-added modules worth knowing

| Area | Code | Docs |
|---|---|---|
| Cooling / heating (`<cooling>`) | `src/srcterms/cooling.*`, `ismcooling.hpp` | `docs/source/modules/cooling.md` |
| Turbulence driving (`<turb_driving>`) | `src/srcterms/turb_driver.*` | `docs/source/modules/turbulence_driving.md` |
| Initial perturbations | `src/srcterms/initial_perturbations.*` | `src/srcterms/README.md` |
| Frame tracking | `src/srcterms/frame_tracker.*` | `docs/source/modules/frame_tracking.md` |
| MC thermo tracers | `src/particles/particles_lagrangian_mc.cpp`, `tracer_fields.*` | `docs/source/modules/lagrangian_mc_tracers.md` |
| Outputs (pdf, spherical slice, tracer history) | `src/outputs/pdf.cpp`, `spherical_slice.cpp`, `prtcl_thermo_history.cpp` | — |
| Per-node/rank restart sharding | `src/file_sharding.hpp`, `src/main.cpp` | `docs/file_sharding_explained.md` |

Python readers for the custom outputs live in `vis/python/` (`bin_convert.py`,
`read_pdf.py`, `read_sphslice.py`) and `scripts/`.

## Conventions

- Branch off `development`/`multiphase_turb`; don't commit to `main` (it tracks upstream).
- `.hpp`+`.cpp` pairs; new `.cpp` files must be added to the `add_executable` list in
  `src/CMakeLists.txt` or they simply won't compile in.
- Prefer editing an existing module over adding a parallel one; check
  `docs/source/modules/` for the intended design first.
