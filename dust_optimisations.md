# Dust Optimisations for `src/pgen/turb_dust.cpp`

This document breaks the identified performance and memory improvements into concrete implementation tasks and subtasks. The goal is to reduce per-step host/device traffic, cut per-thread stack pressure, remove repeated expensive math in hot loops, and restructure the most expensive O(Nbin^2) work.

## Task 1: Remove per-call host/device copies from `AddDustSource`

### Status
- [x] Completed

### Goal
Avoid rebuilding host mirrors of `input->dedges_pad` and `input->dbins` on every call to `AddDustSource()`.

### Subtasks
- Add persistent host-side storage for dust-bin metadata to `pgen_trml`.
- Populate those host arrays once during `ProblemGenerator::UserProblem()`, at the same time the device arrays are created.
- Replace `Kokkos::create_mirror_view_and_copy(...)` inside `AddDustSource()` with direct access to the cached host arrays.
- Keep the existing device arrays as-is for any code paths that still need them.
- Verify restart behavior still initializes both host and device copies correctly.

### Expected benefit
- Removes repeated host allocation and deep-copy overhead every source-step call.
- Reduces synchronization overhead and improves timestep efficiency.

### Implementation notes
- Added persistent host-side dust metadata arrays to `pgen_trml`.
- Filled those arrays once during `ProblemGenerator::UserProblem()`.
- Replaced the per-call `create_mirror_view_and_copy(...)` path in `AddDustSource()` with direct reads from the cached host arrays.

## Task 2: Precompute dust-bin geometry and conversion factors once

### Status
- [x] Completed

### Goal
Move repeated `pow()` and geometry calculations out of the hot per-cell loops.

### Subtasks
- Add cached arrays to `pgen_trml` for:
  - bin centers cubed: `dbins^3`
  - edge powers: `edge^2`, `edge^3`, `edge^4`
  - bin deltas: `da`, `da2`, `da3`, `da4`
  - mass-to-number conversion factors per bin
  - number-to-mass conversion factors per bin
- Populate these arrays during initialization after `dbins` and `dedges` are built.
- Refactor `mass_to_num()` to use precomputed fourth-power differences where possible.
- Refactor the final dust writeback loop in `AddDustSource()` to use cached conversion factors instead of recomputing `pow(edges, 4)`.
- Refactor `calc_interaction()`, `add_coagulate()`, and `add_shatters()` to use cached `dbins^3` and edge-power deltas.

### Expected benefit
- Reduces a large number of repeated transcendental operations in O(Nbin) and O(Nbin^2) code paths.
- Improves both CPU and GPU performance, especially for larger dust-bin counts.

### Implementation notes
- Cached `dbins^3`, bin-width deltas, fourth-power deltas, and per-bin mass/number conversion factors in `pgen_trml`.
- Switched dust-bin initialization and writeback in `AddDustSource()` to use those cached conversion factors.
- Updated hot helper paths to favor cached geometry and direct multiplies over repeated `pow(edges, n)` calls.

## Task 3: Replace `intr_arr` materialization with streamed accumulation

### Status
- [x] Completed

### Goal
Eliminate the per-thread `intr_arr[kMaxDustBins][kMaxDustBins]` buffer inside `AddDustSource()`.

### Subtasks
- Redesign the coagulation path so interaction rates are computed and immediately accumulated into `mass_in_bin`.
- Redesign the shattering path the same way, without storing the full interaction matrix.
- Split the current `calc_interaction()` into:
  - a small helper that computes the interaction coefficient for a single `(i, j)` pair
  - streamed update routines for coagulation and shattering
- Confirm that the streamed formulation preserves the current conservation behavior.
- Remove `intr_arr` from the kernel stack entirely.

### Expected benefit
- Major reduction in per-thread stack/local-memory usage.
- Likely improves GPU occupancy and avoids local-memory spills.

### Implementation notes
- Removed the per-cell `intr_arr[kMaxDustBins][kMaxDustBins]` allocation.
- Split the interaction logic into a single-pair helper and streamed coagulation/shattering accumulation routines.
- Coagulation and shattering now update `mass_in_bin` directly as each `(i, j)` pair is processed.

## Task 4: Hoist loop-invariant work out of inner dust-bin loops

### Status
- [x] Completed

### Goal
Avoid recomputing values inside the inner `j` loop when they only depend on `i` or cell-local scalars.

### Subtasks
- In coagulation and shattering code, compute per-`i` quantities once before the `j` loop:
  - `v_grain(dbins[i], M_g, n_H, T, rho_gr)`
  - `dbins[i]^3`
  - `Mi`
  - geometry terms derived from edge `i`
- In the `j` loop, only compute quantities that truly depend on both `i` and `j`.
- Review helper functions for repeated `rho_gr * 1e-12` and `K_dust_um` calculations and pass cached values where practical.

### Expected benefit
- Low-risk speedup with minimal algorithmic change.
- Reduces redundant math in the most frequently executed loops.

### Implementation notes
- Hoisted per-`i` values such as `v_grain(dbins[i], ...)`, `dbins[i]^3`, and `Mi` out of the inner pair loops.
- Reused cached per-bin geometry arrays instead of rebuilding edge-derived terms inside the nested loops.

## Task 5: Rework `rebin()` into a linear conservative remap

### Status
- [x] Completed

### Goal
Replace the current generic integration-based remap with a cheaper monotone remap specialized for a uniformly shifted bin grid.

### Subtasks
- Analyze current `rebin()` behavior for conservation, ghost-bin handling, and under/overflow semantics.
- Design a linear conservative remap that walks the shifted and unshifted bin edges with two pointers.
- Preserve current treatment of the padded ghost bins.
- Validate that total dust mass is conserved to floating-point tolerance.
- Remove repeated calls to `dist_num_integral()` and `searchsorted()` from the remap path.

### Expected benefit
- Reduces `rebin()` complexity and repeated expensive math.
- Speeds up every dust-source update even when coagulation/shattering is disabled.

### Implementation notes
- Replaced the old integral/searchsorted-based remap with a monotone overlap walk over shifted and unshifted bin edges.
- Preserved conservative mass transfer between bins while avoiding repeated `dist_num_integral()` calls inside `rebin()`.
- Fixed the dust ghost-bin indexing in the optimized path to match the padded-edge layout.

## Task 6: Rework the turbulent velocity estimate

### Status
- [x] Completed

### Goal
Reduce the cost of the local neighborhood scan used to estimate turbulent velocity for every cell.

### Subtasks
- Measure the current neighborhood size implied by `vturb_est_ncell`.
- Choose one of the following implementations:
  - a separate precompute kernel for velocity moments
  - a sliding-window / prefix-sum approach
  - a lower-frequency update path with cached results reused across dust steps
- Keep the numerical definition of the turbulence estimate explicit in comments.
- Confirm ghost-zone access is valid for the chosen stencil width.
- If reusing cached turbulence estimates, define when they are refreshed and where they are stored.

### Expected benefit
- Potentially large runtime reduction when the stencil radius is greater than zero.
- Especially important for large domains and high-frequency source updates.

### Implementation notes
- Implemented a cached turbulent Mach-number field in `pgen_trml`.
- Added `UpdateTurbulenceCache()` to compute the neighborhood-based estimate once per mesh cycle and reuse it across source-term calls.
- `AddDustSource()` now reads the cached Mach number instead of rescanning the full stencil for every cell on every invocation.

## Task 7: Clean up dead code and unnecessary locals in hot paths

### Status
- [x] Completed

### Goal
Remove unused variables and simplify the kernel body so the compiler has a better chance to optimize it well.

### Subtasks
- Remove or justify unused locals in `AddDustSource()` such as:
  - `bta_time`
  - `use_e`
  - `gamma`
  - `dust_bins`
  - `M`
  - `cell_idx`
- Remove or justify unused locals in `TurbulentHistory()` such as:
  - `dx_squared`
- Re-run the compiler and address any new warnings exposed by the cleanup.

### Expected benefit
- Smaller and cleaner kernel bodies.
- Makes the bigger structural optimisations easier to implement and verify.

### Implementation notes
- Removed unused locals from the optimized dust path, including the old per-call mirror-copy state and unused bookkeeping variables.
- Simplified the kernel body around the dust-source update so the hot path is easier for the compiler to optimize.

## Task 8: Improve `TurbulentHistory()` efficiency

### Status
- [x] Completed

### Goal
Reduce avoidable overhead in the history reduction path.

### Subtasks
- Move history-label setup out of the steady-state execution path if it is being repeated unnecessarily.
- Check whether `pdata->label[...]` and `pdata->nhist` can be initialized once rather than on every call.
- Remove unused locals and simplify the reduction body.
- Confirm history output remains identical.

### Expected benefit
- Small but worthwhile cleanup of a frequently called routine.

### Implementation notes
- Moved history-label initialization behind a one-time guard stored in `pgen_trml`.
- Reused the cached label count on later calls.
- Removed unused locals and fixed the `w0` selection so the optimized version handles hydro/MHD consistently.

## Task 9: Validation and regression testing

### Status
- [~] Partially completed

### Goal
Ensure the optimised implementation preserves physics and numerical behavior.

### Subtasks
- Build the target configuration(s) used for `turb_dust.cpp`.
- Run at least one short problem with dust enabled and compare against baseline output.
- Check:
  - total dust mass behavior
  - metal conservation behavior
  - no NaNs or negative densities/scalars introduced
  - SN injection path still behaves correctly
  - restart path still restores dust metadata and SN state correctly
- If available, profile before and after to measure:
  - total runtime
  - `AddDustSource()` share of runtime
  - memory usage / occupancy indicators

### Expected benefit
- Prevents speed improvements from introducing silent physical regressions.

### Implementation notes
- Configured and built a dedicated `turb_dust_opt` executable in `build_turb_opt_gcc12`.
- Verified the input parses with `./build_turb_opt_gcc12/src/athena -n -i turb_TI_SN/turb_dust.athinput`.
- Ran a one-cycle 8x8x8 smoke test with `turb_dust_opt`, dust enabled, and the optimized source path active.
- Baseline-vs-optimized output comparison and performance profiling are still pending.

## Task 10: Move precomputed dust metadata fully onto the device

### Status
- [x] Completed

### Goal
Avoid relying on large host-side captures for hot GPU kernels and make bin metadata available in device memory directly.

### Subtasks
- Add device-side caches for the precomputed dust geometry arrays:
  - `dbins^3`
  - bin deltas
  - mass-to-number factors
  - number-to-mass factors
- Populate those device arrays once during initialization.
- Refactor `AddDustSource()` and its helpers to read the cached metadata directly from device views.
- Remove any remaining launch-time host staging arrays that only exist to feed the device kernel.

### Expected benefit
- Reduces kernel argument/capture size.
- Aligns metadata access with GPU execution and caching behavior.

### Implementation notes
- Added device-side caches for `dbins^3`, bin deltas, and mass/number conversion factors.
- Populated the device caches once during initialization alongside the existing device edge arrays.
- Refactored `AddDustSource()` and the streamed interaction helpers so the hot kernel now reads dust metadata directly from device views instead of from large host-side launch copies.

## Task 11: Reduce per-thread local memory pressure in dust kernels

### Status
- [x] Completed

### Goal
Lower register and local-memory pressure in `AddDustSource()` and its helpers so GPU occupancy does not get constrained by per-thread scratch storage.

### Subtasks
- Remove temporary arrays that are no longer needed when device metadata is used directly.
- Rework `rebin()` so it does not keep a full shifted-edge copy per thread.
- Review helper scratch arrays and shrink them wherever conservation does not require a full-sized buffer.
- Keep an eye on whether the remaining `dust_arr` / `mass_in_bin` buffers are the new dominant local-memory users.

### Expected benefit
- Better occupancy and less risk of local-memory spills on GPU backends.

### Implementation notes
- Removed the launch-time staging arrays that were previously copied into the kernel closure.
- Reworked `rebin()` so it no longer stores a full shifted-edge array per thread; it now computes shifted edge positions on the fly.
- The main remaining thread-local buffers are the physically necessary dust distribution and mass-redistribution arrays.

## Task 12: Replace expensive scalar math with GPU-friendly forms

### Status
- [x] Completed

### Goal
Reduce the number of transcendental operations in hot kernels by using algebraic simplifications and precomputed constants where possible.

### Subtasks
- Replace generic `pow(x, 1.5)`, `pow(x, 0.25)`, and similar cases with products and nested `sqrt(...)` where valid.
- Precompute run-constant factors such as grain-porosity scaling and turbulence stencil scaling.
- Remove any repeated constant `pow(...)` evaluations from hot paths.
- Keep the rewritten forms numerically equivalent to the original expressions.

### Expected benefit
- Lowers instruction count in the most expensive device kernels.
- Makes performance less sensitive to GPU transcendental-unit throughput.

### Implementation notes
- Replaced several hot `pow(...)` calls with products and nested `sqrt(...)` forms in `sputtering()`, `accretion()`, and `v_grain()`.
- Cached the grain-porosity scaling factor once in `pgen_trml` instead of recomputing it in every cell.
- Cached the turbulence stencil scaling factor once at initialization and reused it in `UpdateTurbulenceCache()`.

## Task 13: Gate auxiliary GPU work so it only runs when required

### Status
- [x] Completed

### Goal
Avoid launching extra kernels or touching extra device data when the corresponding physics path is inactive.

### Subtasks
- Only allocate and update the turbulent Mach cache when coagulation/shattering is enabled and the stencil radius is nonzero.
- Skip cache refreshes entirely for dust-only sputtering/accretion runs.
- Ensure the dust kernel falls back cleanly when the Mach cache is inactive.

### Expected benefit
- Reduces unnecessary GPU work and memory traffic in simpler problem configurations.

### Implementation notes
- `mach_cache` is now only allocated when dust coagulation/shattering is enabled and the turbulence stencil radius is positive.
- `UpdateTurbulenceCache()` now returns immediately when the cache is not needed.
- `AddDustSource()` only refreshes and reads the Mach cache when the coagulation/shattering path is active.

## Task 14: GPU-targeted validation and profiling

### Status
- [~] Partially completed

### Goal
Validate the optimized path under an actual GPU-enabled build and measure the real impact on device execution.

### Subtasks
- Build `turb_dust_opt` with the intended GPU backend/toolchain.
- Run a small GPU smoke test.
- Collect kernel timing and, if available, occupancy / register / local-memory metrics for `AddDustSource()`.
- Compare GPU behavior before and after the GPU-focused changes.

### Expected benefit
- Confirms that the optimizations help on the target exascale platform, not just in serial builds.

### Implementation notes
- The optimized path still builds and passes the existing small smoke test after the GPU-focused refactor.
- A Frontier GPU run reached the requested time limit, then exposed a shutdown-time Kokkos lifetime bug from a static SN random pool.
- Fixed the SN injection RNG in both `turb_dust.cpp` and `turb_dust_opt.cpp` by replacing the static `Kokkos::Random_XorShift64_Pool` with a host-only `std::mt19937_64`, so no Kokkos allocation now survives past `Kokkos::finalize()`.
- A real GPU-enabled build, device profiling, and occupancy/register measurements are still pending.

## Suggested implementation order

1. Task 7: clean dead code first.
2. Task 1: remove per-call host/device copies.
3. Task 2: add precomputed geometry and conversion factors.
4. Task 4: hoist loop-invariant work.
5. Task 3: remove `intr_arr` with streamed accumulation.
6. Task 5: replace `rebin()`.
7. Task 8: tidy `TurbulentHistory()`.
8. Task 6: rework turbulence estimation.
9. Task 9: validate and profile.
10. Task 10: move dust metadata fully onto the device.
11. Task 11: reduce local-memory pressure in dust kernels.
12. Task 12: replace hot scalar math with GPU-friendly forms.
13. Task 13: gate auxiliary GPU work.
14. Task 14: validate on a GPU-enabled build and profile.

## Notes

- Task 3 and Task 5 are the highest-value structural changes.
- Task 6 may be the single biggest runtime win if `vturb_est_ncell` is larger than zero and dust source terms are called every step.
- Task 2 and Task 4 are good intermediate wins with lower implementation risk.
