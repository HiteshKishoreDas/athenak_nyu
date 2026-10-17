Fresh target material-curve I/O evidence, captured 2026-10-08.

Runtime source commit: 313977256c95e0aee17433afed27379f0b5579bc.
Pinned Kokkos: 08ceff92bcf3a828844480bc1e6137eb74028517 (clean checkout).
Build: Debug, MPI enabled, Serial Kokkos, double precision, bounds checking enabled.

commands.txt records outer commands and working directory. capture.py reproduces
these checks against the recorded source/build paths. results.json records both
suite exits/assertion outcomes. reader.stdout/reader.stderr and native.stdout/
native.stderr preserve suite output. native-subprocesses.jsonl records each exact
compiler/MPI command, cwd, exit, stdout/stderr, and duration, including expected
collective malformed-input failures.

Fresh results: 10 independent reader tests pass; native persistence tests pass.
The native log has 26 subprocess calls: compile, save, and three rank-changing
loads succeed (five exit0 calls); 20 malformed sidecars and one mixed fluid pair
fail as expected (21 exit1 calls), with material-curve error text verified.
Native save uses two ranks; loads use four, one, and two ranks, including ranks
with no resident curves. Exact IDs, metadata, and AMR clocks are compared.

The native harness links the production I/O implementation but stubs mesh
constructors, sampling, and structural Audit. It tests persistence preflight and
distribution; it does not establish evolving solver/checkpoint-hook correctness.
That is covered by the separate production validation worker.

provenance.json records source HEAD/status, SHA256 of tested sources/headers,
actual runtime athena executable, matching Kokkos static library, launcher,
compile_commands, and config.hpp. source/ preserves exact tested I/O inputs;
build/ preserves generated compiler flags, full compile_commands, config.hpp,
CMakeCache and linker command. toolchain.json records compiler/MPI launcher
versions. restart-port.diff preserves the native hook/guards diff against the
multiphase base526a187f. Existing fluid metadata/counters/native physics payload
code remains unchanged when curves are disabled.

The raw build/git-state.txt is 08ceff92b-CLEAN (configure-time build metadata),
which differs from current runtime source HEAD and matches the Kokkos prefix.
It is preserved verbatim; runtime/source provenance uses HEAD plus exact hashes,
without claiming this raw embedded configure value identifies the port commit.

Schema2 is retained. Fluid checksum binds the actual native time/dt/cycle
checkpoint metadata. dtold/dt_last_completed are not curve state and are used by
physics excluded from this bounded port. AMR clocks/leaf identity/configuration
remain mandatory and preflight before restoration.

No old-branch test outputs were reused and no repository sources were edited for
this evidence capture. The capture helper and archives live only under /private/tmp.
