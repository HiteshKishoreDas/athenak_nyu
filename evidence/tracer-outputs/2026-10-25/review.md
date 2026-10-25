# Independent output-integration reviews

Two existing agents reviewed read-only under the standing numerical timing and
persistence gates. Both applied ponytail-review. Production edits remained with
the coordinator; no overlapping agent builds or edits occurred.

## Persistence/simplicity

`/root/review_curve_id_persistence` checked literal THP v2 byte ordering, signed
integer slots, empty population metadata, immutable curve chunk replay, historical
CSV compatibility, unavailable NaN estimates and partial field validation.

Found and fixed:

1. Empty blocks skipped finite/nonnegative header clock validation. Validate cycle
   and time immediately after block-header unpacking.
2. A same-cycle empty block could overwrite the time of real observations. Reject
   conflicting same-cycle clocks before rows, including empty blocks.
3. Partial field selections could accept an unavailable estimate marked valid.
   Every selected valid=1 flag must satisfy each available prerequisite; full
   equivalence is checked when all prerequisites are selected. A finite estimate
   with valid=0 stays allowed when parent validity was not selected.

One immutable curve file now requires exactly one normal THP block. Generic
particle streams retain legitimate multi-block/identical terminal replay, including
NaNs restricted to the two named shock Mach estimates. Rank-completeness checking
uses observed rank-set size, avoiding allocation proportional to an untrusted
claimed rank count. Six representative malformed binary inputs and expected
rejections are retained under reader-regressions/; the complete 21-case generator
is tst/tracer-output/thp_reader_check.py.

Reviewer independently reran all 21 reader cases and git diff --check: passed.
Final conclusion: no remaining persistence blocker; no useful complexity cuts.

## Numerical timing/simplicity

`/root/review_curve_thermo_timing` checked post-AMR synchronization, containing-cell
sampling, cached versus curve-only mesh detection, unreduced membership, strength
validity, passivity and birth/history semantics.

Found and fixed: zero-step restarts restore encounter generation/time but not the
derived mesh detector cache. Same-generation ObserveShocks returns deliberately,
so generation/time checks alone admitted uninitialized signals. The unchanged
cache detector kernel is now RefreshShockSignals, called after restart primitive
initialization separately from encounter updates. It changes no counters,
observations or journals. SHOUT04 compares particle/curve properties and all six
mesh arrays before/after a zero-step active-AMR restart at cycle13, then checks the
advancing continuation. Both pass exactly.

Also fixed the passivity harness: dcycle=-1 still constructs initial/final native
outputs. The mesh-off control uses dcycle=0 and dt=-1, which disables construction.
Controls also omit optional traces. All fluid arrays, curve geometry and complete
histories remain unchanged in planar, live-AMR and short accepted forced-box runs.

Reviewer confirmed the cache fix and regression read-only, without builds/native
reruns. Final conclusion: no remaining timing blocker. The coordinator ran final
source-bound native suites after the fixes and after including native mesh-output
code in the existing build fingerprint.

## Remaining review boundary

Native MPI compilation passes; multi-rank runtime remains unqualified because the
prior local PMIx/PRRTE mismatch has not been repaired/retested. GPU/cluster costs,
production budgets and human scientific approval remain pending. These automated
and code-review results do not replace scientific sign-off.
