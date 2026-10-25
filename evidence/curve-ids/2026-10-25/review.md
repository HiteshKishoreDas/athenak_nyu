# Independent signed-ID persistence/simplicity review

Reviewer: existing bounded review agent `/root/review_curve_id_persistence`.
Read-only review against `ecd36eff`, following the user's standing independent
persistence-review gate and ponytail-review. No reviewer file edits or broad
native reruns.

## Finding and fix

SaveShockRestart retained a temporary `vector<uint64_t>` of geometry-map keys.
Signed negative IDs widened to large unsigned integers, then implicitly narrowed
when looked up in maps keyed by int32. Native local tests recovered the expected
bits, but the narrowing was implementation-defined and bypassed the explicit
SignedID/IDBits conversion policy.

Changed the temporary vector and its iteration variable to NodeID. No encounter
state ABI or physics changed. Source-bound Serial and MPI builds and all native
suites were rerun after this fix; earlier passes are retained as pre-review
records, not the final evidence.

## Reviewer conclusion after fix

No remaining blocker found. Verified unsigned packing, representable conversion
through int64, bitwise ancestry chronology, matching MPI_INT32_T payloads, signed
thermo/event/trace parsing and restart identity checks. Negative IDs remain valid;
zero alone marks missing links. Ponytail verdict: lean already; no extra
abstraction suggested.

Independently reran 12 geometry-reader tests, strict curve IO/adversarial checks
and git diff --check: all passed.

Optional regression suggested: a legacy unsigned v4/B0 checkpoint with
UINT32_MAX and next_sequence=2^32. Implemented as ID08 in id_check.py, independently
parsed and resumed through the production C++ path; -1 and counter2^32 persist.
The final zero-bit boundary suite therefore has8 passing groups. Multi-rank
negative-ID round trips remain unqualified; compile success is not runtime proof.
Human scientific review remains pending.
