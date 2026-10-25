# Independent configurable curve-bit review

Reviewer: `/root/review_curve_id_persistence`, read-only; ponytail-review applied.
No reviewer file edits. Scope: runtime shifts, allocation/exhaustion, checkpoint
metadata/compatibility, reader decoding, consistency and isolated trace behavior.

Confirmed zero-bit code avoids shifting a uint32 by 32. Limits use uint64, and the
per-curve allocator stores the full $2^{32}$ exhausted sentinel without wrapping.
Curve-node/link storage and MPI ID payloads stay uint32. Counts/digests/generations
and shared encounter state retain their existing widths.

Required strict string parsing because AthenaK's general GetOrAddInteger uses
atoi: garbage/overflow could become valid zero. The coordinator had replaced it
with the existing checked unsigned parser before initial builds. Eight negative
input checks include text, fractions, negative, overflowing and out-of-range bits,
plus multiple/wrong-ID curves in zero-bit mode.

Two reader findings were reproduced and fixed:

1. A standalone B=4 node `0xF0000001` could be guessed as curve 240 instead of 15.
   Current traces require a manifest or an explicit curve_id_bits argument for
   32-bit records; historical 64-bit decoding is unambiguous. Explicit arguments
   must agree with manifests. Regression covers correct decoding, missing metadata,
   explicit arguments and disagreement.
2. Empty direct event readers skipped the bit parser because it ran only inside
   the row loop. Manifest parsing now validates the declared layout for every
   history/suspended/trace API, including empty journals. Regression rejects B=31
   in an empty journal. Chunk parsing also computes layout once, outside its loop.

Final result: no remaining blockers. Reviewer independently reran the 11 geometry
reader tests, standalone curve IO/adversarial checks and diff checks successfully.
Geometry, thermo and journals declare layouts; cross-stream/within-stream and
restart split changes reject. Historical analysis and default-eight-bit geometry
compatibility remain explicit. No material abstraction/dependency bloat found.
The coordinator's native source-bound gates and hardware limits are in the ledger;
this review does not supply scientific approval or MPI runtime qualification.
