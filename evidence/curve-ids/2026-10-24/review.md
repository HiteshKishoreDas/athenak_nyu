# Independent compact-ID persistence review

Reviewer: `/root/review_curve_id_persistence`, read-only; used ponytail-review.
Scope: curve ID allocation, persistent links, reader/version boundaries, MPI ID
payloads and shared shock-state lookup. No reviewer file mutations.

Initial required checks: allow the allocator sentinel without wrap; parse config
IDs into a wide temporary before checking; parse persisted links directly at their
new width; reject wide shared encounter IDs before compact-map lookup; change
only MPI ID datatypes, preserving wide counts/digests/generations; do not treat
unsigned compact IDs as signed THP tags. Implemented and tested.

One correctness finding: using a trailing substring to select compact encoding
could misclassify a historical contract containing the same text in a user
parameter. Reviewer reproduced this with an actual synthetic event. Fixed at the
contract root: current code emits leading `curve-node-history-v2`, and readers
inspect that prefix. Added historical misleading-text and current wide-ID
rejection tests. Reviewer confirmed the fix.

Reviewer requested a valid-checksum checkpoint mutation adding one unsigned-32
span to a real node ID. Added to CSH06: it repairs FNV and requires the exact new
range-guard diagnostic. Native reviewed CSH20 suite passes with that case included.

Final result: no remaining blockers. Reviewer independently ran all 11 Python
geometry-reader tests, standalone curve IO/adversarial checks and diff whitespace
checks successfully. No abstraction/dependency bloat found. Native suite results
are separately recorded by the coordinator; MPI runtime was not qualified.
