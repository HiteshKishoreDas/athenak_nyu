# Independent native-output reviews

2026-10-23. Two read-only agents used the requested ponytail-review skill.
Coordinator owned all edits and integrated sequentially. Shared interface:
`job/tracer_output_format=legacy|binary`, common `.tbin` v1 named typed tables,
unchanged `.thp` v2 and checkpoint state layouts.

## C++ persistence review

Reviewer: `review_native_output_persistence`.

- No blocking C++ persistence bug found. Explicit u64/i64 serialization,
  type/count validation, immutable commits and format-bound restart contracts
  preserve identities and digest reconciliation. Empty tables have headers/footer.
- Requested removal of discarded CSV/text construction in binary paths.
  Applied by gating legacy formatting while reusing geometry calculations and
  validation. No second geometry estimator was introduced.
- Flagged reduced geometry analysis metadata. Added boundaries and adaptation
  thresholds; documented that full AMR topology/refinement clocks stay in mesh
  output and existing checkpoint sidecars. No restart layout was changed.
- Static review does not qualify MPI empty ranks, hardware or production scale.

## Reader and clock review

Reviewer: `review_native_output_readers`.

- P2: declared binary rank count could exceed the surviving manifest set and
  still be accepted. Fixed by comparing each event/suspended table's declared
  rank count with the validated manifests in both carrier readers.
- P2: common reader silently skipped existing shock chunks if every manifest
  was missing. Fixed by rejecting orphan event/trace/birthtrace/suspended files.
- TO12 retains independent reproductions for both carriers. Reviewer rechecked
  fixes for events and suspended views and reported no remaining blocker.
- Exact integer widths, NaN masks, semantic-parser reuse and extracted output
  clock logic appeared sound. Clock extraction preserves the previous thermo
  advancement behavior; it is reused by the new scheduled geometry output.

## Coordinator verification

TO01–TO21 pass on the frozen source. They cover event/trace/geometry/thermo
record parity, conserved-fluid table byte equality, uint64 IDs above 2^63,
AMR births, restart and zero-step journal replay, both review failures, schema
corruption and independent clocks. The final 18 CT, 20 CSH, 14 LT and 9 startup
groups also pass. Human scientific review and deployment qualification remain
pending; these reviews do not supply that approval.
