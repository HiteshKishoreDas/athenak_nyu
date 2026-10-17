# feature_dump integration record

## 2026-10-17 — material curves

- Request: create `feature_dump` from `multiphase_turb`; merge
  `codex/multiphase-material-curves`; preserve the source branch.
- Worked directly in `/Users/hitesh/hitesh/git/athenak_nyu`. Initial checkout was clean;
  previous branch `codex/staircase-ordered-commits` remains intact.
- Refreshed `origin` with `git fetch origin`. Local and remote base/source tips agree.
- Base: `526a187f938e38634e89a88192ebe9731c4d36bb`.
- Source: `a0744a0d6dbc93364201bd8f5248e6be35b938b9`.
- Created with `git switch -c feature_dump multiphase_turb`.
- Merged with `git merge --no-ff codex/multiphase-material-curves`.
  Explicit merge records the feature integration and its original base.
- Merge commit: `b13ca7d65fc821653d04cf006b6bfcf6e4b61e9e`. No conflicts or code edits.
- Verified merge parents and exact tree equality with the source feature tip
  before adding this integration record. Source branch remains at its original tip.
- Check: `PYTHONDONTWRITEBYTECODE=1 python3
  tst/scripts/material_curves_io/test_reader.py` — all 10 tests passed.
- Existing numerical evidence: [material curves](../doc/material_curves.md) and
  [qualification record](../doc/material-curves-evidence/README.md).
  No new build, simulation, GPU qualification, push or PR performed in this merge.
