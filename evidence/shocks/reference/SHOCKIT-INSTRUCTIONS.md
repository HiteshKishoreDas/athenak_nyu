# AMR development documentation

User requirement: maintain aggressive development documentation in
`amr_support_log.md` at the repository root.

- Read the log's current checkpoint and pending work before AMR development.
- Update the log during every AMR task and before reporting completion.
- Record requests, decisions and rationale, changed files/behavior, experiments,
  exact checks and results, failures/corrections, limitations, and next steps.
- Distinguish historical evidence, fresh verification, hypotheses, and plans.
- Preserve older entries; document corrections explicitly. Do not silently
  rewrite failed experiments or unresolved accuracy issues as successes.
- Link evidence using repository-relative paths. Keep user data and prior output
  directories intact. Do not weaken tests or tolerances merely to obtain a pass.
- Keep the log synchronized with the relevant plan and acceptance documents.

This is a documentation requirement, not an authorization to run unrelated work
or create background automations.

# Human verification is the release gate

- Prefer concise, purpose-revealing names; retain familiar scientific notation
  when its meaning is documented. Preserve public keyword compatibility.
- Document function contracts and explain non-obvious numerical steps nearby:
  coordinate frames, units, assumptions, boundaries, limiting and mask semantics.
- Keep control flow explicit. Concision must not conceal the scientific argument.
- Passing automated tests is necessary, never scientific approval. Release or
  scientific use requires recorded human review of the exact code and evidence.
- Maintain `docs/human_verification.md`; leave sign-off pending until a human
  reviewer supplies it. Never infer approval from a request to implement changes.
