# Independent read-only review

The three established reviewers continued their respective equation, integration
and diagnostic audits. No reviewer changed source or launched native simulations.

- Equation reviewer: the global pressure-work sign and half-face normalization
  are correct. Independent Fourier compression and DC-shear residual formulas
  establish initial native amplitude. Fixed-space smooth refinement is a numerical
  reference and does not isolate turbulent RK heating. Suggested prespecified
  order/difference conditions were frozen and met.
- Integration reviewer: actual reconstructed pressure is captured in x/y/z,
  including fallback. Existing final RK quadrature, MPI reduction and precise
  restart/header persistence are used. No fluid RHS modification found. Their
  history-capacity finding was fixed by increasing both capacities from 32 to 33;
  the previously supported forced 12-scalar native case passes.
- Diagnostic reviewer: independently reproduced both complete 32³ 2T–4T budgets
  and four short-CFL fractions exactly. Actual native checkpoint metadata must
  be used for continuation labels/T, not unused generated input; this was fixed.
  Both 64³ checkpoint hashes and exact origin inventories/counters match protocols.
  Original checkpoints lack pressure counters, correctly defining new counter 0 at
  restart; no pre-origin pressure-free history may be concatenated. Undefined
  fractions map to omitted plot bars. No blocking schema/code issue found.

The short-CFL probe establishes only observed 2T–2.1T sensitivity at 32³ with
coupled state evolution and interpolated endpoints. No full-window or 64³ temporal
convergence claim is supported. Final 64³ fraction review is recorded below after
successful bounded completion. The competing numerical channel is combined
spatial/RK/shock-associated heating, not an LLF-only measurement.

## Completed 64³ review

Independent raw native-history calculations including the verified origin reproduce
all seven 2T–4T increments exactly. The combined report has identical 64³ rows.
PLM SGS/numerical fractions are .5663938910981392/.43360610890186074;
PPM4 fractions are .7005373918142387/.29946260818576126.
Represented/forcing window residuals are -8.75e-14/+2.29e-13 and SGS reservoir
residuals 3.66e-15/6.33e-15. The reviewer independently compared all active fields,
forcing amplitudes/full metadata and RNG with the original final checkpoints:
exact agreement, time 2. Both bounded native lanes pass all recorded checks within
1800 seconds. PPM4 approaches the approximate paper 75/25 reference; one-seed
2T–4T statistics, unmatched forcing and the combined spatial/RK/pressure-work/shock
residual preclude full reproduction or an LLF-only claim. No blocking reporting
issue found. This review remained read-only; no additional native jobs were run.
