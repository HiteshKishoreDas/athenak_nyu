# Bounded stage-consistent forcing verification

Local serial-double verification completed. The revised forcing restart reproduces the uninterrupted solution exactly in three tested lanes. These pilots do not qualify a stationary turbulence reference or demonstrate scientific closure benefit.

## Reproducibility

- Repository base commit: `36be17da5c3335795c287986cbd437ce786a6390`; production changes were supplied by the parent/source owners.
- Executable: `/private/tmp/athenak-sgs-release/src/athena`.
- Executable SHA-256: `ba80e7d18d590956a27bd940e6f79afc1f40c6c3ec8d0993f3a8cdafa40f9105`.
- Tracked source-diff SHA-256 at packet start: `0e45c52634d7f3867c563ebfbf19911ff330ad9b7e33f9158b8f11d631890dbd`.
- Exact source file hashes, command arrays, per-input hashes, logs, native checkpoints and measured checks are retained under `/private/tmp/sgs-forcing-verification/pilot-02`; actual constant-acceleration kernel oracle is under `oracle-02`.
- Same compiler/configuration as the serial Release build: Apple clang 17, native ARM64 double, Kokkos Serial, PROBLEM=sgs_verify. Native restart decoding is qualified only for unrefined serial double and validates block layouts, byte counts, active cell coordinates and variable count.
- All forcing simulations are 16³, blocks 8³, RK2/PLM/LLF, CFL 0.1, gamma 5/3, density 1, pressure 1, velocity zero, two passive scalars. SGS seed k=1e-6 when enabled. Inputs: dedt=0.1, tcorr=0.1, dt_turb_update=0.01, seed=42, nlow=1, nhigh=3, constant_edot=true, stage_consistent=true; sol_fraction either 1 or 0.
- Pilots reached requested t=1 under 5000-cycle/60-second caps; all new attempts use unique directories. Revised packet simulation and rejection runtime was about 26.4 seconds. Entire forcing evidence currently occupies about 65 MB. No broad campaign was run.

## Frozen checks

- Native represented inventory E_resolved+K changes equal actually applied `force-work` at scaled absolute 1e-10. Worst measured residual 3.0532e-16.
- Periodic momentum changes equal actual integrated impulses at scaled absolute 1e-11. Worst measured residual 8.6737e-18. Mass and passive-scalar inventories also conserved.
- SGS reservoir change equals its independent production minus dissipation ledger at scaled absolute 1e-8. Worst residual 1.3553e-20. Forcing work stays separate from these SGS source counters.
- Stored states strictly finite, rho>0, thermal energy>0, K>=0; history inventories agree with native-double states. No threshold was loosened.
- Restart active-cell fields, raw OU memory and precise force counters require scaled Linf 1e-12; RNG native bytes require exact equality.

## Actual kernel constant-acceleration oracle

`problem/forcing_oracle=true` exercises the actual StageForcing production kernel twice with RK2 blending, acceleration (0.7,-0.2,0.3), dt=0.01. Starting rho=1 and zero velocity, the independent expected momentum increments are (0.007,-0.002,0.003), expected kinetic and total-energy increments are 3.1e-5, and thermal energy, SGS reservoir and passive scalar fields stay unchanged in the isolated kick. The pgen asserts every conservative field and global work/impulse counters at 1e-12, then restores state. Reported max error was 5.09141e-16. The no-requested-work follow-on step also passes ordinary inventory checks. This is a source-local manufactured kernel oracle, not an independent turbulent-box reference.

## Pilot measurements

| Lane | Actual work | Final resolved kinetic | Final thermal | Final SGS | Late KE fractional drift |
|---|---:|---:|---:|---:|---:|
| sol_sgs1_seed42 | 0.0998990828 | 0.0352356632 | 1.564597238 | 6.718140129e-05 | 0.22032 |
| sol_sgs0_seed42 | 0.09973038869 | 0.03529243009 | 1.564437959 | 0 | 0.22234 |
| comp_sgs1_seed42 | 0.09842576503 | 0.01788023187 | 1.58052348 | 2.305295368e-05 | 0.26016 |
| comp_sgs0_seed42 | 0.09563473661 | 0.01701476601 | 1.578619971 | 0 | 0.22227 |

Late trends are straight-line fits over the last six history samples, approximately t=0.5 to 1.0, and describe evidence sufficiency only. Kinetic energy is still increasing substantially. No stationarity, decorrelation, convergence, spectral fidelity or closure benefit is inferred. Same nominal dedt and random seed do not imply equal realized work across SGS settings: actual work differs because the state and normalization moments evolve differently. These paired lanes are a controlled implementation smoke test, not a controlled equal-work physical comparison.

## Restart and rejection results

- forcing.restart:sol_sgs1_restart: all fields/raw OU/RNG/work/impulses exactly equal (measured difference 0). Checkpoint t=0.25128827005726406, cycle=178, saved OU updates=25; final t=[1.0, 1.0], cycles=[753, 753], OU updates=['100', '100'].
- forcing.restart:sol_sgs0_restart: all fields/raw OU/RNG/work/impulses exactly equal (measured difference 0). Checkpoint t=0.25236553238038451, cycle=63, saved OU updates=25; final t=[1.0, 1.0], cycles=[276, 276], OU updates=['100', '100'].
- forcing.restart:comp_sgs1_restart: all fields/raw OU/RNG/work/impulses exactly equal (measured difference 0). Checkpoint t=0.2510391622746172, cycle=175, saved OU updates=25; final t=[1.0, 1.0], cycles=[756, 756], OU updates=['100', '100'].
- SGS on with stage_consistent=false is rejected with the required diagnostic. Changed dedt, kpeak, rseed, stage mode, OU counter and integrated work are rejected before evolution. An unknown new CLI key is also rejected by ParameterInput before the immutable-contract check.
- Revised packet ledger: 26 pass, 2 retained harness parsing/diagnostic failures, 0 blocked. The retained failure expected the immutable-contract diagnostic for an unknown CLI key, but the correct earlier ParameterInput rejection occurred. The diagnostic expectation was corrected, and all eight rejection cases were rerun under unique `_diagnostic02` names and pass. The revised physical runs and three restart comparisons have no failures. A supplemental offline check initially expected thermal-E from standard SGS-off histories; this parsing failure is retained as post_inventory_parser01. Reading thermal as tot-E minus component KE for that standard layout passes for all lanes, and both passive scalar native inventories are exactly conserved in every pilot/restart lane (post_inventory_parser02).

## Retained failures and correction

`pilot-01` retained four physically admissible pilot budgets and one failed normal intermediate OU restart: native active-cell scaled Linf 0.002514200399265, raw OU scaled Linf 0.02683886403526, RNG bytes still exact. Its checkpoint was t=0.25128827005726406, cycle178, OU count25; both final runs reached t=1, cycle753, OU count100. No-advance readback retained under `noadvance-01` demonstrated exactly identical state/raw OU/RNG payload at the checkpoint. Saved used dt=0.0013536942146648895 versus initialization-recomputed next dt=0.0013535879480367143 is expected from output-before-next-timestep ordering and by itself does not identify a restart defect.
The source owner identified a derived forcing parameter, default kpeak=4*pi, being printed as 12.5664 by ordinary parameter serialization. Rebuilding the forcing spectrum from rounded values changed the resumed trajectory despite exact raw OU and RNG payload. The approved source fix writes every actual stage-forcing Real parameter with max_digits10. `pilot-02` repeats all affected physics lanes and now achieves exact restart equivalence in solenoidal SGS on/off and compressive SGS on. Neither original failure nor diagnostic failure was deleted.

## Scope still open

This packet verifies the selected stage forcing, actual work/impulse accounting, separate SGS source budgets and normal restart memory on a serial 16³ startup pilot. It does not qualify forced MPI/decomposition, mixed solenoidal fractions, anisotropic/density-weighted forcing, a second seed, 32³/refinement, long-time stationarity, strong-forcing safety or force cutoff timing. A validated turbulence reference with adequate statistics and resolution remains missing. Scientific M7/M8 reference and accuracy/benefit acceptance remain incomplete; controlled non-turbulent SGS verification is reported separately in `/private/tmp/sgs-verification-report.md`.

## Harness files

- `/Users/hitesh/hitesh/git/athenak_nyu/tst/sgs/verify.py`: SHA-256 `0f64ad1510f455282a64ac231447856a23fcfbf2feb4778ffcc806eb558207a2`.
- `/Users/hitesh/hitesh/git/athenak_nyu/tst/sgs/verify_forcing.py`: SHA-256 `575375f5adfc8870cbc5705fe231e237a9130175efe9e72d852ac18c786cb6cc`.
- `/Users/hitesh/hitesh/git/athenak_nyu/tst/sgs/thresholds.json`: SHA-256 `5f4dd868584a213c263a7bf6e3c150f4d7eb23c5c58d07bf4734b019d34fa2f8`.

## Fresh official frozen-harness packet

The parent reran the complete current reusable harness after the diagnostic/parser
corrections under `/private/tmp/sgs-forcing-final`. Independent readback confirms
**18 pass, 0 fail, 0 blocked**: four physical pilot lanes, three resumed physical
lanes, three exact endpoint restart comparisons and eight expected rejections.
All seven physical records contain both distinct passive scalar inventory checks;
worst scalar difference is zero. All previously reported failing attempts remain
in their original directories/ledgers. The official packet executable SHA-256 is
`f8b29a46a2e7d6014a1e73227a322dcac6dce056208b25688a5a8b955830c3c5`. Exact source and harness hashes are in its
manifest. This fresh packet supersedes the corrected replay for clean acceptance,
without deleting the failure history.

## Final Debug bounds check and source review

Latest committed source `d6c5512edbb70fb38b650c78e28044a73a564325` was checked
with `/private/tmp/athenak-sgs-debug/src/athena`, PROBLEM=sgs_verify, Debug,
Kokkos Serial and Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON. The actual StageForcing
constant-acceleration oracle passed at max error 5.09141e-16, and the single
no-requested-work follow-on step passed native inventory, positivity, thermal and
passive scalar checks. Evidence is `/private/tmp/sgs-debug-forcing-oracle`; exact
Debug binary SHA-256 is `cae1a61ff497e1d24062c0e4aa416096e014f7da00f926dd2db90a78612b7135`.

Final independent read-only audit reconfirmed StageForcing precedes boundary
exchange, conserved-to-primitive refresh and timestep calculation; legacy
post-step forcing is omitted in the new mode. Raw OU write/read selection,
explicit version/count metadata, precise derived spectrum defaults and immutable
restart forcing parameters agree with the successful exact endpoint tests. No
new material source/runtime defect was found in the tested lane. Strong forcing,
cutoff timing and alternate weighting/MPI layouts stay unqualified as recorded
above.

## Bounded extension to t=2

The same frozen binary and all production source hashes match the clean t=1 packet. This fresh extension used the unchanged four 16³, seed42 pilot lanes, three normal intermediate restarts and eight rejection checks. **18 passed, 0 failed, 0 blocked** under 5000-cycle/60-second-per-case caps; aggregate execution time was 61.6 seconds. No source or harness was edited for this extension. Exact commands/input/source hashes and measured checks are in `/private/tmp/sgs-forcing-extended`.

| Lane | t=2 actual work | Final kinetic | Final thermal | Final SGS | Late KE fit drift, t≈1–2 | Fit fractional rate: earlier → extended | Late KE range / mean |
|---|---:|---:|---:|---:|---:|---:|---:|
| sol_sgs1_seed42 | 0.1997975026 | 0.04207958783 | 1.657422124 | 0.0002967907684 | 0.17450 | 0.44176 → 0.17453 | 0.17820 |
| sol_sgs0_seed42 | 0.1994670879 | 0.04217417149 | 1.657292916 | 0 | 0.17551 | 0.44574 → 0.17581 | 0.17887 |
| comp_sgs1_seed42 | 0.1974086899 | 0.01768824726 | 1.679611282 | 0.0001101601999 | 0.12878 | 0.52100 → 0.12883 | 0.91088 |
| comp_sgs0_seed42 | 0.192847412 | 0.01734053015 | 1.675506882 | 0 | 0.13980 | 0.44766 → 0.14014 | 0.91958 |

Fractional fit rate has units of inverse code time and divides the fitted fractional drift by its actual window length. The earlier t=1 window is roughly 0.5–1 (six samples); this t=2 window is roughly 1–2 (eleven samples). Rates decrease, but the fitted t=2 drift is still 12.9–17.6%. Compressive kinetic energy fluctuates strongly: late-window range is 91–92% of the window mean, so the positive fitted slope does not imply monotonic growth. No stationarity acceptance threshold was invented after observing these values.

All applied-work/represented-energy checks pass, with worst residual 6.384e-16; momentum/actual-impulse worst residual 1.055e-17. Both passive scalars are exactly conserved in every physical lane. Three continued/restarted endpoint states, raw OU memory, RNG bytes and integrated work/impulse counters again match exactly. The nominal requested injection remains 0.1, while actually applied work differs between SGS-on/off lanes; a physical equal-work comparison is still unqualified.

Only seed42 was executed. A three-seed campaign, 32³ expansion, long-time independent statistics and a validated turbulence reference remain not run: the M7 reference and stationarity/statistics gates are unmet. The additional run is a bounded local pilot, and scientific M7/M8 acceptance is unchanged.

Extended manifest commit: `14620bb33d580246e32d80f4bd584c6c021bbd1e`; executable SHA-256 `f8b29a46a2e7d6014a1e73227a322dcac6dce056208b25688a5a8b955830c3c5`.
