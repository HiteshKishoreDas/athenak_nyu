# Paper Section 3 filtered-snapshot closure test — 2026-10-12

Ran fresh 32³/64³ SGS-off periodic Newtonian HD turbulence with the existing
supersonic setup: gamma 1.001, RK2/LLF/PPM4, CFL .1, initial rho0=V=1,
seed 42, solenoidal modes 1–3, fixed acceleration RMS 3. Runs reach 10T with
the forcing law unchanged. Nine nominal integer-T snapshots from 2T through 10T
are analyzed; exact timestamps and hashes are in full-window-selection.json.
The snapshot comparison has no adjustable forcing parameter: reference and
model evaluate the identical flow.

## Method

Following [Schmidt & Federrath §3](https://arxiv.org/html/1010.4492#S3), filter
density, momentum and all six second moments with the identical operator:

    rho_bar = filter(rho)
    u_Favre = filter(rho*u)/rho_bar
    tau_ij = rho_bar*u_Favre_i*u_Favre_j - filter(rho*u_i*u_j)
    K_reference = -trace(tau)/2
    measured_transfer = tau_ij * derivative_j(u_Favre_i)

Supply this measured K, rho_bar and gradient to AthenaK's actual C++ Stress
kernel with fixed C1=.02/C2=.7. Compare its tensor and transfer with the measured
quantities. Do not supply evolved SGS K. Gaussian filters use the explicitly
chosen Delta_G=sqrt(12)*sigma convention, with positive periodized sampled
weights and unit normalization. Box filtering provides a sensitivity check.
Physical widths are box/16, box/8 and box/4. Sampling remains on the fine grid,
and primary derivatives are centered differences. Spectral derivatives provide
a separate sensitivity check at box/16. No coefficient fitting or clipping.

## Results: nine-snapshot means over 2T–10T

| Reference | Filter | Width / box | Model / measured mean transfer | Mean correlation | Mean relative RMS error | Mean K_reference |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| 32³ | gaussian | 1/16 | 0.9351 | 0.9966 | 12.3% | 0.059623 |
| 32³ | gaussian | 1/8 | 0.9275 | 0.9892 | 15.4% | 0.146410 |
| 32³ | gaussian | 1/4 | 0.9419 | 0.9646 | 25.1% | 0.296176 |
| 32³ | box | 1/16 | 0.9608 | 0.9963 | 12.7% | 0.049823 |
| 32³ | box | 1/8 | 0.9396 | 0.9842 | 18.0% | 0.140439 |
| 32³ | box | 1/4 | 0.9284 | 0.9446 | 31.4% | 0.297545 |
| 64³ | gaussian | 1/16 | 0.9292 | 0.9949 | 13.1% | 0.082437 |
| 64³ | gaussian | 1/8 | 0.9042 | 0.9848 | 17.3% | 0.172061 |
| 64³ | gaussian | 1/4 | 0.9014 | 0.9568 | 27.6% | 0.315709 |
| 64³ | box | 1/16 | 0.9383 | 0.9929 | 14.2% | 0.078678 |
| 64³ | box | 1/8 | 0.9074 | 0.9761 | 20.9% | 0.170516 |
| 64³ | box | 1/4 | 0.8799 | 0.9322 | 34.8% | 0.318966 |

The primary 64³ Gaussian box/16 comparison has mean transfer
-7.1% bias, correlation
0.99489 and relative local RMS error
13.1%. Full correlation ranges
0.99276–0.99687 across
the nine snapshots. Removing the identical isotropic pressure contribution
leaves mean trace-free production correlation
0.98972. This is supporting evidence
for the production closure on these flows, and does not show a factor-of-two
production defect. It does not validate the evolved SGS reservoir or the native
face-flux discretization. Correlation alone does not establish amplitude accuracy.

Mean ratio is the ratio of nine mean fluxes, not the mean of nine ratios.
Correlation and relative RMS error are means of nine per-snapshot metrics, not
pooled fits or independent-seed uncertainty. Full-volume fields supply metrics;
maps show central-z slices. The plotted joint-distribution/PDF range is limited
for readability; all metrics include the complete volume and tails.

## Checks actually run

- 1134/1134 comparison arithmetic checks across 108 filter/snapshot cases pass.
  Native-versus-independent-Python stress disagreement is at rounding level;
  direct positive real-space centered covariance is checked at two cells per
  case. Checks also cover kinetic-energy decomposition, mass/momentum,
  covariance positivity, trace/symmetry and Galilean invariance.
- 590/590 independent native inventory checks cover all 82 stored checkpoints,
  including applied forcing work/impulse, scalar mass and final histories.
- Both native continuation reports pass 11/11 checks and finish exactly10T.
- 12/12 existing/extended filter unit tests pass. Native hand-case, Gaussian
  attenuation, spectral derivative-axis and truncated-record checks pass.
- The existing native closure oracle passes; C++ style, focused Python static
  checks and whitespace checks pass. A read-only independent scientific and
  ponytail review found no blocking issue in the pilot comparison.

The initial generic report harness used basename reference while looking for
iles32/iles64 history filenames. Both simulations actually completed and saved
valid native data, but those original report attempts remain FAIL in
native/results.json. Independent validation reads their actual history names;
subsequent native continuations use the existing forcing-aware harness. No
scientific input was changed to repair the report. Initial/pilot/full analysis
packets are retained. Final figures remove overlapping log-axis minor labels;
the corrected rendering retains identical scientific metrics.

## Limits and unresolved energy question

The primary Gaussian width is four fine cells at 64³, but only two at 32³.
The latter is explicitly flagged as undersampled. These small grids do not
establish an inertial range or reproduce the paper's 768³–1024³ reference data.
Reference mean K changes 27.7% between 32³ and 64³ at the
same explicit physical width. This fails the earlier proposed 10% reference
sensitivity screen, so reference convergence remains unestablished.

**K_reference=0.082437 is energy removed at filter width box/16.**
It is not a direct target for the evolved64³ reservoir at grid spacing box/64,
and must not be compared to .0367689 as though both had the same cutoff.
Gaussian-paper width equivalence and the actual finite-volume filter remain
unconfirmed. Supplying exact K deliberately removes reservoir evolution from
the production test. Dissipation calibration and native spatial energy transfer
still require separate validation against an adequately resolved reference.

No solver equations, closure defaults, source splitting, forcing or restart
format changed. New code is a Gaussian analysis extension, a host utility that
calls the existing kernel, an analysis runner, independent checks and usage docs.
No new dependency. Native athena SHA256 remains
b5f0fc19de082e07bc45c989cc24d5cc0c46dfd4d62e7994d8dd1efaac901197.
The source branch is codex/multiphase-schmidt-federrath-hd. Focused local code
commit is recorded in commits.json; the wiki commit is reported in the chat.
Nothing was pushed.

## Reproduce on the Mac

Build the test helper and run its smallest standalone checks from the SGS checkout:

```sh
cmake -S . -B /private/tmp/athenak-sgs-heating -DBUILD_TESTING=ON
cmake --build /private/tmp/athenak-sgs-heating --target sgs_filtered_stress -j2
/private/tmp/sgs-plot-venv/bin/python -B tst/sgs/test_filter_reference.py
/private/tmp/sgs-plot-venv/bin/python -B tst/sgs/verify_filtered_snapshot.py --self-test \
  --kernel /private/tmp/athenak-sgs-heating/src/sgs_filtered_stress
```

analysis-command.json contains the complete argument list for all 18 archived
snapshots. Change its --output value to a new directory, then run the list with
Python subprocess.run(...,check=True). This regenerates every metric and the
fourteen plot panels from saved data, without rerunning native simulations.
See tst/sgs/README.md for the portable CLI command with your own snapshot paths.

For fresh native reproduction, use the saved native/iles32/input.athinput and
native/iles64/input.athinput in new working directories, overriding time/tlim=5.
Their grid, SGS-off, gamma, forcing and output parameters are explicit. Use two
OpenMP threads and retain the original native binary or build matching source.
The original segments and continuations together used about 15 min of wall time
in this run; timing includes output and validation and is not a benchmark.

## Plots

![Measured versus native transfer and covariance](full-window-checked/snapshot-comparison.png)

![Time, trace-free production, derivatives and backscatter](full-window-checked/time-diagnostics.png)

![Filter and reference-resolution sensitivity](full-window-checked/filter-sensitivity.png)

Each figure also has a PDF in the same directory. Raw checkpoints remain in
the original artifact tree; compact copies and SHA256 provenance are kept in
the wiki evidence packet. Source, binary, input and selected snapshot hashes
appear in the comparison protocol/results and manifests.
