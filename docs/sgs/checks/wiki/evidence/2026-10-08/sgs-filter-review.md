# Independent M7 preparation review

2026-10-08. Reviewed and edited only `tst/sgs/filter_reference.py` and
`tst/sgs/test_filter_reference.py` under the assigned ownership. No C++ changes.
Runtime: `/Users/hitesh/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3`.

## Scientific convention review

The positive normalized periodic separable top-hat filter preserves constants.
Favre velocity is filter(rho u)/filter(rho); stress is the negative covariance
filter(rho) uFavre_i uFavre_j-filter(rho u_i u_j), with SGS K=-trace(tau)/2.
Velocity gradient storage is A_ij=partial_j u_i. Symmetric stress contracts with
A exactly as with S. The independent Python mixed-stress expression matches
published Eq. 8 after cancellation of the paper's gradient-norm factor two.
It retains the isotropic exact-zero-gradient convention and rejects negative K
or norm-range failure rather than silently applying a denominator epsilon/floor.

The filter's even-width offsets center the output half a cell below the input
sampling point. All fields use the same shift. Added an independent width-two
Fourier oracle cos(pi/n)*sin(phase-pi/n) and the discrete filter/centered-gradient
commutation identity, both passed. Delta=width*dx is the declared top-hat width,
not a calibrated equivalence to the paper's Gaussian-filter convention.

## Added preparation diagnostics

Added signed Pearson production correlation; exact/model forward and backscatter
fractions; false-sign fraction over nonzero exact production; separate forward,
backscatter and filtered-density subsets. Density subsets split at their median,
with ties in the lower subset, and include sample counts, mean density and signed/
RMS errors. Undefined constant-field correlation, empty subsets and zero-reference
ratios are explicit null values. A zero prediction counts as a sign error for a
nonzero exact reference; exact-zero cells are excluded from that fraction.

Independent tiny-array metrics oracle: exact=(-2,-1,1,2), model=(-1,-2,2,1)
has correlation 0.8 and no sign errors. Model=(1,-2,2,-1) has correlation 0 and
sign error fraction 0.5. With densities=(1,2,3,4), the lower-density signed error
is +1, upper-density signed error -1, and both RMS errors sqrt(5). Constant and
empty-subset behavior is independently asserted.

## Executed evidence

- Original eight synthetic tests passed before edits.
- Ten tests passed after edits, with frozen absolute tolerance 1e-12 unchanged.
- CLI synthetic n=8,width=3 output serialized finite JSON and retained
  `validated_turbulent_reference=false`.
- Synthetic trace residual 2.220446049250313e-16; production signed correlation
  0.6949144293221277, relative RMS error 1.6947466754679696, sign-error fraction
  0.41015625 and signed production bias +0.1888103036788796. These errors are
  retained; they are not turbulent-reference accuracy measurements.
- Focused flake8 attempt could not run because the bundled runtime has no flake8
  module. No dependency installed solely for this packet.

Durable evidence: `/private/tmp/sgs-filter-review/tests.log`,
`synthetic-n8-width3.json` and `manifest.json` (runtime, frozen tolerance and
source hashes). Both test and CLI exit codes are zero.

## Boundaries

This completes bounded synthetic pipeline preparation only. It does not establish
M7 closure accuracy, turbulent reference convergence, sampling uncertainty or
filter-scale transferability. Direct covariance subtraction can lose accuracy under
very large Galilean boosts; current moderate-boost test passes, while strongly
boosted data would require stable covariance/conditioning work. No true turbulence
dataset was introduced, no coefficients were tuned, and no scientific validation
status was upgraded.
