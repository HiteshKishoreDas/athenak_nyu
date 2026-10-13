# Independent review, 2026-10-09

Read-only specialist audit of 7cdff361: monotone PPM4 reconstructs all variables
including specific SGS energy and ordinary scalars; nonnegative inputs have
nonnegative interface values. LLF receives SGS-dependent wave speeds in all
directions. Existing whole-pack donor-cell fallback recomputes fluxes, stresses
and rates. Native constructor rejects fewer than three PPM ghost zones. PPMX,
WENOZ, HLLC and FOFC remain rejected.

Qualification: uniform, zero/tiny K, contact, diffusion, three-dimensional box,
block decomposition, native restart, collision and forcing checks all pass with
frozen thresholds. Hard collision activated 9 fallback stages. Configuration
34/34 checks include PPM4/ng2 rejection. Independent C++ closure oracle passes.

Independent 32-cubed review: actual inputs differ only in reconstruction;
initial active fields equal and final OU/RNG equal. Time-weighted averages use
2T–4T. K +21.66%, SGS Mach +11.08%, SGS/resolved Mach ratio +9.03%, density-
weighted dissipation +34.62%. Thermal growth 4.13%/4.16%, so the temperature
change does not account for this reconstruction difference. This is one seed,
below the paper resolution floor, and a short window.

Stationary shear is a valid separate Euler oracle, with K=0 and SGS coefficients
zero. At64 squared, PLM numerical kinetic loss0.05608% versus PPM4 0.001679%.
This does not measure a turbulence numerical-heating fraction.

Reviewer recommended deleting a duplicate title and asserting exact reconstruction
labels; both changes were applied and the 32 comparison regenerated in
sgs-scheme-32-comparison-checked. No blocking numerical/scientific defect found.
Full paper forcing equivalence and LLF-specific dissipation remain unresolved.

## Final 64-cubed independent review

Original 2400-second timeout results remain FAIL/incomplete and their hashes
match recovery provenance. PLM resumes at 3.8002895546T; PPM4 at 3.3000522031T,
with fixed 1200-second continuation caps and four threads each. No threshold,
executable, source or physical parameter changed. Native -r commands use only
basename/time-limit overrides; generated harness inputs are unused.

Recovered PLM 532/532 and PPM4 492/492 checks pass. Every prefix checkpoint
passes admissibility and independent budget checks; exact-load fields, OU/RNG
and counters match. Both assembled views contain 41 native/history records
covering 0–4T; maximum physical gaps .051204/.051206 are below .055. Parent,
probe and continuation duplicates agree before deduplication. Raw archive
dereferences view symlinks and retains all target artifacts.

Initial active fields are identical; actual initial configurations differ only
in hydro/reconstruct. Final OU amplitudes, RNG and full forcing metadata match
exactly at physical time 2. Independent clipped-trapezoid means, without using
the reporting statistics function, reproduce all six diagnostic means.

2T–4T means: K .03216780045/.03713052930 (+15.4276%); resolved Mach
5.49394/5.51484 (+.3804%); SGS Mach 1.04634/1.15413 (+10.3015%); mean
instantaneous SGS/resolved Mach ratio .1904476/.2092991 (+9.8985%).
Raw density-weighted dissipation .7364406/.9186201 (+24.7379%); normalized
by L=.5 the plotted means are .3682203/.4593101. Applied SGS dissipation
means .7379248/.9187340; production .7359435/.9163612. Thermal growth
4.1239%/4.1106%, almost identical. No fallback or diffusion limiting.

Plot is readable; no blocking defect found. This is one-seed reconstruction
sensitivity with LLF retained, averaged over a short 2T–4T window. It does not
establish forcing equivalence, full paper reproduction or scientific release.
