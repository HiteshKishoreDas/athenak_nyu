# Fixed-forcing reconstruction comparison

Local source commit `7cdff3618726c3ac386a8bb1efa7ec32f0b7d6a2` on
`codex/multiphase-schmidt-federrath-hd`. Only the reconstruction support guard
changes production code: existing native monotone PPM4 is allowed, alongside
DC/PLM. LLF, RK2, the SGS stencil/flux/source equations and closure coefficients
remain fixed. PPMX/WENOZ/HLLC/FOFC remain rejected. Three ghost zones are required
by the existing PPM implementation; both comparison lanes use three. The
paper runner gains explicit reconstruction/ghost options and a bounded native
qualification check. Independent review found no blocking defect. No push/merge.

## Forcing audit and correction

[Schmidt–Federrath §4](https://arxiv.org/html/1010.4492) fixes V=L/T.
For this run L=.5 and actual OU correlation time T=.5 imply V=1. The earlier
algebraic V=.661 rescaling that aligns normalized K cannot reinterpret the same
forcing: it would require T=.756. No existing results were renormalized.
The pure-solenoidal Table 4 does not explicitly state V/c0 or an absolute modal
amplitude. [Schmidt 2009 equations 7–9](https://arxiv.org/pdf/0809.1321) give
asymptotic acceleration RMS=3V²/L, or 6 for our scales, whereas the native driver
imposes instantaneous volume RMS 3. Adoption of that exact convention for the
SGS paper's Table 4 case is not established. Other cited forcing methods use
different conventions. Thus the reproduction remains Mach-matched with incomplete
forcing equivalence; the current coefficients and forcing were not retuned.
The table dissipation header remains inconsistent with nearby rho-epsilon prose;
the density-weighted interpretation is supported by prose, not a confirmed erratum.
[Detailed primary-source audit](sgs-forcing-convention-audit.md).

## Executed controls and results

Matched fresh-from-rest 32³ and 64³ lanes use identical active fields, physical
inputs, seed 42, Fourier OU process, signed radial modes 1–3, pure solenoidal forcing,
instantaneous accel_rms 3, correlation T=.5, gamma 1.001, c0=1/5.5, k0=1e-8,
C1=.02,C2=.7,Ceps=1.5,Ckappa=.65, RK2/LLF, CFL .1 and block 16³.
Actual input files differ only in hydro/reconstruct. Exact first fields and final
modal OU/RNG comparisons pass. Final physical time is 2=4T in both lanes.
Snapshots are spaced .1T; averages use actual timestamps clipped 2T–4T.
The completed campaigns pass native admissibility, inventories, applied forcing work,
paired SGS exchange, ordinary scalar, momentum and history checks with the
existing frozen thresholds. These are one-seed short-window pilots, not the
paper's full 2T–10T statistical campaign, seed-ensemble or convergence evidence.

| 64³ quantity, average 2T–4T | PLM | PPM4 | Relative change |
|---|---|---|---|
| Resolved RMS Mach | 5.49394 | 5.51484 | +0.38% |
| SGS RMS Mach | 1.04634 | 1.15413 | +10.30% |
| Mean instantaneous SGS/resolved Mach ratio | 0.190448 | 0.209299 | +9.90% |
| Mean SGS energy / (rho0 V²) | 0.0321678 | 0.0371305 | +15.43% |
| Density-weighted dissipation L/(rho0 V³) | 0.36822 | 0.45931 | +24.74% |

At 64³, thermal energy grows 4.1239% with PLM and 4.1106% with PPM4 over
2T–4T. Applied SGS dissipation increases from .737925 to .918734 in raw
code units (+24.50%); the independently sampled rho-epsilon average changes
by +24.74%. Higher SGS energy therefore accompanies higher transfer, while
resolved Mach and thermal evolution remain close. This supports numerical
method sensitivity; it does not assign the whole paper discrepancy to numerics.

32³ pilot changes: K +21.66%, mean instantaneous Mach ratio +9.03%, density-weighted
dissipation +34.62%; resolved Mach +1.87%. This grid lies below the paper's resolution
floor. Do not infer a resolution law by comparing the two grids.

![Controlled 64³ comparison](64-comparison.png)

Blue=PLM, orange=PPM4. The first five panels compare the same driven case and
fixed coefficients. The bottom-right panel is a separate stationary shear oracle,
with zero SGS energy/coefficients and no forcing. It quantifies numerical loss;
it does not measure a turbulent heating fraction. Thermal growth in turbulence
includes resolved pressure work, shocks and numerical dissipation.

The smooth Euler shear has rho=1, uniform pressure, vx=sin(2pi y), vy=vz=0 and
periodic boundaries. Its exact continuum solution is stationary. All nine native
runs pass. At 64² and t=1, PLM loses 0.05608% of kinetic energy versus PPM4 0.001679%,
about 33.4 times less; heat gain equals kinetic loss within 1e-10. All SGS transfers
are exactly zero. Refinement 16²/32²/64² and DC control are retained. This supports
a reconstruction-dependent numerical-loss channel, while turbulence sensitivity
alone does not isolate LLF's separate effect or prove the paper discrepancy solved.

## Qualification, cost and evidence

The sixteen native PPM qualification records pass: uniform decay, zero/tiny K,
contact, diffusion, 3D box/decomposition, restart, collision and forcing. The
collision activates 9 fallback stages. Configuration 34/34 records pass, including expected rejection of PPM4
with two ghost zones, PPMX/WENOZ/HLLC/FOFC. Independent
C++ closure oracle, cellwise diagnostics and irregular-time averaging selfchecks
pass. The PLM ghost-layout check uses inherited tolerances. No threshold changed.

64³ lanes run concurrently with four OpenMP threads each (eight total); each
has 2400 s native cap and 100000 cycle cap. 32³ caps 600 s with eight threads per run;
small checks/reporting overlapped. Timings are not controlled performance benchmarks.

The original 64³ PLM attempt hit its 2400 s cap and remains FAIL/incomplete. It resumed from t/T=3.80029 with an additional fixed 1200 s cap, four threads, and identical physics. Exact-load, every-prefix budget, continuation and assembled full-window checks pass. The original results/workflow files remain failed; [PLM recovery ledger](64-plm-recovery.json) records the completed recovered campaign separately. The executed continuation uses only `-r` plus basename/time-limit overrides; generated default harness input files are unused. Native checkpoint headers are the actual config.

The original 64³ PPM4 attempt hit its 2400 s cap and remains FAIL/incomplete. It resumed from t/T=3.30005 with an additional fixed 1200 s cap, four threads, and identical physics. Exact-load, every-prefix budget, continuation and assembled full-window checks pass. The original results/workflow files remain failed; [PPM4 recovery ledger](64-ppm4-recovery.json) records the completed recovered campaign separately. The executed continuation uses only `-r` plus basename/time-limit overrides; generated default harness input files are unused. Native checkpoint headers are the actual config.

Every run retains native inputs, checkpoints, histories, hashes and logs. Large
raw data/executable remain outside wiki Git in `/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-scheme-raw-evidence.tar.gz`; SHA256
`d18bb20898ba68825626dff961c580c37e7647ce085f93e407c44ad21a91f09c`. All 700 archived regular members
independently match source hashes and sizes. Size 7217122021 bytes.
The prior paper and decay records remain frozen.

[64³ six-panel PDF](64-comparison.pdf), [32³ comparison](32-comparison.pdf),
[64³ PLM 18-panel report](64-plm-report.pdf), [64³ PPM4 18-panel report](64-ppm4-report.pdf),
[comparison statistics](64-comparison.json), [independent review](sgs-scheme-independent-review.md),
[exact commands](commands.json), [local runner instructions](reproduce.sh).

## Run locally

Run from the feature checkout at the source commit recorded above:

```sh
cd /Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu
PYTHON=/private/tmp/sgs-plot-venv/bin/python sh \
  /Users/hitesh/Documents/das-wiki/science/athenak-subgrid-turbulence/evidence/2026-10-09-scheme-comparison/reproduce.sh
```

The existing local Python environment supplies NumPy and Matplotlib. A different
Python with those packages can be selected through PYTHON. The script builds
and qualifies PPM, then runs the matched 64³ cases sequentially with eight
threads. It produces each native 18-panel report. Runners refuse to overwrite
existing output directories; use fresh paths for repeated runs. This is a
reproduction recipe, not an additional executed campaign.
