# Forced supersonic turbulence: a particle-encounter showcase

Open supersonic-showcase.pdf for the12-page report. Its10figure pages put explanations
next to each plot. index.html is the gallery; individual vectorPDF/SVG and250dpiPNG
exports are ready to edit. episodes.csv contains every committed and terminal-open
identity; evolution.csv contains the field diagnostics. flow-and-mask-data.npz holds
the exact selected primitive volumes and native detector arrays used in the maps.
This is a reproducible numerical implementation showcase for paper development;
human physical calibration/approval remains pending.

## Experiment and measured results

Uniform periodic unit cube, baseline64³/8192native mass-seeded MCparticles,
initialrho1/P1/1.4/Fourier rmsMach5. Ideal gamma1.4, RK2/PLM/HLLC/CFL.30,
continuous native edot50forcing(modes1..3/npeak2,sol_fraction.5,tcorr.1,
update.01,seed2718), tlim.30. Initial driving length.5 and initialvrms5 define
an initial turnover unit.10; the run spans3such initial units, without stationarity
being asserted. Velocity/particle/schedule seeds24680/11/903. Width4/minMach1.1/
capacity32768/tracefalse. These are experimental showcase settings, not approved
physical-error tolerances. No cooling/source/mesh/runtime scope was broadened.

| Case | Complete | Partial | Open | Complete / particle | Final volume rms Mach |
| --- | ---: | ---: | ---: | ---: | ---: |
|64³, CFL.30,width4|57013|0|1212|6.9596|1.31598|
|64³, CFL.15,width4|60813|0|1169|7.4235|1.31686|
|32³, CFL.30,width4|35045|0|1198|4.2780|1.26339|
|32³, CFL.30,width2|44472|0|1529|5.4287|1.26339|

Baseline median first-inside pressure-Mach estimate1.59865;10th/90th percentiles
1.13878/2.83846. Final61.7043percent of baseline cell volume is locally supersonic.
All saved baseline volume rms Mach values exceed1. All8192IDs survive, no invalid
observations occur, canonical identities are unique and per-particle sequences
contiguous. Strict-reader counts agree with native tracker counters. All episodes
retain no_plateau/cadence_limited;5703baseline identities have merged-flicker and
1520have direction ambiguity.No extra shock counts are inferred from suppressed
attempts: the actual surviving-state suppressed-entry counter is75455, while emitted
downstream_reentry flags are0. Counter and flags describe different populations.

Five actual simulations complete: baseline95.343s, tracking-off62.317s,
half-CFL186.398s, coarse6.723s, matched-coarse7.300s. These are individual local
wall times with startup/output; matched-coarse partly overlapped half-CFL. They are
not replicated production performance benchmarks. Baseline tracking-on/off entire
native thermo file bytes and all54saved primitive/conserved/force volumes match
exactly. The two coarse runs' particle bytes and27saved volumes match exactly,
isolating detector width on that flow. Periodic mass error stays<=1e-12.

## Read the statistics correctly

Counts refer to numerical encounter episodes, not a certified census of distinct
physical shocks. Rate timestamps are first-inside observations; no interpolated
crossing time is invented. Duration span bounds are[last-first,post-pre]for uncensored
completed episodes. Consecutive complete sequence gap bounds are[max(0,next.pre-
previous.post),max(0,next.first-previous.last)]. Merged spans may contain brief
outside samples. Open outcomes are separately right-censored and can resume;
no finite complete duration/survival fit is assigned. All raw flag populations,
zero-encounter particles and censored identities remain in CSV/report data.

The map outlines are native full-region detector recomputations from actual saved
primitive states at generations150/400/802, not a compression surrogate or inferred
live surface. The existing fixture executable supplies genuine periodic neighbors.
All262144cells are valid; inside counts29183/46119/51957. Fields use native zyx
storage, transposed to the fixture's xyz contract and back. Snapshot histories align
by generation/full-double native particle clocks, not six-digit binary metadata.
Fixed endpoints/temperature normals and first-inside strengths are not independently
plateau-certified; existing ShockIt parity is port fidelity, not physical truth.

HalvingCFL changes completed counts+6.665percent, jointly through fluid/MC RNG/
cadence. Coarse width4 doubles physical detector half-span. Coarse width2 matches
64³width4half-span(1/16box length), but fluid/MC resolution still differs. In the
identical coarse flow, width2 produces26.900percent more completed identities than
width4. These are sensitivity results, not convergence acceptance or uncertainty
bars. No independent-seed ensemble or stationary waiting-time law is claimed.

The known unchanged native RK2 gas-energy forcing oracle still FAILS; its independent
shear evidence is in ../forcing-enabled-20261013/. Evolving adiabatic temperature
and strengths can reflect that source bias. Actual forcing work, conservative shock
heating and source bias are not equated. No new physical sign-off, GPU/MPI/AMR or
RUSTY/Frontier qualification is implied. Use these figures with the adjacent methods
and limitations as implementation evidence; production scientific results need their
separate calibration/ensemble/resolution/hardware gates.

## Reproduce

Coreproductionad4225edea78d255b67df97fe8f593c0c268fa1e, fingerprint
f6195adafb7e6a554c23a7537053bb00e6ae8c618a10ff1951a487d6981cb70b.
Showcase input/analysis/PDF tools commitbac7241d59d56fee565bf9c25e3a0d03d24a8ee5.
No production source was modified for this showcase. These commits are local,
not newly pushed. The actually executed launcher is retained; the portable runner
performs the equivalent5cases sequentially. Logs bind the actual configuration
(including older comments); the current input correctly labels the actual Mach5.

From /Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu:

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DAthena_SHOCK_TESTS=ON
cmake --build build-release -j8
MPLCONFIGDIR=/tmp/shock-turb-mpl \
/Users/hitesh/hitesh/git/shockit/.venv/bin/python tst/shocks-mc/supersonic_showcase.py \
  --executable build-release/src/athena \
  --detector build-release/shock_detector_test \
  --workdir /tmp/showcase-new --output /tmp/showcase-figures-new
/Users/hitesh/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  tst/shocks-mc/supersonic_showcase_pdf.py /tmp/showcase-figures-new
```

Use fresh work/output directories. The first Python uses already installed
NumPy/Matplotlib, native binary/thermo/strict event readers. The second uses bundled
ReportLab/pypdf and existing explained_pdf helpers. Each simulation is capped600s;
three native detector snapshots are capped180s each. --reuse reanalyzes a completed
matching-executable execution.json in an existing work directory, with a new output
directory. No rerun occurs during PDF composition. Full raw runs/native histories,
thermo data and detector text/CSV fixtures remain/private/tmp/supersonic-showcase-20261013.
Reproduce those full native histories by rerunning; the committed paper packet holds
CSV statistics, selected numerical arrays, figures, inputs/logs and full raw hashes.
Previous layout/color preflight is retained, all measurements unchanged after polish.
Scientific timing/statistics and integration reviewers checked the definitions and
independently recomputed data. All12final PDF pages were rendered/visually inspected.
