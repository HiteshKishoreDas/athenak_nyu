---
title: SGS test catalogue execution ledger
class: living
created: 2026-10-08
last_touched: 2026-10-08
---

# SGS test catalogue execution ledger

Every catalogue family is listed. **Pass means the concrete instance named in the scope column passed; it does not close every combination or a milestone.** Partial means useful evidence exists but the stated oracle or coverage is incomplete. Blocked identifies an execution or dependent reference gate. Not run means no relevant executed oracle was assigned. Earlier failed attempts remain failed in their raw result ledgers; a later passing instance does not erase them.

See [implementation and limits](implementation-progress.md), [the original catalogue](test-matrix.md), [machine-readable family ledger](coverage-ledger.json) and [raw evidence index](evidence/2026-10-08/README.md). Detailed per-instance commands, inputs, thresholds, source hashes and metrics live in the archived results, rather than pretending a family row is a run manifest.

Coverage count: 230 families; 65 pass, 55 partial, 22 blocked, 88 not run.

| ID | Class | Status | Executed scope or remaining gate |
|---|---|---|---|
| C01 | B | pass | Clean built-in SGS-absent regressions; configuration omission. |
| C02 | B | pass | Omission versus explicit disable: bitwise fields/layout. |
| C03 | B | partial | Effective coefficients printed; full defaultconfiguration coverage open. |
| C04 | B | partial | Multiple true scalars andcleanzero-scalar seed; wider combinations open. |
| C05 | B | pass | Unknown model rejected. |
| C06 | B | pass | Negative coefficient/nonfinite seed rejected; finite endpoint kernel cases. |
| C07 | B | pass | C2>1 rejected. |
| C08 | S | partial | Isothermal early rejection tested; no heat deposition support. |
| C09 | S | partial | Source rejects multiphysics; not every rejected block executed. |
| C10 | S | pass | AMR rejected with valid mesh setup. |
| C11 | S | pass | Stretched cells rejected. |
| C12 | S | pass | One ghost zone rejected. |
| C13 | S | pass | RK3/HLLC/WENOZ rejected. |
| C14 | B | not run | No executed instance assigned. |
| C15 | S | pass | Mode/coeff/scalar restart overrides rejected. |
| C16 | S | partial | Cooling/viscosity/legacyFOFC rejection tested; remaining hook/options sourceaudited. |
| K01 | B | pass | Exact-zero matrix and finite K oracle. |
| K02 | B | pass | Nonzero shear/zero K oracle. |
| K03 | B | partial | Independent component matrix/random oracle; complete hand-matrix family open. |
| K04 | B | pass | Isotropic compression stress/production hand oracle. |
| K05 | B | pass | Traceless strain hand oracle. |
| K06 | B | partial | One analytic shear orientation plusgenericrandomtensor oracle; all orientations notqualified. |
| K07 | B | pass | Antisymmetric rotation hand oracle. |
| K08 | B | pass | 1000 random tensors, symmetry/trace/component checks. |
| K09 | B | pass | State/width powers and random positive state scales. |
| K10 | B | pass | C2=0 endpoint. |
| K11 | B | pass | C1=0 endpoint. |
| K12 | B | pass | C1=C2=0 retains isotropic pressure. |
| K13 | B | pass | Ce=Ckappa=0 endpoint. |
| K14 | B | partial | Boostedcontact conditioning tested; full transformed stress/work solution open. |
| K15 | B | pass | One independently specified orthogonal rotation. |
| K16 | B | not run | No executed instance assigned. |
| K17 | B | pass | Hand forward and negative production. |
| K18 | B | pass | Gradients down to1e-280 normalized safely. |
| K19 | B | pass | Density/K combinations1e-100 and1e100. |
| K20 | B | pass | Gradients up to1e150 in declared finite range. |
| K21 | B | not run | No executed instance assigned. |
| K22 | B | pass | 1025-case independent oracle versus Kokkos Serial. |
| K23 | S | blocked | SGS precision support restricted todouble. |
| K24 | S | blocked | CUDA/HIP unavailable; noGPU kernelexecution. |
| K25 | B | partial | Reduced PDE cases and3D rawtrace; complete reduced-gradient matrix coverage open. |
| K26 | D | not run | No executed instance assigned. |
| L01 | B | partial | Conversion through varieddensity runs; standalone face/conversion roundtrip open. |
| L02 | B | pass | Distinct true scalar masses/fractions in independent native-double reads. |
| L03 | B | not run | No executed instance assigned. |
| L04 | B | pass | Clean built-in uniform seed .03 exact bookkeeping. |
| L05 | B | pass | Exactly zero K remainszero. |
| L06 | B | pass | Same-pgen normal checkpoint state not reseeded. |
| L07 | B | partial | Transport smoke/conservation; exact translation convergence open. |
| L08 | B | partial | Contact cases pass; allconstant-specific variants/convergence open. |
| L09 | B | partial | ConstantK densitycontrast cases andnonlinear transport smoke; isolatedoracle open. |
| L10 | B | partial | Difficultprofiles remainnonnegative; deliberate negativefacestate trigger unexercised. |
| L11 | B | partial | Periodicity exercised; individualaxis crossing oracles open. |
| L12 | B | pass | 16^3 block8^3 versus16^3 active state bitwise. |
| L13 | B | not run | No executed instance assigned. |
| L14 | B | partial | Multiblock3D stencils boundschecked; targeted corner-crossing oracle open. |
| L15 | S | not run | No executed instance assigned. |
| L16 | S | not run | No executed instance assigned. |
| L17 | S | not run | No executed instance assigned. |
| L18 | S | not run | No executed instance assigned. |
| L19 | S | not run | No executed instance assigned. |
| L20 | B | not run | No executed instance assigned. |
| P01 | B | pass | Uniform zeroK supported1/2/3D. |
| P02 | B | partial | Uniformoff-rate seed cleanbuiltin; full moving variants open. |
| P03 | B | pass | Independent decay and heat law. |
| P04 | B | pass | Three CFL decay orders1.95053/2.00308. |
| P05 | B | partial | Decaysteprestriction used; deliberately largeststable source cases open. |
| P06 | B | not run | No executed instance assigned. |
| P07 | B | not run | No executed instance assigned. |
| P08 | B | partial | Conditional shear component timeconvergence; exact discrete/PDE momentumoracle open. |
| P09 | B | partial | Rawproduction handstates andcoupledbudgets; prescribed-gradient localevolution notrun. |
| P10 | B | partial | Nonlinear transport smoke budgets; isolated analytic spatial convergence open. |
| P11 | B | partial | Filterunits reviewed; specificconstant diffusion operator oracle notrun. |
| P12 | B | partial | Densitycontrast smoke; isolated variablek flux oracle open. |
| P13 | B | not run | No executed instance assigned. |
| P14 | B | not run | No executed instance assigned. |
| P15 | B | pass | Phase-shifted longitudinal C2 branch, four reservoir errors16/32; C1=Ce=Ckappa=0. |
| P16 | B | partial | HighK/lowP wave stability/budget smoke; waveformoracle notthisregime. |
| P17 | B | partial | Compatible semidiscrete identity reviewed + separate component timeconvergence; full exactpartition open. |
| P18 | B | not run | No executed instance assigned. |
| P19 | B | partial | Decay/shear timeorder; fullmixed manufactured temporal order open. |
| P20 | B | partial | Longitudinalbranch16/32 refinement; full mixedspatial asymptotics open. |
| P21 | B | pass | Heun-weighted applied exchanges independently close K change. |
| P22 | B | pass | Periodic E+K inventories over executed bounded steps. |
| P23 | B | pass | History versus native-double active-cell sums. |
| P24 | B | pass | Actual isolated StageForcing constant-acceleration oracle plus measured box work. |
| P25 | S | not run | No executed instance assigned. |
| P26 | D | partial | TinyK contacts boosted0,.3,100; fullGalilean transformedsolution open. |
| R01 | B | partial | Active highK coarsebackscatter cap; tinyK insufficientdonor variants open. |
| R02 | B | partial | Hardcollisions andweakthermal wave admissible; independent intenseproductionpartition open. |
| R03 | B | not run | No executed instance assigned. |
| R04 | B | partial | Actual combinedpredictor fallback; deliberately isolatednegativeK trigger open. |
| R05 | B | not run | No executed instance assigned. |
| R06 | B | not run | No executed instance assigned. |
| R07 | B | partial | SGS-owned globaldonorcell fallback observed; legacyFOFC explicitlyrejected. |
| R08 | B | partial | Globalrebuild sharedfacesconsistent in testedblocks; targeted adjacenttrigger case open. |
| R09 | B | partial | Multiblock collisions; explicitrank-interface trigger blockedMPI. |
| R10 | B | partial | Multi-stage fallbackcounts observed; isolatedstage-target instrumentation open. |
| R11 | B | pass | Fresh clean built-in Sod smoke and baseline wave regression unchanged. |
| R12 | B | partial | Colliding shocksfiniteK passbudgets; SGSshocktube convergence open. |
| R13 | B | partial | v5/v20 collisions activate9/19 fallbacks; detailedshockaccuracy open. |
| R14 | B | not run | No executed instance assigned. |
| R15 | B | partial | Boost100 contacts characterizeconditioning; near-cold highMach translation open. |
| R16 | B | not run | No executed instance assigned. |
| R17 | B | partial | CFL.1/.05/.025 andhardcollision.2 tested; broaderregime sweep open. |
| R18 | B | partial | Receivingdensity rowbound anddiffusion CFLcases; resolution-ladder limit oracle open. |
| R19 | B | partial | Stress/decay limits +trial checks; generalnonlinear stability remainsunproved. |
| R20 | B | partial | HighSGS pressure firststep smoke passes; isolatedconstraint audit open. |
| R21 | B | partial | Initial/restart timestep recomputed andnormalcontinuation passes; changedlocalconstraint scenario open. |
| R22 | B | blocked | No multi-rank launcher execution. |
| R23 | B | pass | ZeroK has noextraSGS timestep limit in supported cases. |
| R24 | B | partial | Combined source/transport predictors passselectedcases; explicitminimum oracle open. |
| R25 | D | not run | No executed instance assigned. |
| R26 | D | partial | Contact/shear laminar responses recorded; broadcharacterization open. |
| R27 | D | not run | No executed instance assigned. |
| R28 | B | pass | Three temporary sign/factor mutations detected by independent oracle. |
| I01 | B | pass | Selected clean built-in SGS-absent regressions, not fullsuite. |
| I02 | B | pass | Selected Debug bounds unforced, forcing and wave checks. |
| I03 | B | partial | Bothbuilds wave/sourceoracles agree; systematicall-field Debug/Release comparison open. |
| I04 | B | pass | Two serial block decompositions. |
| I05 | S | partial | MPIcompile/directrank1 passes; launcherABI prevents2/4ranktests. |
| I06 | S | blocked | Multi-rank launcher ABI mismatch. |
| I07 | S | blocked | Threadedbackend not configured/qualified. |
| I08 | S | blocked | CUDA/HIP unavailable on arm64Mac; SGS rejects unqualifiedGPU. |
| I09 | S | blocked | SGS singleprecision explicitlyunqualified; no claimednumericalpath. |
| I10 | B | pass | Normal unforced checkpoint bitwise. |
| I11 | B | pass | Three forcing continued lanes bitwise fields/rawOU/RNG/counters. |
| I12 | S | not run | No executed instance assigned. |
| I13 | B | pass | Native-double SGS state and precise metadata retained. |
| I14 | B | pass | Clean legacy scalar checkpoint restart bitwise. |
| I15 | S | pass | SGS off/on restart mismatch rejected. |
| I16 | B | partial | Native-double restart andhistory independentlyread; binaryselectors codeadded fullreader/formatqual open. |
| I17 | B | pass | Independent inventories and actual forcing/work ledgers. |
| I18 | B | not run | No executed instance assigned. |
| I19 | B | blocked | No multi-rank global extrema/counter tests. |
| I20 | S | partial | Selectedlegacyhydro/scalar tests; cooling-index/fullbranch regressions notrun. |
| I21 | B | partial | Focusedchanged-line C++ checks andPython syntax pass; flake8 unavailable. |
| I22 | B | pass | Normal checkpoint after3 fallback stages reproduces9 finalstages and fields. |
| I23 | B | partial | Exactnormalrestart andsourceoracle repeats; fullstochastic startup repeats open. |
| I24 | S | pass | Independent clean archive build at14620bb3 + pinnedKokkos. |
| F01 | B | pass | Synthetic constant top-hat. |
| F02 | B | pass | Odd/even sinusoid transfer/phase oracle. |
| F03 | B | pass | Hand density-weighted tiny arrays. |
| F04 | B | pass | Independent negative covariance and trace. |
| F05 | B | pass | Independent tensor/strain production contraction. |
| F06 | B | pass | Conditioned independent model evaluated on synthetic filtered fields only. |
| F07 | B | partial | Evenfilter center andgradientcommutation synthetic test; turbulent sampling sensitivity open. |
| F08 | B | blocked | No referencewithdemonstrated converged turbulent scales. |
| F09 | B | blocked | No validated turbulentreference regimes. |
| F10 | S | not run | No executed instance assigned. |
| F11 | B | partial | Synthetic signedmetricstests; turbulent forward/backscatter accuracy open. |
| F12 | B | partial | Synthetic densityconditionedmetrics; turbulent Mach/compression/shear subsets open. |
| F13 | B | blocked | No turbulent independentrealization/time referencecampaign. |
| F14 | B | blocked | No referenceconvergence ladder. |
| F15 | D | not run | No executed instance assigned. |
| F16 | B | blocked | No held-out turbulentdata. |
| F17 | B | pass | Synthetic sign/density weighting defect checks. |
| F18 | D | not run | No executed instance assigned. |
| T01 | B | pass | Capped16^3 measured-work SGSon/off pilots; implementation smoke only. |
| T02 | B | not run | No executed instance assigned. |
| T03 | B | pass | Ordinarygamma5/3 heating/partition recorded in bounded pilots. |
| T04 | B | blocked | Reference/stationarity gates incomplete. |
| T05 | B | blocked | Reference/stationarity gates incomplete. |
| T06 | B | blocked | Reference/stationarity gates incomplete. |
| T07 | B | blocked | Stationarity andreference gates incomplete; onepilot seedonly. |
| T08 | B | blocked | Pilotstationarity insufficient; no frozenstatistical sampling. |
| T09 | B | pass | Actual sourcework independently closes E+K. |
| T10 | B | not run | No executed instance assigned. |
| T11 | B | not run | No executed instance assigned. |
| T12 | B | not run | No executed instance assigned. |
| T13 | B | blocked | No qualifiedcommonreference scale/campaign. |
| T14 | B | not run | No executed instance assigned. |
| T15 | B | blocked | No stationarymulti-resolution uncertaintyfit. |
| T16 | B | not run | No executed instance assigned. |
| T17 | B | partial | Pilotcounts/budgets recorded; longcampaign localfrequency/volume qualification open. |
| T18 | S | not run | No executed instance assigned. |
| T19 | D | not run | No executed instance assigned. |
| T20 | D | not run | No executed instance assigned. |
| T21 | B | blocked | No validatedfilteredreference or equal-work/cost comparison. |
| T22 | B | pass | Impulseledgers close measured momentum, not assumedzero forcing. |
| T23 | D | not run | No executed instance assigned. |
| T24 | B | blocked | No finalscientificcomparison to independentlyreproduce. |
| H01 | B | not run | No executed instance assigned. |
| H02 | B | not run | No executed instance assigned. |
| H03 | B | not run | No executed instance assigned. |
| H04 | S | not run | No executed instance assigned. |
| H05 | D | not run | No executed instance assigned. |
| H06 | B | not run | No executed instance assigned. |
| H07 | B | not run | No executed instance assigned. |
| H08 | B | not run | No executed instance assigned. |
| A01 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A02 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A03 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A04 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A05 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A06 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A07 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A08 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A09 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A10 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A11 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A12 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A13 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A14 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A15 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A16 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A17 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A18 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A19 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A20 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A21 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A22 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A23 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A24 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A25 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A26 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A27 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A28 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A29 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A30 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A31 | A | not run | Deferred extension; noimplementation authorized before M9. |
| A32 | A | not run | Deferred extension; noimplementation authorized before M9. |
| X01 | X | not run | Deferred extension; noimplementation authorized before M9. |
| X02 | X | not run | Deferred extension; noimplementation authorized before M9. |
| X03 | X | not run | Deferred extension; noimplementation authorized before M9. |
| X04 | X | not run | Deferred extension; noimplementation authorized before M9. |
| X05 | X | not run | Deferred extension; noimplementation authorized before M9. |
| X06 | X | not run | Deferred extension; noimplementation authorized before M9. |
| X07 | X | not run | Deferred extension; noimplementation authorized before M9. |
| X08 | X | not run | Deferred extension; noimplementation authorized before M9. |


## Scoped addition: filtered-snapshot closure test, 2026-10-12

Source `b513bc9d` adds the native Section 3 comparison. Fresh SGS-off 32³/64³
supersonic runs reach 10T with fixed forcing; nine nominal 2T–10T snapshots
produce 108 Gaussian/box comparisons. All 1134 arithmetic, 590 independent
native inventory and 22 continuation checks pass, as do 12 filter unit tests
and the existing native closure oracle. Primary 64³ Gaussian mean transfer is
7.1% low, correlation .995, relative local RMS error 13.1%, with trace-free
production correlation .990. Original report-name failures remain recorded.

This is a useful turbulent production pilot, not closure or reservoir scientific
qualification. Gaussian width is explicit but not certified paper-identical;
reference K changes 27.67% between grids at fixed physical width. F08/F09/F14
reference qualification and evolved SGS-energy validation remain open. Do not
convert those baseline blocked gates to pass based on the arithmetic counts.
[Protocols, results, failures and reproduction](evidence/2026-10-12-filtered-snapshot/README.md).

## Scoped addition: unforced evolved-energy comparison, 2026-10-12

Test-only source `7250a4e4` adds matched initially Mach9 decay through10T:
128³/64³ SGS-off references and8³/16³ SGS-on/off. Six runs and488 numerical
checks pass, including native inventories, exact initialization, admissibility,
periodic conservation and paired SGS transfers. Separate native reader and
short-window reporting checks pass. Three earlier harness failures are retained.

The1–10T model/reference mean K ratios are .75/.79; at10T .45/.44. However,
reference sensitivity reaches17.52%/34.61%, failing the predeclared10% screen.
The fine reference ends subsonic, at Mach.899. This is a useful unforced-energy
pilot, not qualified reference accuracy or a paper reproduction. F14/T21
scientific qualification remains incomplete; historical family gates are not
promoted to pass. No model coefficient or solver equation changes were made.
[Protocol, plots, checks and interpretation](evidence/2026-10-12-unforced-decay/README.md).

## Scoped addition: higher unforced reference, 2026-10-12

Test-only source `921cf832` supports resource-controlled256³ decay without
changing native solver physics. The matched256³/128³ references and8³/16³
SGS-on/off lanes reach10T; six runs and488/488 numerical checks pass. All324
native source hashes and the executable match provenance. Pilot physical fields
are bitwise identical after block/thread/passive-scalar resource changes.

Mean SGS/reference K over1–10T is .69840/.66505; final ratios .41885/.34980.
Maximum reference sensitivity is9.67%/21.52% at the8³/16³ filters. Only the8³
filter meets the predeclared10% screen; this is not full DNS convergence proof.
Its30% mean energy deficit remains larger than this reference sensitivity.
The16³ reference qualification remains open. Late fine flow is subsonic,
Mach.896; one seed, small coarse grids and undeveloped initial conditions remain
limitations. F14/T21 and M9 scientific qualification are not promoted to pass.
[Inputs, metrics, plots, source verification and reproduction](evidence/2026-10-12-higher-reference/README.md).

## Scoped addition:64³/128³ SGS decay, 2026-10-12

Test-only `ab54c37d` adds guarded reference reuse and unavailable same-grid
convergence reporting. New64³ on/off and128³ on trajectories complete10T;
the matched256³/128³ references are revalidated, and128³ off is the shared
reference/control. Six analyses,488 checks and105 unique field snapshots pass.
All324 native source hashes and the executable remain unchanged.

Mean model/reference K ratios1.02266/1.93662 and final .45964/.83733 vary with
resolution and time. Total kinetic ratios are .95410/.98233 on versus
.95532/.98176 off. At64³, reference sensitivity59.06% fails the10% screen;
at128³ the screen is unavailable because no second finer reference exists.
These results neither qualify absolute SGS accuracy nor support a universal
factor-of-two deficit. F14/T21 and M9 scientific gates remain open.
[Metrics,28 plot panels, reuse checks and reproduction](evidence/2026-10-12-sgs64-128-decay/README.md).
