---
title: AthenaK SGS turbulence — comprehensive test catalogue
class: living
created: 2026-10-08
last_touched: 2026-10-08
---

# AthenaK SGS turbulence: test catalogue

This is a planned coverage catalogue, not a result ledger. Every case starts
**not run**. Tests must have an independent oracle or a clearly identified
diagnostic question. Passing conservation does not establish closure accuracy;
passing a turbulent spectrum does not establish conservation or AMR correctness.

See the [development plan](development-plan.md) for equations, thresholds,
milestones and the support contract, and the [agent runbook](agent-runbook.md)
for execution and evidence requirements.

## Implementation status on 2026 10 08

Restricted HD code and independent verification are now recorded in [implementation progress](implementation-progress.md) and the [execution ledger](coverage-ledger.md). M0–M9 scientific release is incomplete. Original planning statements and checklists below are preserved as their historical record.

## Coverage labels and scheduling

- **B:** mandatory baseline coverage for the supported HD configuration.
- **S:** required when the associated solver, boundary, backend or mode is supported;
  otherwise test its explicit rejection and mark the numerical case not applicable.
- **D:** diagnostic/characterization; a limitation or null result can be the conclusion.
- **A:** later AMR extension gate, after baseline HD release.
- **X:** future physical extension; no implementation promise in the current scope.

An ID describes a test family; its execution ledger lists the actual parameter
instances. Kernel families run quickly. Smooth verification/robustness usually use
small grids. Reference/turbulent statistics are batch campaigns, not per-commit CI.
No unsupported case is silently marked passed. Every B family needs at least one
appropriate concrete instance; S families require every claimed execution path to
receive coverage, using pairwise selection where full combinations are unnecessary.

For each row, the gate means the stated oracle passes with the threshold frozen
for its precision, timestep, scale and conditioning. Differential comparison with
the same buggy implementation is insufficient on its own.

## C — configuration and support contract

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| C01 | B | Omit SGS configuration from an existing input | Identical legacy layout and supported baseline behaviour; no new required parameters |
| C02 | B | Explicitly disable SGS | Same state/counts/output meaning as omission; no closure kernels or tensor storage |
| C03 | B | Enable default SGS on the supported HD configuration | Effective coefficients, filter width and activation are visible in the run manifest |
| C04 | B | Enable with zero and multiple true passive scalars | SGS appended without changing passive identities; independent scalar initialization |
| C05 | B | Unknown closure name or malformed option | Actionable early failure, rather than fallback to another model |
| C06 | B | Negative/nonfinite coefficients, width or seed | Reject invalid physical/numerical inputs; accepted zero endpoints follow the declared contract |
| C07 | B | Out-of-range mixed coefficient | Reject unless a separately reviewed parameter-study contract permits it |
| C08 | S | Strict isothermal EOS with SGS heating requested | Explicit rejection or qualified documented mode; no claim to validate heat deposition |
| C09 | S | MHD, relativistic, radiation, two-fluid or incompatible coordinates | Reject unsupported coupling before evolution |
| C10 | S | Static refinement / AMR before A0–A2 are available | Explicit rejection; no silent passive-scalar LES remap |
| C11 | S | Noncubic cells | Reject or report the tested filter convention and associated support limits |
| C12 | S | Insufficient ghost zones for chosen reconstruction/SGS stencil | Early validation catches each stencil radius, including fallback and corner gradients |
| C13 | S | Unqualified integrator or Riemann solver | Reject or mark explicitly experimental; never imply it passed another solver's gate |
| C14 | B | Duplicate/overridden settings and command-line input changes | One effective setting with reproducible precedence and printed provenance |
| C15 | S | Restart changes SGS mode/coefficients/state counts | Defined compatibility policy, explicit rejection or reviewed migration; no silent reinterpretation |
| C16 | S | User hooks, viscosity, conduction, cooling or shearing settings | Detect unsupported combinations; qualified combinations have separate coupling tests |

## K — closure algebra and device kernels

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| K01 | B | Zero velocity gradient with finite positive SGS energy | Declared exact-zero stress convention; zero production; trace remains minus twice SGS energy |
| K02 | B | Positive gradient with zero SGS energy | Zero SGS stress/production/diffusivity/dissipation; finite output |
| K03 | B | General hand-computable 3×3 gradient | Independent componentwise stress and production calculation |
| K04 | B | Pure isotropic compression/expansion | Correct isotropic stress, compression sign and exponent; no duplicate compression term |
| K05 | B | Traceless diagonal strain | Analytic eddy and nonlinear contributions, including signed production |
| K06 | B | Simple shear in each coordinate orientation | Correct gradient index ordering and tensor contraction; no assumed nonlinear shear production |
| K07 | B | Antisymmetric solid-body rotation gradient | Symmetric stress; production vanishes because strain vanishes; stress need not be isotropic |
| K08 | B | Deterministic random tensors over a moderate range | Stress symmetric and trace equals minus twice SGS energy; deviatoric trace zero |
| K09 | B | Density/SGS-energy/filter-width scaling | Correct independent powers and units in stress, dissipation and diffusivity |
| K10 | B | Eddy-only coefficient endpoint | Matches independently coded eddy closure with the proper trace |
| K11 | B | Nonlinear-only coefficient endpoint | Matches independent normalized nonlinear formula with the proper trace |
| K12 | B | Both cascade coefficients zero | Isotropic pressure remains; no accidental disappearance of compression work |
| K13 | B | Zero dissipation or zero transport coefficient | Only that term vanishes; pressure/production/other transport remains |
| K14 | B | Add a uniform velocity offset | Gradient-based stress and production invariant; only transport/work-frame terms transform |
| K15 | B | Rotate coordinates by an orthogonal matrix | Tensor covariance, scalar production invariance and consistent component transformation |
| K16 | B | Reflect/permute coordinate axes | Expected parity and permutation of components; scalar invariants unchanged |
| K17 | B | Positive and negative production examples | Signed result agrees with oracle; negative transfer not unconditionally clipped to zero |
| K18 | B | Tiny positive gradient norms approaching zero along different directions | Finite, documented branch behaviour; trace preserved; no false continuous-limit assumption |
| K19 | B | Very small positive SGS energy and density | No invalid square root/division; dimensionless normalized results remain conditioned |
| K20 | B | Large gradients and energy within declared safe range | Scaling-safe evaluation; no overflow caused by avoidable intermediate products |
| K21 | B | NaN/Inf injected into input data in debug verification | Deliberate detection/reporting at the defined boundary; no misleading successful metric |
| K22 | B | CPU scalar oracle versus Kokkos kernels | Precision-aware component agreement; oracle derived separately from production helper |
| K23 | S | Single versus double precision | Conditioned algebra accuracy and meaningful relaxed thresholds; unsupported precision rejected |
| K24 | S | Available GPU backend kernel | Same hand/random cases and trace/production checks without host-only code assumptions |
| K25 | B | Reduced 1D/2D gradients with three velocity components | Three-dimensional stress trace convention preserved; missing derivatives exactly zero |
| K26 | D | Add time-dependent rotating-frame terms | Characterize model limitations; do not require full material objectivity absent from baseline |

## L — storage, transport, initialization and boundaries

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| L01 | B | Conservative/primitive conversion with nonuniform density | SGS primitive equals specific energy; reconstructed conservative energy restored correctly |
| L02 | B | Multiple ordinary scalars plus SGS state | True scalar values/names/indices unchanged; SGS slot independently transported |
| L03 | B | Existing metallicity/dust scalar input with SGS off | Existing branch regression remains intact; no first-scalar repurposing |
| L04 | B | Uniform state, prescribed positive seed | Initialization contains exactly the stated represented energy; no hidden subtraction/addition |
| L05 | B | Exact-zero seed | Defined startup behaviour reported; no untracked repeated seed injection |
| L06 | B | Restart through the same pgen | SGS state not reinitialized and history hooks correctly re-enrolled |
| L07 | B | Smooth SGS profile translated by constant velocity | L1/L2/Linf convergence and conservation; correct specific-energy mass-flux coupling |
| L08 | B | Density contact with uniform specific SGS energy | Specific energy preserved within transport error; density-weighted conserved SGS profile |
| L09 | B | Constant SGS energy density with varying density | Transport/diffusion responds to the correct variable; not mistaken for constant specific energy |
| L10 | B | Sharp nonnegative SGS profile | No negative reconstructed face states or undocumented cell clipping |
| L11 | B | Advection across periodic boundary in each axis | No seam, loss or duplication; global SGS integral conserved with other SGS terms disabled |
| L12 | B | One block versus many blocks | Same short-run transport accuracy and boundary flux cancellation |
| L13 | B | Uneven pack sizes and final partially populated pack | No uninitialized padded-block reads; correct reductions and allocation extents |
| L14 | B | SGS signal crosses edges/corners of block interfaces | Valid gradient stencil and ghost ordering in all directions |
| L15 | S | Reflecting wall | SGS scalar parity, velocity-gradient parity, zero normal mass/transport flux where prescribed |
| L16 | S | Outflow boundary | Declared SGS extrapolation and outward represented-energy flux; no spurious incoming seed |
| L17 | S | Prescribed inflow | Density, thermal state and SGS state mutually consistent; injection measured as boundary flux |
| L18 | S | User boundary function | Explicit SGS boundary contract exercised; uncovered fields rejected or diagnosed |
| L19 | S | Shearing-periodic/orbital remap | Separate state/gradient/energy-work contract and remap tests before support is claimed |
| L20 | B | Allocation/destruction over repeated small mesh setups | No leaked SGS storage, stale index or view extent mismatch |

## P — smooth PDEs, ODEs and compatible energy transfer

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| P01 | B | Uniform stationary HD, zero SGS energy | State preserved; no production/diffusion/dissipation or new floor event |
| P02 | B | Uniform moving HD, finite SGS energy; all SGS rates disabled | Uniform translation preserved, including pressure without a gradient |
| P03 | B | Uniform SGS energy, dissipation only | Exact constant-density decay law and matching thermal-energy increase |
| P04 | B | P03 at three or more timesteps | Expected temporal order; decay error not hidden by large internal energy |
| P05 | B | Large decay timestep | Nonnegative energy and bounded heat transfer under the selected stable solver; otherwise correct dt restriction |
| P06 | B | Prescribed homogeneous compression without decay/transport | SGS energy proportional to density to power 5/3; specific SGS energy to power 2/3 |
| P07 | B | Homogeneous expansion counterpart | Correct decrease and paired work sign; no negative SGS energy |
| P08 | B | Smooth transverse shear with isolated stress contribution | Momentum diffusion/transfer agrees with independent discrete or manufactured oracle |
| P09 | B | Production-only prescribed-gradient local evolution | Coupled source conserves the sum at every stage; forward/backscatter signs correct |
| P10 | B | Smooth SGS gradient, isolated transport | Flux points down specific-energy gradient; integral conservation and spatial convergence |
| P11 | B | Varying density with constant specific SGS energy | Zero SGS diffusion despite a conserved-energy gradient |
| P12 | B | Varying density and uniform conserved SGS energy | Correct nonzero diffusion from the specific-energy gradient |
| P13 | B | Frozen constant diffusivity sinusoid verification | Analytic exponential decay for the isolated transport operator; labelled operator test, not full nonlinear PDE |
| P14 | B | Full nonlinear diffusion manufactured solution | Independently derived source residual; space/time convergence of actual nonlinear coefficient |
| P15 | B | Smooth thermal-plus-SGS pressure wave | Wave speed/energy partition agrees with linearized chosen numerical contract, including SGS compression |
| P16 | B | Large SGS pressure, weak thermal pressure smooth wave | SGS contribution controls stability/wave bounds; no thermal-only timestep assumption |
| P17 | B | Stress energy-work test with paired SGS production | Resolved kinetic change, internal change and SGS change separately agree with compatible discrete accounting |
| P18 | B | Full smooth coupled manufactured solution | Independent residual exercises advection, stress, work, exchange and transport together |
| P19 | B | Temporal refinement with spatial error suppressed | Correct effective order of the selected integrator/source treatment |
| P20 | B | Spatial refinement with small timestep | Correct stencil/reconstruction order in the asymptotic smooth regime |
| P21 | B | Identical transfer integrated over several stage decompositions | No multiplication by stage count, wrong RK weight or reuse of stale stage data |
| P22 | B | Periodic total-energy budget over many short steps | Sum conserved to declared tolerance; applied sources cancel and face fluxes telescope |
| P23 | B | Independent volume-summed output versus online history | Same resolved, thermal, kinetic and SGS inventories within integration/reduction accuracy |
| P24 | B | Known externally forced smooth flow | Actual measured forcing work closes represented-energy budget, including discrete kinetic kick |
| P25 | S | Supported thermal conduction/viscosity combination | Separate transport and viscous/SGS work/heating budgets; no duplicated source |
| P26 | D | Galilean-boosted smooth solution | Agreement after coordinate transformation to discretization accuracy; energy budget transforms correctly |

## R — positivity, shock handling, fallback and timestep

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| R01 | B | Valid state with tiny SGS reservoir and negative production | Backscatter limited by available energy; paired applied transfer and limiter budget |
| R02 | B | Intense forward production with little thermal energy | All reservoirs admissible; compatible kinetic work; no thermal energy borrowed invisibly |
| R03 | B | Deliberately negative reconstructed SGS face state | Defined face-state repair/fallback; shared-face conservation preserved |
| R04 | B | Deliberately negative SGS cell predictor | SGS admissibility recognized before an ordinary scalar floor conceals it |
| R05 | B | Thermal-pressure floor would activate | Recorded thermal/SGS correction and represented-energy impact; no unreported repair |
| R06 | B | Density floor would activate | Defined SGS-specific-energy consequence and explicit correction budget |
| R07 | B | Actual SGS-aware FOFC event | Corrected mass/momentum/energy/SGS fluxes and production use the applied face state/stress |
| R08 | B | Two neighbouring cells trigger fallback on one shared face | One consistent face flux/limiter; no double application or race-dependent result |
| R09 | B | Fallback localized at block/rank interface | Same conservation and admissibility as a single-block/rank case |
| R10 | B | Fallback triggers on different RK stages | Sources and final primitives correspond to the correct stage; no stale-work leakage |
| R11 | B | Sod-type shock tube with SGS disabled | Existing Euler regression unchanged |
| R12 | B | Shock tube with finite SGS energy | Admissible stable evolution, conserved represented energy and dt convergence; do not use Euler solution as exact SGS oracle |
| R13 | B | Strong colliding-flow/shock case | No NaN/negative reservoir; quantified corrections, shock structure and SGS production |
| R14 | B | Strong rarefaction/near-vacuum case | Stable on declared density floor/support range; corrections and restrictions explicit |
| R15 | B | High Mach translation plus small internal energy | Cancellation conditioning understood; no unexplained thermal/SGS loss |
| R16 | B | Low density next to high SGS energy | Safe coefficient evaluation, valid wave/source limit and conservative transfer |
| R17 | B | CFL sweep below and near proposed stable maximum | A characterized stability range; stability at smaller timestep alone does not validate claimed maximum |
| R18 | B | Diffusion-limited mode | Correct dx-squared restriction and global minimum, tested as resolution increases |
| R19 | B | Stress/source-limited mode | Restriction includes actual nonlinear closure behaviour; eddy-viscosity-only bound not presumed sufficient |
| R20 | B | First step with large initial SGS pressure | Initial timestep already respects SGS constraints before the first update |
| R21 | B | Restart with changed local constraints | First resumed timestep valid; old cached dt not blindly reused |
| R22 | B | Minimum constraint occurs on a nonzero MPI rank / different block | Global timestep reduction catches it |
| R23 | B | Zero SGS energy gives no SGS timestep restriction | No divide-by-zero, zero dt, overflow or false finite constraint |
| R24 | B | Multiple source/transport restrictions compete | Correct global minimum and CFL convention; no applied coefficient twice |
| R25 | D | Planar shock with varied grid resolution/orientation | Characterize closure-generated SGS energy versus numerical shock treatment |
| R26 | D | Laminar/coherent shear, rotation and contact | Characterize non-turbulent responses; no blanket expectation that the baseline model recognizes turbulence |
| R27 | D | Repeated limiter activity in a small driven box | Measure changed transfer and correction energy; determine whether the regime should remain supported |
| R28 | B | Selected localized sign/factor/index mutations in a throwaway branch | Independent tests fail; confirms important gates can catch defects rather than mirror code |

## I — parallel execution, restart, output and regression

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| I01 | B | Legacy input suite with SGS off | No regression beyond pre-existing failures; disabled layout/output semantics preserved |
| I02 | B | CPU debug with Kokkos bounds checking | All small coupled tests without out-of-bounds/uninitialized view use |
| I03 | B | CPU release versus debug | Conditioned deterministic results agree; no debug-only initialization masking a defect |
| I04 | B | Several block decompositions for the same domain | Short-run norms/budgets agree; periodic shared-face conservation |
| I05 | S | 1, 2 and 4 MPI ranks where resources permit | State/SGS ghost exchange and global diagnostics correct; no deadlock |
| I06 | S | Rank decomposition in each coordinate direction | Each directional stencil/flux communication path exercised |
| I07 | S | Available threaded backend with different thread counts | No race or incorrect reduction; deterministic short-run tolerance characterized |
| I08 | S | CPU versus available CUDA/HIP backend | Kernel plus short PDE cases agree within precision/reduction tolerance |
| I09 | S | Different supported precision | Admissibility/conservation/order thresholds meaningful for that precision |
| I10 | B | Short unforced run split by restart | Conserved state and subsequent evolution match uninterrupted reference within defined tolerance |
| I11 | B | Forced run split by restart | Forcing modes, random state, update counters and actual work continuity checked separately |
| I12 | S | Restart with changed MPI ranks/pack layout | Correct repartitioned state, counters and forcing policy; explicit rejection if unsupported |
| I13 | B | SGS-aware restart round trip | State values/counts, coefficients, activation and diagnostics metadata retained |
| I14 | B | Legacy SGS-off restart | Backward compatibility retained without needing extra SGS fields |
| I15 | S | SGS-off/on restart mode mismatch | Reviewed rejection/migration; no offset-shift corruption |
| I16 | B | All claimed field-output formats and existing reader | Correct names, units, dimensions and values; distinguish SGS density/specific energy |
| I17 | B | History resolved-plus-SGS budget | Correct active-cell sums, forcing work, transport/exchange/correction terms and time integration |
| I18 | B | Output intervals varied / output disabled | Numerical solution unaffected by diagnostics; no output-triggered recomputation changes |
| I19 | B | Min/max/counter diagnostics with extreme cell on another rank | Global extrema/counters and their volume-weighted energy cost correct |
| I20 | S | Applicable existing hydro/scalar/cooling branch regressions | No unrelated input breakage or passive-index shift |
| I21 | B | C++ style and focused Python checks | Target repository conventions; no new unnecessary dependency or broad unrelated cleanup |
| I22 | B | Restart after a limiter/FOFC event | Applied state/counters consistent; no replay/double-counting of earlier corrections |
| I23 | B | Repeated short runs with fixed input/seed on same build | Reproducible deterministic startup and documented stochastic/reduction policy |
| I24 | S | Clean-clone independent build/run | Same results from recorded inputs and configuration; no hidden local files |

## F — filtered-reference analysis

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| F01 | B | Filter a constant field | Exact normalized constant preservation within roundoff |
| F02 | B | Analytic periodic sinusoid | Known filter transfer function and periodic handling |
| F03 | B | Nonuniform density and velocity with hand-computable averages | Favre weights and exact filtered SGS kinetic energy correct and nonnegative for positive filter |
| F04 | B | Analytic small array for stress covariance | Correct negative-covariance stress convention and trace-energy identity |
| F05 | B | Exact filtered stress versus independent production contraction | Correct density, derivatives, sign and scale in exact SGS production |
| F06 | B | Model evaluated on filtered reference fields | Same declared filter width and variables; no use of unavailable fine gradients in the LES predictor |
| F07 | B | Derivative/filter order and sampling changes | Sensitivity assessed; filter/derivative commutation assumptions explicitly qualified |
| F08 | B | Several filter scales above numerical cutoff | Bias/error/correlation and transfer sign measured by scale; contaminated scale range identified |
| F09 | B | Solenoidal, compressive and mixed reference forcing | Accuracy and conditional biases measured independently for each sampled regime |
| F10 | S | Supersonic plus lower-Mach/adiabatic reference cases | Transferability assessment; no claim outside actually sampled regimes |
| F11 | B | Forward- and backward-transfer subsets | Signed bias, false sign fraction and conditional magnitude errors, not correlation alone |
| F12 | B | Density, SGS Mach, compression and shear subsets | Identify hidden conditional failures even when aggregate correlation is high |
| F13 | B | Separate times/realizations for reference analysis | Correlated sampling handled; uncertainty in reported metrics |
| F14 | B | Reference resolution changed at fixed physical filter | Establish useful reference convergence range; do not label an untested finer run exact |
| F15 | D | Compare plausible explicit filters/effective cutoff mapping | Characterize filter dependence; separate from baseline published-coefficient assessment |
| F16 | B | Held-out data and fixed published coefficients | No tuning leakage; any refit reported as a new calibration study |
| F17 | B | Offline analysis of known synthetic defect in sign or weighting | Filter pipeline detects wrong density weights and stress sign |
| F18 | D | Resolution-degraded reference with alias-prone features | Quantify aliasing and numerical contamination, declare excluded scales |

## T — evolved turbulence campaign

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| T01 | B | Small pilot, model off/on | Stable execution, sane forcing and complete energy ledger; cost estimate before scaling |
| T02 | B | Pseudo-isothermal ideal-gas box | SGS plus resolved energy minus actual forcing work closes; temperature drift explicitly measured |
| T03 | B | Ordinary-gamma energy-evolving box | Dissipative heat and Mach evolution observed with full budget, not hidden by huge internal energy |
| T04 | B | Solenoidal forcing resolution ladder | SGS/resolved statistics and common-filter comparison with uncertainty |
| T05 | B | Compressive forcing resolution ladder | Same diagnostics; shock/fallback/correction dependence quantified |
| T06 | B | Mixed forcing resolution ladder | Same diagnostics; forcing-mixture definition verified |
| T07 | B | At least three initial independent realizations, extended if needed | Confidence estimates account for between-seed and temporal correlation; pilot sample size not assumed sufficient |
| T08 | B | Longer time sampling and changed burn-in window | Stationarity and sampling sensitivity established |
| T09 | B | Forcing work integrated from actual applied kicks | Energy budget agrees with measured work, independent of input dedt |
| T10 | B | Smaller dt on selected boxes | Mean budgets/statistics converge within uncertainty; limiter effects separated |
| T11 | B | SGS seed varied, including documented zero behaviour | Startup and saturated sensitivity reported; no hidden source |
| T12 | B | Driving switched off after a developed state | Decay, dissipation and heat budgets remain closed; stationary assumptions removed |
| T13 | B | Common-scale spectra and structure functions | Resolution/filter/window conventions explicit; reported ranges supported by reference |
| T14 | B | Density PDFs and compressive/solenoidal kinetic fractions | Definitions/weighting and uncertainty consistent between cases |
| T15 | B | Mean SGS energy, production and dissipation versus resolution | Fit scalings with uncertainty; expected stationarity balance includes transport/source/corrections |
| T16 | B | Cellwise/conditional SGS Mach and pressure | Correct thermal sound speed and weighting; identify large-SGS-pressure regime |
| T17 | B | Limiter/floor/FOFC diagnostics through the campaign | Frequency, volume and energy impact reported; no disappearance into residual |
| T18 | S | Second qualified reconstruction/solver | Numerical-method sensitivity assessed without assuming coefficients are optimal for every method |
| T19 | D | Low-Mach and high-Mach stress regimes | Characterize reliability and numerical damping outside core reproduced cases |
| T20 | D | Matched replayed forcing versus independent/seed-matched forcing | Separate pathwise input control from ensemble statistical comparisons |
| T21 | B | Model-off/on versus filtered reference at fixed cost and resolution | Predeclared observable errors; null/negative improvement retained honestly |
| T22 | B | Mass/momentum and mean-flow drift | No unexplained net injection; forcing momentum-removal convention verified |
| T23 | D | SGS pressure/dissipation correlations with resolved indicators | Characterize the model, without interpreting correlation as causal recovery |
| T24 | B | Held-out run that reproduces the final reported comparison | Independent inputs/build manifest and analysis reproduce the claim |

## H — cost, portability and resource qualification

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| H01 | B | Disabled/on matched short runs | Measured kernel and end-to-end overhead, including changed timestep count |
| H02 | B | Memory allocation by grid/pack size | Peak memory and bytes/cell characterized; no unbounded diagnostic/tensor storage |
| H03 | B | Diagnostics on/off and several cadences | Separate compute from output/reduction overhead |
| H04 | S | Available CPU/GPU representative jobs | Synchronization, device transfers and closure work measured with valid timing methodology |
| H05 | D | Strong/weak scaling where allocated resources permit | Communicated state and global reductions quantified; hardware/job conditions recorded |
| H06 | B | Repeated timings after warmup | Median/spread rather than one lucky run; matched output and compiler settings |
| H07 | B | Pilot projected to planned campaign | Wall-time/storage/memory cap before dispatch; revise allocation transparently |
| H08 | B | Error-versus-cost for selected model-off/on cases | Improvement per cost assessed; finest resolution alone not used as cost proxy |

## A — AMR extension after baseline qualification

The [2014 AMR work](https://arxiv.org/html/1309.3996), especially §2 and
Appendix A.2, is a separate numerical foundation for this phase. Its compensation
choices must be reviewed for this model and AthenaK's mesh workflow; simply carrying
SGS as a passive scalar is insufficient. The tests below are proposed adaptations.

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| A01 | A | Uniform fine cells restricted to one coarse cell | Mass/momentum/represented energy conserved; zero resolved-kinetic variance |
| A02 | A | Random positive-density fine cells | Restriction variance identity and nonnegative kinetic difference checked against independent arithmetic |
| A03 | A | Variable density and large uniform bulk velocity | Galilean-invariant kinetic difference computed stably; no subtraction-induced false turbulence |
| A04 | A | Restriction with coherent shear/shock structure | Chosen partition rule applied; mathematical lost kinetic energy not automatically interpreted as turbulence |
| A05 | A | Uniform parent prolonged into children | Parent averages conserved; no invented velocity dispersion or SGS/thermal energy |
| A06 | A | Nonuniform parent stencil prolonged | Fine mass/momentum/represented energy and admissibility meet reviewed transfer contract |
| A07 | A | Insufficient SGS reservoir for interpolated fine kinetic energy | Positivity-constrained reconstruction; explicit fallback and momentum/energy consequences |
| A08 | A | Low thermal energy plus transfer limits | All reservoirs admissible; no hidden thermal compensation |
| A09 | A | Multiple density/SGS contrasts and near-floor child | Limiter does not privilege one child incorrectly or violate weighted parent invariants |
| A10 | A | Restrict/prolong repeatedly without physical advance | No secular represented-energy drift; thermal/SGS partition effects measured, not assumed reversible |
| A11 | A | Coarse ghost fill versus permanent refinement | Correct filter/stress state without double compensation or double-counted physical volume |
| A12 | A | Constant flow crosses fixed coarse/fine interface | No mesh-interface force or SGS generation |
| A13 | A | Advected SGS pulse crossing interface | Conservative scalar transport and proper filter/partition response |
| A14 | A | SGS diffusion across interface | Correct face flux, coefficient scale, reflux and global conservation |
| A15 | A | Stress/pressure wave crossing interface | Compatible momentum/work/production, reflection quantified against reference |
| A16 | A | FOFC at a coarse/fine interface | Refluxed flux and local sources remain compatible and admissible |
| A17 | A | Refinement-level timestep constraints | Finest dx/filter and stiffest state control global step; audit subcycling only if present/supported |
| A18 | A | Insert nested grids into developed turbulence | Controlled energy compensation and transients; hierarchy and forcing histories recorded |
| A19 | A | Forced refine/coarsen oscillation | Isolate remapping from physics; quantify correction and partition accumulation |
| A20 | A | Threshold dithering/refinement hysteresis cases | No runaway grid thrashing or hidden energy source; deterministic mesh-control behaviour |
| A21 | A | Feature moving across levels and domain periodic boundary | Correct remap, ghosts, stress and transport with no duplicate/inaccessible leaf state |
| A22 | A | Regrid then load balance to different ranks | All state, indices, coefficients and correction counters survive migration |
| A23 | A | Restart immediately before/after regridding | Hierarchy and SGS partition reproduce the uninterrupted supported trajectory |
| A24 | A | MPI ranks change on AMR restart | Reviewed repartition policy preserves invariants; rejected if unsupported |
| A25 | A | Leaf-cell reductions and covered coarse cells | Physical volume/mass counted once; history agrees with independent leaf inventory |
| A26 | A | Common-filter analysis on adaptive hierarchy | Conservative comparison at declared physical scale; no plotting interpolation treated as new physics |
| A27 | A | Static versus dynamic hierarchy with matched selected statistics | Separate interface bias from remap bias and physical-resolution differences |
| A28 | A | AMR versus matched uniform ladder, multiple realizations | Error/energy/SGS statistics with uncertainty; effective resolution defined |
| A29 | A | Different refinement criteria and active volume fractions | Selection bias and refined coverage measured; no criterion assumed universally optimal |
| A30 | A | AMR cost accounting and finest-grid saturation | Cell updates, occupancy, global dt, communication and wall time; no savings claim if most volume refines |
| A31 | A | Available GPU/MPI AMR execution | Reallocation, pack extents, remap and device kernels remain valid |
| A32 | A | Explicitly varied transfer exponent/partition assumption | Sensitivity report; does not count a freely tuned exponent as faithful baseline reproduction |

## X — future physics, only after a new contract

| ID | Class | Setup | Oracle / required outcome |
|---|---|---|---|
| X01 | X | Radiative cooling plus SGS heating | Independent thermal budget, source-order/stiffness tests and declared cooling closure |
| X02 | X | Multiphase cooling/interface problem | Phase mass, mixing and heating convergence; momentum SGS alone does not define thermal mixing |
| X03 | X | MHD closure | Magnetic stress/EMF and kinetic–magnetic budgets with divergence constraint tests |
| X04 | X | SGS scalar/enthalpy transport | Separate turbulent transport closure and conservation/positivity tests |
| X05 | X | Self-gravity/unresolved gravitational effects | New filtered source closure and total-energy/work contract |
| X06 | X | Stellar/feedback SGS injection | Explicit external energy/momentum source and reservoir-partition tests |
| X07 | X | Learned/dynamic coefficients or explicit-filter alteration | Calibration/held-out accuracy, stability and cost evaluated as a distinct model |
| X08 | X | Synthetic line broadening from SGS motions | Declared velocity distribution, tracer weighting and observable validation |

## Concrete suite construction

Build tests in increasing cost and reuse existing AthenaK infrastructure:

1. One closure/reference suite covers K families using hand cases and parameterized
   random tensors, rather than one executable per row.
2. A focused verification pgen/input set covers lifecycle and P/R tests with term
   isolation or manufactured forcing explicitly labelled as test-only mechanisms.
3. Existing regression scripts cover configuration, disabled-mode and scalar compatibility.
4. Integration scripts cover block/MPI/restart/output/backend variants.
5. A batch protocol covers F/T/H scientific campaigns and collects manifests/metrics.
6. A separate AMR suite is built for A0–A2 when that extension begins.

Some term-isolation tests need a test harness, not public knobs to disable every
model term. Coefficient endpoints are valid diagnostics only if the parameter
contract permits them; otherwise use internal verification functions.

## Scheduling tiers and release coverage

| Tier | Trigger | Coverage | Evidence |
|---|---|---|---|
| Fast | Relevant packet/commit | Kernel checks, selected C/L/P cases and disabled-path smoke | Seconds/minutes where feasible; failing assertions and metrics |
| Integration | Coupled milestone / affected path change | Core C/L/P/R/I, selected MPI/restart/backend | Build/run logs, norms, budgets, corrections |
| Qualification | Baseline release candidate | All mandatory B families and claimed S paths, targeted method coverage | Coverage map and independently reproduced run bundle |
| Science campaign | M7–M8, frozen protocol | F/T/H with independent realizations/references | Raw runs, uncertainty and comparison report |
| Refinement qualification | A0–A2 | A catalogue plus relevant baseline regressions | Separate remap/interface/dynamic hierarchy report |

Before selecting concrete instances, classify which rows verify algebra,
discretization, state lifecycle, physical modelling or statistical benefit. Passing
one level does not substitute for another. Repeat tests only after relevant changes,
new failures or unresolved concerns; retain unchanged expensive results with exact
commit/dependency provenance and state why they remain applicable.

## Required interactions, beyond isolated feature tests

Pairwise selection reduces redundant combinations, but these intersections are
mandatory wherever the constituent paths are supported. Add concrete instance
names under the existing test IDs rather than claiming a separate isolated pass
covers their interaction.

| Interaction | Relevant families | Specific defect to expose |
|---|---|---|
| SGS + FOFC + multiple blocks/MPI + restart | R07–R10, I05, I22 | Shared-face inconsistency or correction replay |
| SGS + multiple true scalars + restart/output | L02–L03, I13–I16 | SGS slot aliases metallicity/dust or shifts serialized offsets |
| High SGS pressure + low thermal energy + near-limit dt | P16, R15–R20 | Thermal-only wave bound or wrong work partition |
| Backscatter + small reservoir + multistage update | R01, R10, P21 | Stage-weight positivity failure or unmatched exchange |
| Nonlinear diffusion + strong density contrast + boundary crossing | P11–P14, L11–L14, R18 | Diffusion of the wrong energy variable or ghost coefficient |
| Forcing + SGS + restart + rank change | I11–I12, T09 | Random-state mismatch, forcing discontinuity or double-counted work |
| GPU + MPI + partially filled packs | K24, L13, I05, I08 | Padded memory access, synchronization or reduction defect |
| High-order reconstruction + minimal allowed ghosts + FOFC | C12, R03, R07, I02 | Corner/extended-stencil out-of-bounds access |
| AMR + SGS + FOFC + reflux | A14–A17, R07–R10 | Corrected face work inconsistent with coarse/fine flux |
| AMR + regrid + migration + restart + GPU/MPI | A21–A24, A31 | State loss/reinitialization during topology and backend changes |
| Diagnostic counters + restart + rank reduction | I17–I19, I22, A25 | Discontinuous accumulators, wrong active-cell sums or hidden corrections |

Also include mass and momentum inventories in every closed periodic conservation
case, even when the principal assertion is energy. The budget oracle must account
for applied boundary fluxes in supported open-boundary cases.

## Result ledger template

Create the actual ledger during implementation, not a fictitious table of passes now:

| Test ID / instance | Gate | Commit / input hash | Build / backend | Status | Metric / threshold | Evidence / limitation |
|---|---|---|---|---|---|---|
| Example only | M2 | To be recorded | To be recorded | not run | To be frozen | No execution in this planning session |

Allowed statuses: `not run`, `pass`, `fail`, `blocked`, `not applicable`.
Each blocked/not-applicable row names the missing resource or explicit support
restriction. Every failure links raw evidence and a later fix/retest; preserve
the original failure rather than overwriting it with a green result.
