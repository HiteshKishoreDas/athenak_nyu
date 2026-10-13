# SGS lifecycle/integration audit

Read-only packet, 2026-10-08. Checkout `/Users/hitesh/hitesh/git/athenak_nyu`,
verified HEAD `07c529202befbeb6ef5e694d95286929618ff7c1`. No builds or source edits.
Existing changes at inspection: modified `.gitignore`, untracked `AGENTS.md`.
Read repository AGENTS, user attachment, wiki subject/development-plan/runbook/test-matrix.
Applied ponytail's reuse/minimum-change guidance. Source line numbers below are this HEAD.

## Proposed layout contract

- Keep `nhydro=5`, `nscalars` = true user passive count, unchanged.
- Add `nsgs = enabled ? 1 : 0`, `nvars = nhydro+nscalars+nsgs` and
  `isgs = enabled ? nhydro+nscalars : -1` to Hydro.
- Conserved `u0(isgs)=K` is SGS **energy density**. Primitive `w0(isgs)=K/rho`
  is specific SGS energy. Existing scalar `s_i` remains at `nhydro+i`.
- `u0(IEN)` remains resolved total energy (internal + resolved kinetic), not K.
  Represented energy is `u0(IEN)+u0(isgs)`. This preserves existing `ener`, `eint`,
  temperature and coarsened-moment meanings. Dissipation subtracts K/adds IEN;
  production adds K/subtracts IEN using the same limited exchange.
- When off, nvars equals the old count, all cell/face/communication buffers retain
  exact old extents, no SGS arrays/kernels/diagnostics are allocated or launched.
  Merely adding always-allocated one-cell Kokkos Views is not zero allocation cost.
- Prefer one total-count member rather than reinterpreting nscalars or global
  variable registry abstractions. Every count site below should use that member.

## Exact count and allocation sites

Paths are relative to `/Users/hitesh/hitesh/git/athenak_nyu`.

| File:line | Action |
|---|---|
| `src/hydro/hydro.hpp:78-79` | Declare true passive count separately from total count and SGS slot. |
| `src/hydro/hydro.cpp:72` | Parse enable, counts and coefficients before allocation; reject negative nscalars and incompatible configurations. |
| `src/hydro/hydro.cpp:104-105` | u0/w0 second extent -> nvars. |
| `src/hydro/hydro.cpp:114-115` | coarse_u0/coarse_w0 count -> nvars; baseline rejects multilevel, so no new AMR path claimed. |
| `src/hydro/hydro.cpp:120` | boundary buffer count -> nvars. |
| `src/hydro/hydro.cpp:124-125` | orbital/shearing count -> nvars for correctness of generic total-count convention; reject shearing_box while SGS on. |
| `src/hydro/hydro.cpp:275-278` | u1/all three uflx counts -> nvars. |
| `src/hydro/hydro.cpp:283` | utest is currently nhydro only. SGS-on FOFC needs at least nhydro plus the SGS slot; simplest total extent nvars when enabled, old nhydro when off. |
| `src/hydro/hydro_tasks.cpp:93,99` | InitRecv/InitFluxRecv MPI counts -> nvars. |
| `src/hydro/hydro_tasks.cpp:142` | CopyCons rk4 average loop -> nvars (rk4 remains rejected initially if unsupported). Other copy is deep_copy and already transports whole extent. |
| `src/hydro/hydro_fluxes.cpp:46` | reconstruction scratch count -> nvars. |
| `src/hydro/hydro_update.cpp:36` | RK flux update count -> nvars. |
| `src/outputs/restart.cpp:70,152` | host allocation and serialization count -> nvars. |
| `src/pgen/pgen.cpp:552` | restart input count -> nvars; feeds both partitioned and shared readers. |
| `src/mesh/load_balance.cpp:148-149,412-413,553,864` | four count sites -> nvars; important for offsets after hydro even though AMR baseline rejected. |
| `src/pgen/gr_torus.cpp:1761` | currently total hydro count; SGS+GR rejected, therefore unreachable. Do not alter unrelated pgen while implementing baseline. |

The remaining `nhydro+nscalars` occurrences in passive initialization/output are
intentional passive bounds and must **not** become nvars: `slotted_cyl.cpp:117`,
`tests/advection.cpp:144`, `cgm_cooling_flow_no_disk.cpp:581`, and
`outputs/basetype_output.cpp:309,340`. `tests/lw_implode.cpp:78,82` initializes only
the first true scalar; its scalar identities remain unchanged.

## EOS and floors

- `src/eos/ideal_hyd.cpp:35-36,96-103`: the true scalar loop clips any negative
  conserved scalar to zero silently, then divides by density. Leave its bounds
  `nhydro..nhydro+nscalars`; add dedicated SGS conversion after it. SGS must never
  enter that clipping path.
- `ideal_hyd.cpp:70-75` only checks density/thermal floors in `only_testfloors`;
  add SGS admissibility checks there when enabled, before any clipping. NaN/Inf,
  negative K and represented/thermal partition should be deliberate failures or
  trigger the declared conservative fallback, not disappear into max(0,K).
- `ideal_hyd.cpp:127-128,152-155`: dedicated `K=rho*k` primitive-to-conserved entry.
- `src/eos/ideal_c2p_hyd.hpp:31-64`: density/thermal floors can alter represented
  energy or kinetic partition independently of SGS. Existing counters count
  events, not correction energy. SGS diagnostics must track actual IEN/K correction
  energies and density correction consequences. Do not change the off path.
- `ideal_c2p_hyd.hpp:59-64` entropy floor changes w.e without updating u.e, unlike
  pressure/temperature floors. This is pre-existing and needs explicit audit in
  SGS enabled admissibility; default sfloor is the likely mitigation but not proof.
- SR/GR/isothermal EOS files also have scalar-only conversions; reject those
  physics, so do not make the SGS state silently reachable through them.
- `src/bvals/prolong_prims.cpp:53-54,208-209` implements its own hydro conversion
  with scalar bounds; coarse SGS conversion/energy compensation is unimplemented.
  Reject `mesh.multilevel` (both static refinement and AMR), not just adaptive.

## Reconstruction/advection

- DC/PLM infer variable count from `q.extent_int(1)`:
  `src/reconstruct/dc.hpp:23,42,61`, `plm.hpp:48,66,84`; appended slot is
  automatically reconstructed.
- PPM (`ppm.hpp:195,243,291`) and WENOZ (`wenoz.hpp:93,129,165`) also use total
  extent, but their floors only apply to IDN/IEN, not SGS. SGS face positivity
  must be handled explicitly, or reject these methods until qualified.
- `hydro_fluxes.cpp:133-142,234-243,329-338` transports all entries after nhydro
  with mass-flux upwinding. Appended primitive specific energy can reuse this
  advection formula. The routine currently computes scalar flux only on active
  faces, whereas hydro FOFC fluxes are computed on extended faces. New SGS FOFC
  must fill every face used in its extended predictor, or restrict the predictor.
- Hydro Riemann solver structs contain only five hydro variables. SGS pressure
  is not included automatically just because reconstruction carries its slot;
  that is a separate numerical-contract integration decision.

## FOFC and stage ordering: highest-risk integration point

- `src/hydro/hydro_fofc.cpp:50,68-77` trial update currently advances only nhydro.
  To assess K, extend predictor to include the SGS slot and paired stage sources.
- `hydro_fofc.cpp:82` calls only EOS floor detection. Negative K cannot trigger
  correction without the dedicated EOS/SGS check above.
- `hydro_fofc.cpp:155-159,191-195,236-240,273-277,319-323,356-360`
  replaces only five hydro entries on six surrounding faces. It leaves passive
  fluxes at their original values despite changing mass flux. For SGS-on, recompute
  K advection using corrected mass flux, replace SGS stresses/work/diffusion as
  one face operation, and recompute matched production using final face values.
  Passive flux consistency with corrected mass flux should also be checked in
  SGS+multiple scalar tests; preserve legacy-off behavior unless fixing separately.
- `hydro_tasks.cpp:184-194` adds diffusion fluxes before FOFC; FOFC can discard
  momentum/energy diffusion without discarding other fields. Do not add SGS fluxes
  there and assume existing FOFC preserves their paired effects.
- `hydro_tasks.cpp:55-75`: Copy -> Fluxes(+FOFC) -> SendFlux/RecvFlux -> RKUpdate
  -> source terms -> ghost communication -> BC -> ConsToPrim -> dt.
  Flux/source SGS coupling must use final shared face flux and its stage weight.
  Sources added at `HydroSrcTerms:238-255` are not in existing FOFC predictor.
- A cell-by-cell post-update clamp or per-cell flux replacement can give different
  decisions on a shared face and break MPI conservation. SGS-aware fallback must
  have a single deterministic/shared face result, and ledger its paired exchange.
- RK weights come from Driver, especially `driver.cpp:99-111` RK2. Freeze RK2
  initially; derive positivity for exact stage weights, not full dt alone.

## Boundaries and communication

- `bvals/bvals.cpp:111-116`: u_in is nvar x 6; increasing Hydro InitializeBuffers
  count carries the SGS slot to inflow storage, but pgen must populate it.
- `bvals/bvals_cc.cpp:46,304` packs/unpacks all variables using array extent.
  Correct allocations and InitRecv count suffice for same-level periodic state.
- `bvals/physics/hydro_bcs.cpp:33` loops over full array extent. Reflect parity
  is even for appended K, and outflow extrapolates it. These are plausible state
  semantics but SGS stress/diffusion boundary work still requires qualification.
- Inflow reads u_in(slot,face) (`hydro_bcs.cpp:59-62,98-101` and analogues), user
  BCs can ignore new slot, vacuum sets density/K to zero, diode changes momentum.
  Reject these initially unless explicitly implemented/tested with energy budget.
- Minimal initial support: all physical faces periodic. This includes multi-block
  uniform grids; do not confuse multiple blocks with multiple refinement levels.
- For centered multidimensional gradients and face gradients, map required corner
  ghost depth against PLM/FOFC extents before freezing minimum ng.
- `bvals/prolongation.cpp`, `mesh/mesh_refinement.cpp:1200-1232`, and AMR
  pack/unpack kernels use variable extents; this transports K numerically but
  does not supply the required filter-change energy compensation.

## Restart contract and missing protection

- Existing parameter dump (`restart.cpp:190-193`) stores GetOrAdd parameters;
  adding enable/coefficients/version with GetOrAdd naturally persists them.
- Existing payload size (`restart.cpp:285-307`) and reader size check
  (`pgen.cpp:659-686`) use total variable byte count. Updating three count sites
  above serializes appended K without a new binary record and leaves old off
  restarts readable.
- Critical: `src/main.cpp:322,326-333` loads restart params, then -i and CLI may
  override them **before** Hydro construction. Byte-count check alone cannot
  reject every incompatible layout: changing nscalars while toggling SGS can
  yield the same total count and silently reinterpret a passive slot as K.
- Freeze and compare restart metadata before overlays: SGS enable/model version,
  nscalars, EOS/nhydro, SGS coefficients/filter definition, and any support-critical
  numerical parameters. A small restart-specific snapshot/check in main is simpler
  than adding binary-header machinery that breaks legacy layout. Legacy missing
  SGS metadata means SGS off. Reject off->on/on->off, changed scalar count,
  unsupported model version, and changed frozen model coefficients; no migration.
- Same-layout SGS roundtrip should restore u0 including K and ghost state, then
  pgen restart hook must enroll callbacks without reinitialization.
- `pgen.cpp:1054-1056` calls UserProblem(pin,true) after state load.
- Cumulative diagnostics need explicit restart persistence, separate from state.
  Existing event counters are not automatically SGS cumulative energy budgets.
- Forcing RNG already serializes in restart and reader `pgen.cpp:603-623`; forcing
  arrays serialize separately. Preserve those offsets and verify restart continuity.

## Output and pgen hooks

- `src/outputs/basetype_output.cpp:303-360` scalar output loops/names remain true
  passives (`r_0/s_0`, etc.); preserve them unchanged.
- Add explicit SGS output choices in `outputs/outputs.hpp` (NOUTPUT_CHOICES must
  match table count) and mapping in `basetype_output.cpp` with enable validation.
  Suggested `hydro_u_sgs` -> `sgs_energy` and `hydro_w_sgs` -> `sgs_specific`.
- **Existing `hydro_sgs` name collision:** `basetype_output.cpp:778-787` and
  `derived_variables.cpp:662-709` already output 14 raw moments for postprocessing
  filtering. Do not rename or repurpose these fields as the new closure tensor.
  Retaining IEN as resolved energy also preserves their internal-energy formula.
- `history.cpp:78-137` hydro history measures resolved energy/mass/momentum/KE;
  extend enabled-only diagnostics with SGS K and represented total, production,
  dissipation, actual forcing work, floors/fallback adjustments and thermal partition.
  Report signed cumulative exchanges and instantaneous rates distinctly.
- `src/pgen/turb.cpp:29` returns immediately on restart, before enrolling history
  at :41. Its existing history `:194-196` dereferences pmhd and therefore is not
  safe for HD-only user_hist. Enroll callbacks before the return and provide an
  HD-safe branch/dedicated SGS pgen history instead of invoking MHD-only history.
- Turb pgen initializes hydro five entries at `:61-72`, no true scalars or K.
  Generic SGS initialization after ordinary pgen should set K=rho*k_init only
  on fresh runs; restart must never seed again. If verification pgen supplies its
  own K profile, the initialization contract must prevent overwriting it.
- Existing pgen scalar identities are hard-coded throughout CGM metallicity/dust
  and turbulent_scalars generators (`nhydro`, `nhydro+1`, etc.); nscalars must stay
  true count. `srcterms.cpp:299` and `srcterms_newdt.cpp:157` read first scalar as Z.
  Cooling combinations are rejected while SGS enabled but SGS off must be intact.

## Support rejection sites and initial policy

Reject at Hydro construction before arrays/expensive cooling tables if possible.
Use input block presence rather than ppack pointers for later-built modules:
`meshblock_pack.cpp:109-137,159-176,190-228` constructs Hydro before many modules.

Mandatory reject while SGS enabled: non-ideal EOS; SR/GR/dynamical coordinates;
mesh.multilevel/static/adaptive refinement; MHD/ion-neutral/radiation/ADM/Z4c;
cooling flags `ism_cooling`, `cgm_cooling`, `rel_cooling`; shearing_box/orbital;
kinematic/static evolution. Also reject particles or other coupled physics until
qualified rather than implying published HD baseline supports them.

Initial qualification recommendation: RK2, PLM, chosen HLLC/LLF-compatible pressure
treatment, periodic uniform Cartesian mesh, double precision, true passive scalars,
CPU serial. Mark MPI/GPU conditional until executed. Reject unsupported integrators,
Riemann solvers/reconstruction, viscosity/conduction, constant acceleration/user
source/work-in-loop hooks unless a reviewed test-only interface or qualified budget
contract enables them. Forcing is within scope but must qualify task ordering/work.
Do not reject strict isothermal as a published-model judgment: reject it because
the requested baseline requires a thermal energy equation/heating verification.

## Timestep and initialization sites

- `hydro_newdt.cpp:29-32` computes only final-stage estimates; Driver Initialize
  calls it with nexp_stages on both fresh/restart (`driver.cpp:309-323`). SGS bound
  integrated into dtnew therefore covers first/restart step automatically after
  ghost/primitive initialization (`driver.cpp:597-614`).
- Newtonian hydro presently bounds waves using thermal EOS sound speed only;
  extend enabled path with derived SGS pressure/stress/diffusion/source bounds.
- `mesh.cpp:594-605` applies cfl_no to hydro and diffusion/source dtnew; decide
  whether SGS limits include CFL factors themselves to avoid double application.
  Existing global MPI reduction then propagates minimum; no new global dt path
  should be necessary if SGS contributes to Hydro dtnew.

## Concrete exit checks for lifecycle owner

No tests were executed by this audit. At implementation exit select finite small
instances from C/L/R/I guidance, recording frozen thresholds and exact inputs:

1. SGS-off existing scalar regression and legacy restart (same extents/output names).
2. Multiple true scalar fields + appended nonuniform K at varying density: EOS
   roundtrip and passive output identities, explicit K/rho expected values.
3. Uniform positive/zero K fresh run, restart through same pgen without reseeding.
4. Periodic K advection across every axis, one vs many blocks, ghost edges/corners.
5. Negative K predictor, strong SGS pressure/low thermal state and FOFC: shared
   corrected mass/K/stress/work consistency, no scalar-floor concealment.
6. SGS restart on/on, legacy off/off, off/on, on/off, nscalars collision with equal
   total count, changed coefficients, broken payload, changed rank layout where supported.
7. Existing hydro_sgs filtering output byte/field meanings unchanged; new fields
   recover energy density/specific energy independently.
8. Every rejected configuration fails at initialization with the rejected option
   named; model off remains accepted for previously valid legacy inputs.

## Memory used

Memory registry `/Users/hitesh/.codex/memories/MEMORY.md:258-286` reinforced that
AMR state transport alone is not energy compensation and that earlier notes were
unimplemented framing. Live checkout/wiki evidence above supersedes stale context.
Associated rollout: `01a113a0-cf1a-7ab3-bed7-d1da0ba79c75`.
