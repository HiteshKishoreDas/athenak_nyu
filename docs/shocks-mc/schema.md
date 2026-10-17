# Event and restart contract

All quantities use Athena code units. Event identities are the tuple
(lineage, unsigned64 particle_id, unsigned64 sequence). Existing native tags are
signed32 and are checked before casting/incrementing; event output never converts
IDs through floating point. No exact physical crossing time is inferred.

## Numerical detector

Pressure=(gamma−1)*primitive IEN, temperature=P/rho, entropy=logP−gamma log(rho).
Use anisotropic centered derivatives on active axes and real synchronized ghosts.
Temperature gradient supplies normals; largest absolute component chooses nearest
signed axis, ties x/y/z; endpoints are separated by width cells. Pressure orders
upstream/downstream. Strict compression/gradient/alignment/jump comparators,
inclusive upstream floors, density saturation and inherited25percent temperature
consistency filter follow the frozen ShockIt source. Pressure Mach is uncapped;
temperature analytic inversion has its own cap. Either estimator meeting the
inclusive threshold accepts unreduced region membership. Invalid touched states or
arithmetic mark observations invalid (documented deviation from global oracle abort).
Physical-outflow extrapolation is not a valid sampled upstream/downstream endpoint.

## Encounter states and evidence

Outside: preserve latest valid prior observation. Enter active on valid region
membership if the previous rho/P state is nearer simultaneous upstream than
 downstream in absolute log-distance; otherwise retain suppressed-reentry evidence.
No prior valid observation means left-censored. Representative is first valid
inside sample; last tracks the most recent inside evidence. Normal reversal is
flagged ambiguous. First valid outside enters pending_exit; two consecutive valid
outside observations emit and rearm. Earlier reentry merges and flags flicker.
Invalid observations cannot close an episode. Same-generation retry is a no-op.
Births initialize only new state; source-cell MC movement uses the original RNG.

Domain exit emits a partial record before deletion, with actual timestep endpoint
and out-of-domain position in invalid post evidence. End-of-run suspended records
are right-censored views, keep their identity and may resume; they do not consume
another sequence or duplicate the canonical committed-event stream.

Every event includes `pre`, `first`, `last`, `post`, `representative` observations.
For valid pre/post, entry mask bracket is [pre.time,first.time] and exit mask bracket
[last.time,post.time]. These bracket sampled mask evidence, not a proven material
physical-shock crossing. Fixed sampler endpoints are not certified plateaus.
MC jumps or stochastic returns cannot be replaced by interpolated exact times.
All records are cadence_limited; missed passages remain a qualification limit.

## CSV schema1

`manifest.rankR`: three lines: `AKSHOCK-EVENTS 1`, lineage, opaque resolved contract.
Every generation has `cycleG.rankR.csv`, including empty generations. Header:
`# AKSHOCK-EVENTS 1 lineage G count`; second line contract; footer `# END count`.
Rows have124columns: particleID,sequence,status,flags + five24column observations.
Each observation: time,generation,valid,inside,x1,x2,x3,image1,image2,image3,
rho,P,T,MachP,MachT,normal1..3,upstream(rho,P,T),downstream(rho,P,T).
Images are signed64 global periodic wraps. NaN is permitted only for explicitly
invalid evidence or invalid Mach estimates, never positive required state data.

Status0 means completed numerical episode; status1 partial/censored.

| Flag value | Meaning |
| --- | --- |
|1|left censored: no valid pre-entry sample|
|2|invalid observation during tracking|
|4|pending exit merged with an early reentry|
|8|no certified upstream/downstream plateau|
|16|endpoint cadence limitation|
|32|temperature Mach invalid while region accepted through another estimator|
|64|right censored|
|128|domain exit|
|256|normal reversal/direction ambiguity|
|512|downstream reentry suppressed|

`vis/python/shock_history.py` validates integers, generation continuity, metadata,
counts/footer and canonical identity uniqueness. `read_history()` returns committed
records. `read_suspended()` returns latest censored views separately. Optional
trace files have28columns (ID,sequence,phase,flags+observation).
Optional transfer files carry version1 header with lineage,generation,count,start/
end time,RNG seed;17columns: ID,sourceGID,i,j,k,start density,sixoutgoingprobabilities,
x,y,z,draw,direction0..6. Existing .thp is legacy thermo output, not an event log.

## Persistence scope and failure model

Native MC `ATHKPRTCLMC` v1 already saves7real and5integer rows, schedule progress and
next tag. `AKSHOCK2` sidecar v2 saves complete trivially copyable encounter states,
rank/count/partition/native ABI/time/generation, bound source/config contract,
diagnostic counters and prefix watermark. FNV64 checks native checkpoint bytes,
sidecar bytes and the chained complete event prefix; it detects accidental damage,
not hostile collisions. Same exact implementation/ABI/rank layout is required.
The old velocity-tracer sidecar format is incompatible and is rejected.

Per-generation chunks are sorted and committed using exclusive temporary files,
write/fsync/close, non-overwriting hard link, cleanup and directory fsync. Existing
chunks are accepted only for byte-identical replay. Sidecar commits only after
native payload/metadata writes and close checks, native file flush and checksum.
Checkpoint preflight rejects an existing fluid/sidecar pair before truncation;
MPI peers collectively test conflicts. Missing/torn/incompatible files and missing
or changed prefix chunks reject; replay requires retained prefix through checkpoint.
Pending domain-exit events drain before checkpoint; active/pending state persists.

This is same-layout checkpoint replay with explicit failures, not an all-rank atomic
transaction or a claim about every filesystem/power-loss fault. Storage short-write/
close checks were source-reviewed; actual storage-failure injection remains pending.
RUSTY/Frontier filesystem durability and true cross-rank recovery need target tests.


Current state stores four observations (previous/first/last/post). Event representative
is derived from first when emitting, retaining all five public observation groups.
Three-review source3de18d97 reduces local double EncounterState1000to816bytes;
Event remains944bytes and native MC header48bytes, with unchanged field offsets.
Earlier1000-byte sidecars are rejected, even for zero particles. Replay always
requires the matching executable; no cross-version checkpoint conversion is supplied.
