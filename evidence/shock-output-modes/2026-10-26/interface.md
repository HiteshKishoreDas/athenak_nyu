# Frozen binary encounter interface (implementation target)

Inputs: existing file_type=crv (curve nodes) or file_type=prtcl_thermo_history
(particles), mode=encounters; default mode=snapshot. Exactly one active history
block per carrier, enabled by native dcycle!=0 or dt>0 semantics. Event blocks
contain no thermo/geometry snapshots. Keep global shock_tracking/enabled as
sampler/detector configuration, independent of explicit history activation.

Files: shock_events/<basename>.<id>.cycle<G>.rank<R>of<N>.thp and
shock_events/curves/<basename>.<id>.cycle<G>.rank<R>of<N>.crv. Terminal views use
.suspended<G> in place of .cycle<G>. No completed-event CSV journals. Optional
qualification trace CSVs remain explicitly opt-in. Historical readers remain.

All scalars packed individually, no native struct padding, declared byte order.
Magic16: ATHK_EVENTS followed by5 zero bytes. Then uint32 version1, endian marker
0x01020304, Real byte width, carrier(0 particles/1curve), curve_id_bits (0 for
particles), kind (0 completed batch,1censored terminal view,2pending checkpoint).
Then uint64 generation, previous_flush_generation (UINT64_MAX before first flush),
previous_file_digest (FNV whole predecessor including its footer), prefix_digest (FNV stream digest before this batch), count;
Real time; int32 rank/nranks; uint32-length UTF8 strings lineage, contract, stem.
Final uint64 FNV-1a checksum of every preceding byte.

Each event: uint64 emission_generation; int32 signed identity; uint64 sequence;
uint32 status/flags; 5 observations (pre,first,last,post,representative); Real birth
(-1 for particles), int32 first_level/last_level (-1 for particles), Real
first_spacing[3]/last_spacing[3] (zero for particles).
Observation: Real time; uint64 generation; uint32 valid/inside; Real x[3]; int64
image[3]; Real rho,pressure,temperature,mach_pressure,mach_temperature,normal[3],
upstream[3],downstream[3]. This is the existing observation/event contract, typed.
Curve IDs interpret int32 exactly; no float conversion. Curve contract prefix
continues curve-node-history-v4 <bits>. Particle identities remain nonnegative.

Completed batches cover (previous_flush_generation,generation], with initial
coverage starting0. Per-file final checksum and previous_file_digest detect
missing/modified batches; prefix_digest binds entire stream. After a flush,
rolling digest = FNV(prefix_digest, whole batch bytes including checksum).
Pending checkpoint packets save unflushed records, last-flush generation/hash and
current rolling digest without clearing the queue. Terminal censored views do
not consume identities or advance the committed stream. Output cadence controls
flush only; current every-step physical encounter rules remain unchanged.

Reader: parse into same dictionaries/identities as historical shock_history and
curve_shock_history APIs. Emission_generation is original completion clock, not
flush clock. Validate file lengths/checksums/types/flags/brackets, complete rank
sets, stream metadata, previous links/digests and duplicate identities. Censored
views remain separate. A through_cycle query between flushes must not silently
claim coverage: reject unless committed batches cover that requested clock.
