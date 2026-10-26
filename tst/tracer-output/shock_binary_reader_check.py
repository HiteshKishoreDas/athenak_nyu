"""Stdlib reader checks against native batches; retain every mutated fixture.

Run after the native suite, optionally comparing its unchanged physics to an
earlier per-step CSV run. Synthetic domain exits test serialization semantics,
not live outflow physics. No builds or simulations are launched here.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import struct
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'vis/python'))
import shock_binary as binary
import shock_history as particles
import curve_shock_history as curves


def layout(data):
    endian = '<' if data[20:24] == b'\x04\x03\x02\x01' else '>'
    width, = struct.unpack_from(endian+'I', data, 24)
    offset = 88+width
    for _ in range(3):
        count, = struct.unpack_from(endian+'I', data, offset)
        offset += 4+count
    return endian, width, offset


def repair(data):
    endian = layout(data)[0]
    data[-8:] = struct.pack(endian+'Q', binary.digest(data[:-8]))
    return data


def normalized(events):
    # CSV/binary source contracts can differ. Every physical/provenance field,
    # including original emission clocks, must match; NaNs compare canonically.
    return json.dumps([{k:v for k, v in event.items() if k not in ('lineage', 'identity')}
                       for event in events], sort_keys=True)


def domain_exit_packet(event, header, kind, flags):
    """Independent typed encoder for the committed/pending exit regression."""
    data = bytearray(binary.MAGIC+struct.pack('<6I5Qd2i', 1, 0x01020304, 8, 0, 0,
        kind, header['generation'], binary.MASK, 0, binary.OFFSET, 1, header['time'], 0, 1))
    for text in (header['lineage'], header['contract'], 'synthetic.domain_exit'):
        encoded = text.encode('utf-8')
        data += struct.pack('<I', len(encoded))+encoded
    data += struct.pack('<QiQ2I', header['generation'], event['particle_id'],
                        event['sequence'], 1, flags)
    for key in ('pre', 'first', 'last', 'post', 'representative'):
        o = event[key]
        data += struct.pack('<dQ2I3d3q14d', o['time'], o['generation'], o['valid'],
            o['inside'], *o['x'], *o['image'], o['rho'], o['pressure'], o['temperature'],
            o['mach_pressure'], o['mach_temperature'], *o['normal'], *o['upstream'], *o['downstream'])
    data += struct.pack('<d2i6d', -1, -1, -1, 0, 0, 0, 0, 0, 0)
    return bytes(data)+struct.pack('<Q', binary.digest(data))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', type=Path, required=True,
                        help='Native complete run, or suite containing complete/')
    parser.add_argument('--reference', type=Path,
                        help='Optional prior CSV run or its shock_events directory')
    parser.add_argument('--output', type=Path, required=True, help='New evidence directory')
    args = parser.parse_args()
    run = args.runs.resolve()
    if not (run/'shock_events').is_dir():
        run = run/'complete'
    assert (run/'shock_events').is_dir(), 'native complete run missing'
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    reference = args.reference.resolve() if args.reference else None
    if reference and (reference/'shock_events').is_dir():
        reference = reference/'shock_events'
    results, fixture_hashes, failure = [], {}, None

    def passed(name, carrier, **details):
        results.append(dict(id=name, carrier=carrier, status='pass', **details))

    def reject(name, carrier, action, diagnostic=None):
        try:
            action()
        except particles.HistoryError as error:
            if diagnostic:
                assert diagnostic in str(error), (name, str(error))
            passed(name, carrier, rejected=str(error))
        else:
            raise AssertionError('accepted corrupt '+name)

    try:
        for carrier, module, folder in ((0, particles, run/'shock_events'),
                                        (1, curves, run/'shock_events/curves')):
            files = binary.binary_files(folder, carrier)
            decoded = {path:binary.read(path) for path in files}
            assert all(p['nranks'] == 1 for p in decoded.values()), 'requires local Serial fixtures'
            committed = {packet['generation']:path for path, packet in decoded.items()
                         if packet['kind'] == 0}
            generations = sorted(committed)
            assert len(generations) >= 3, 'fixture needs at least three flushes'
            for path in files:
                fixture_hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
            actual = module.read_history(folder)
            assert actual, 'fixture needs completed events'
            passed('native_public_api', carrier, events=len(actual),
                   censored=len(module.read_suspended(folder)))
            if reference:
                previous = reference/'curves' if carrier else reference
                for name, reader in (('exact_csv_parity', module.read_history),
                                     ('exact_censored_csv_parity', module.read_suspended)):
                    assert normalized(reader(folder)) == normalized(reader(previous)), name
                    passed(name, carrier, records=len(reader(folder)))
            sampled = run.parent/'sampled'/('curves' if carrier else 'prtcl_thermo_history')
            snapshot = next(iter(sorted(sampled.glob('*.crv' if carrier else '*.thp'))))
            fixture_hashes[str(snapshot)] = hashlib.sha256(snapshot.read_bytes()).hexdigest()
            mixed = output/f'{carrier}-snapshot_same_directory'
            mixed.mkdir()
            for path in files:
                shutil.copy2(path, mixed/path.name)
            shutil.copy2(snapshot, mixed/('sampled-snapshot.'+snapshot.name))
            assert normalized(module.read_history(mixed)) == normalized(actual)
            passed('recognized_snapshot_ignored', carrier, snapshot=str(snapshot))
            unknown = output/f'{carrier}-unknown_magic'
            unknown.mkdir()
            for path in files:
                shutil.copy2(path, unknown/path.name)
            suffix = 'crv' if carrier else 'thp'
            data = bytearray(files[0].read_bytes())
            data[:7] = b'BROKEN!'
            (unknown/f'unknown.cycle0.rank0of1.{suffix}').write_bytes(data)
            reject('unknown_magic_not_ignored', carrier, lambda:module.read_history(unknown),
                   'Unknown/truncated binary encounter header')
            terminal = next(path for path, packet in decoded.items() if packet['kind'] == 1)
            ancestor_case = output/f'{carrier}-terminal_ancestor'
            ancestor_case.mkdir()
            for path in files:
                shutil.copy2(path, ancestor_case/path.name)
            data = bytearray(terminal.read_bytes())
            endian = layout(data)[0]
            previous_digest, = struct.unpack_from(endian+'Q', data, 56)
            struct.pack_into(endian+'Q', data, 56, previous_digest ^ 1)
            (ancestor_case/terminal.name).write_bytes(repair(data))
            reject('terminal_ancestor_corruption', carrier,
                   lambda:module.read_suspended(ancestor_case),
                   'Terminal view has missing/modified committed ancestor')
            for name in ('missing_middle', 'modified_prefix', 'modified_predecessor',
                         'duplicate_filename', 'missing_rank', 'query_coverage'):
                case = output/f'{carrier}-{name}'
                case.mkdir()
                if name == 'missing_rank':
                    initial = committed[generations[0]]
                    packet = decoded[initial]
                    assert packet['count'] == 0, 'rank fixture needs initial empty flush'
                    data = bytearray(initial.read_bytes())
                    endian, width, _ = layout(data)
                    struct.pack_into(endian+'2i', data, 80+width, 0, 2)
                    target = case/initial.name.replace('rank0of1', 'rank0of2')
                    target.write_bytes(repair(data))
                    reject(name, carrier, lambda:module.read_history(case),
                           'Missing binary encounter rank packet')
                    continue
                for path in files:
                    shutil.copy2(path, case/path.name)
                if name == 'missing_middle':
                    (case/committed[generations[1]].name).unlink()
                elif name in ('modified_prefix', 'modified_predecessor'):
                    path = case/committed[generations[1]].name
                    data = bytearray(path.read_bytes())
                    endian, width, _ = layout(data)
                    if name == 'modified_prefix':
                        value, = struct.unpack_from(endian+'Q', data, 64)
                        struct.pack_into(endian+'Q', data, 64, value ^ 1)
                    else:
                        real = 'd' if width == 8 else 'f'
                        value, = struct.unpack_from(endian+real, data, 80)
                        struct.pack_into(endian+real, data, 80, value+0.00001)
                    path.write_bytes(repair(data))
                elif name == 'duplicate_filename':
                    path = case/committed[generations[1]].name
                    shutil.copy2(path, case/path.name.replace('.cycle', '.cycle0'))
                else:
                    cutoffs = sorted({0, generations[-1],
                        (generations[0]+generations[-1])//2,
                        *(g for g in (50, 90, 100) if g <= generations[-1])})
                    for cutoff in cutoffs:
                        expected = [e for e in actual if e['emission_generation'] <= cutoff]
                        assert normalized(module.read_history(case, cutoff)) == normalized(expected)
                    # Isolate coverage from terminal-ancestor validation: canonical
                    # final views can link the latest flush that this case removes.
                    for path, packet in decoded.items():
                        if packet['kind'] == 1:
                            (case/path.name).unlink()
                    (case/committed[generations[-1]].name).unlink()
                    reject(name, carrier, lambda:module.read_history(case, generations[-1]),
                           'not covered by committed batches')
                    results[-1]['covered_queries'] = cutoffs
                    continue
                reject(name, carrier, lambda:module.read_history(case))
            event_path = next(path for path, packet in decoded.items()
                              if packet['kind'] == 0 and packet['count'])
            native = event_path.read_bytes()
            endian, width, record = layout(native)
            observation_size = 40+18*width
            mutations = {
                'unknown_flags': (record+24, 'I', 1024),
                'invalid_status': (record+20, 'I', 2),
                'invalid_validity': (record+28+width+8, 'I', 2),
                'nonfinite_inside_time': (record+28+observation_size,
                                          'd' if width == 8 else 'f', float('nan')),
                'unknown_real_width': (24, 'I', 16),
            }
            for name, (offset, fmt, value) in mutations.items():
                data = bytearray(native)
                struct.pack_into(endian+fmt, data, offset, value)
                # Invalid width prevents layout() from walking strings; checksum
                # uses the already independently decoded native byte order.
                data[-8:] = struct.pack(endian+'Q', binary.digest(data[:-8]))
                path = output/f'{carrier}-{name}.packet'
                path.write_bytes(data)
                reject(name, carrier, lambda:binary.read(path))
            for name, data in (('truncated_header', native[:20]),
                               ('truncated_record', bytes(repair(bytearray(native[:-16])))),
                               ('checksum_corruption', native[:-9]+bytes([native[-9] ^ 1])+native[-8:])):
                path = output/f'{carrier}-{name}.packet'
                path.write_bytes(data)
                reject(name, carrier, lambda:binary.read(path))
            if carrier == 0:
                event = copy.deepcopy(actual[0])
                header = decoded[committed[generations[-1]]]
                event['post'] = dict(time=header['time'], generation=header['generation'],
                    valid=False, inside=False, x=[1.1, 0.5, 0.0], image=[0, 0, 0],
                    rho=0, pressure=0, temperature=0, mach_pressure=0, mach_temperature=0,
                    normal=[0, 0, 0], upstream=[0, 0, 0], downstream=[0, 0, 0])
                for kind in (0, 2):
                    path = output/f'domain_exit_kind{kind}.packet'
                    path.write_bytes(domain_exit_packet(event, header, kind, 192+8+16))
                    restored = binary.read(path)['events'][0]
                    assert restored['post'] == event['post'] and restored['flags'] == 216
                    passed(f'domain_exit_kind{kind}', carrier, fixture='synthetic',
                           emission=restored['emission_generation'])
                    future = bytearray(path.read_bytes())
                    endian, width, record = layout(future)
                    post_generation = record+28+3*(40+18*width)+width
                    struct.pack_into(endian+'Q', future, post_generation, header['generation']+1)
                    path = output/f'future_domain_exit_generation_kind{kind}.packet'
                    path.write_bytes(repair(future))
                    reject(f'future_domain_exit_generation_kind{kind}', carrier,
                           lambda:binary.read(path), 'Observation after emission generation')
                    path = output/f'terminal_only_kind{kind}.packet'
                    path.write_bytes(domain_exit_packet(event, header, kind, 64+8+16))
                    reject(f'terminal_only_kind{kind}', carrier, lambda:binary.read(path))
                initial = copy.deepcopy(next(e for e in actual if e['status'] == 0))
                initial['pre'] = dict(event['post'], time=-1, generation=0)
                path = output/'left_censored_sentinel.packet'
                path.write_bytes(domain_exit_packet(initial, header, 0, 1+8+16))
                restored = binary.read(path)['events'][0]
                assert restored['pre']['time'] == -1 and restored['pre']['generation'] == 0
                passed('left_censored_sentinel_preserved', carrier, fixture='synthetic')
    except Exception as error:
        failure = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        sources = [Path(__file__), *(ROOT/'vis/python'/name for name in
                   ('shock_binary.py', 'shock_history.py', 'curve_shock_history.py'))]
        report = dict(status='fail' if failure else 'pass', failure=failure,
            command=shlex.join(['python3', *sys.argv]), run=str(run),
            reference=str(reference) if reference else None, checks=results,
            source_sha256={str(path):hashlib.sha256(path.read_bytes()).hexdigest()
                           for path in sources}, fixture_sha256=fixture_hashes,
            limitation='Reader/serialization checks; no live outflow, MPI or GPU qualification')
        (output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(f'{len(results)} binary encounter reader checks PASS; {output}/report.json')


if __name__ == '__main__':
    main()
