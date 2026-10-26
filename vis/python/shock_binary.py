"""Strict typed encounter batches; sampled THP/CRV snapshots are separate files.

Footer checksums hash the body. Previous-file links hash the whole prior file,
including its footer; prefix links roll that whole file into the prior prefix.
Emission clocks are preserved independently of output flush clocks. Historical
CSV readers supply the shared event semantic validator used here.
"""
import math
from pathlib import Path
import re
import struct

from shock_history import HistoryError, parse_event

OFFSET = 14695981039346656037
PRIME = 1099511628211
MASK = (1 << 64)-1
MAGIC = b'ATHK_EVENTS'+b'\0'*5


def _require(condition, message):
    if not condition:
        raise HistoryError(message)


def digest(data, seed=OFFSET):
    for byte in data:
        seed = ((seed ^ byte)*PRIME) & MASK
    return seed


def binary_files(directory, carrier):
    """Event candidates; recognized snapshot formats can share this directory."""
    _require(carrier in (0, 1), 'Invalid binary encounter carrier')
    result = []
    for path in sorted(Path(directory).glob('*.crv' if carrier else '*.thp')):
        try:
            with path.open('rb') as stream:
                magic = stream.read(32)
        except OSError as exc:
            raise HistoryError('Cannot inspect binary encounter candidate') from exc
        if ((carrier and magic[:16] == b'ATHKCRV'+b'\0'*9) or
                (not carrier and magic == b'ATHK_PRTCL_THERMO_HISTORY'.ljust(32, b'\0'))):
            continue
        # Unknown/truncated magic stays a candidate: corruption must fail loudly.
        result.append(path)
    return result


def read(path, max_events=1000000):
    """Read one completed, terminal or pending packet, with bounded allocation."""
    _require(type(max_events) is int and max_events >= 0, 'Invalid event reader cap')
    path = Path(path)
    try:
        _require(path.stat().st_size <= 21000+1012*max_events,
                 'Binary encounter file exceeds reader cap')
        data = path.read_bytes()
    except OSError as exc:
        raise HistoryError('Cannot read binary encounter packet') from exc
    _require(len(data) >= 112 and data[:16] == MAGIC,
             'Unknown/truncated binary encounter header')
    endian = ('<' if data[20:24] == b'\x04\x03\x02\x01' else
              '>' if data[20:24] == b'\x01\x02\x03\x04' else None)
    _require(endian is not None, 'Invalid binary encounter endian marker')
    footer = struct.unpack_from(endian+'Q', data, len(data)-8)[0]
    _require(digest(data[:-8]) == footer, 'Binary encounter checksum mismatch')
    offset = 16

    def take(fmt):
        nonlocal offset
        parser = struct.Struct(endian+fmt)
        _require(offset+parser.size <= len(data)-8, 'Truncated binary encounter records')
        values = parser.unpack_from(data, offset)
        offset += parser.size
        return values

    version, marker, width, carrier, bits, kind = take('6I')
    _require(version == 1 and width in (4, 8) and carrier in (0, 1) and
             kind in (0, 1, 2) and (bits <= 30 if carrier else bits == 0),
             'Unsupported binary encounter schema/carrier/split')
    real = 'd' if width == 8 else 'f'
    generation, previous, previous_digest, prefix, count = take('5Q')
    time, = take(real)
    rank, nranks = take('2i')
    _require(math.isfinite(time) and time >= 0 and 0 <= rank < nranks and
             count <= max_events and (previous == MASK or previous <= generation),
             'Invalid binary encounter clock/rank/count')

    def text(cap):
        nonlocal offset
        size, = take('I')
        _require(0 < size <= cap and offset+size <= len(data)-8,
                 'Invalid binary encounter string extent')
        try:
            value = data[offset:offset+size].decode('utf-8')
        except UnicodeDecodeError as exc:
            raise HistoryError('Invalid binary encounter UTF8') from exc
        offset += size
        _require('\0' not in value and '\n' not in value and '\r' not in value,
                 'Invalid binary encounter text')
        return value

    lineage, contract, stem = text(128), text(16384), text(4096)
    _require(re.fullmatch(r'[A-Za-z0-9_-]{1,128}', lineage) is not None and
             '/' not in stem and '\\' not in stem, 'Invalid binary encounter lineage/stem')
    if carrier:
        from curve_shock_history import id_layout, _node_id, _spacing
        _require(contract.startswith('curve-node-history-v4 ') and
                 id_layout(contract) == (32, bits, True),
                 'Binary curve encounter contract/split mismatch')
    event_size = 28+5*(40+18*width)+7*width+8
    _require(offset+count*event_size == len(data)-8,
             'Binary encounter record count/trailing bytes mismatch')
    events, identities = [], set()
    for _ in range(count):
        emission, identity, sequence, status, flags = take('QiQ2I')
        _require(emission <= generation and
                 (kind == 1 or previous == MASK or emission > previous),
                 'Event emission outside flush interval')
        row = [str(v) for v in (identity, sequence, status, flags)]
        for _ in range(5):
            clock, = take(real)
            observed, valid, inside = take('Q2I')
            x = take(3*real)
            image = take('3q')
            thermo = take(14*real)
            row.extend(str(v) for v in (clock, observed, valid, inside, *x, *image, *thermo))
        birth, = take(real)
        first_level, last_level = take('2i')
        first_spacing, last_spacing = take(3*real), take(3*real)
        event = parse_event(row, lineage, emission, signed_id=bool(carrier))
        _require(status == int(bool(flags & (1 | 64 | 128))),
                 'Event partial status/censoring flags disagree')
        _require(all(not event[k]['valid'] or event[k]['time'] <= time
                     for k in ('pre', 'first', 'last', 'post', 'representative')),
                 'Observation after flush time')
        _require(event['last']['valid'] and event['last']['inside'] and
                 event['representative']['valid'] and event['representative']['inside'] and
                 (not event['pre']['valid'] or not event['pre']['inside']) and
                 (not event['post']['valid'] or not event['post']['inside']),
                 'Invalid event inside/outside evidence')
        first, last = event['first'], event['last']
        representative = event['representative']
        _require(first['generation'] <= last['generation'] and
                 (not event['pre']['valid'] or
                  event['pre']['generation'] <= first['generation']) and
                 (not event['post']['valid'] or
                  last['generation'] <= event['post']['generation']) and
                 first['time'] <= representative['time'] <= last['time'] and
                 first['generation'] <= representative['generation'] <= last['generation'],
                 'Reversed event generation/representative bracket')
        if carrier:
            node_id, curve = _node_id(str(identity), (32, bits, True))
            _require(math.isfinite(birth) and 0 <= birth <= event['first']['time'] and
                     min(first_level, last_level) >= 0 and curve % nranks == rank,
                     'Invalid binary curve birth/refinement/owner metadata')
            event.update(carrier='curve_node', node_id=node_id, curve_id=curve,
                         identity=('curve_node', lineage, node_id, sequence), birth=birth,
                         first_level=first_level, last_level=last_level,
                         first_spacing=_spacing([str(v) for v in first_spacing]),
                         last_spacing=_spacing([str(v) for v in last_spacing]))
        else:
            _require(identity >= 0 and birth == -1 and first_level == last_level == -1 and
                     first_spacing == last_spacing == (0, 0, 0),
                     'Invalid binary particle identity/metadata')
        if kind == 1:
            _require(emission == generation and status == 1 and flags & 64,
                     'Terminal view is not current right-censored evidence')
        else:
            _require(not flags & 64 or flags & 128,
                     'Committed event carries terminal-only censoring flag')
        _require(event['identity'] not in identities, 'Duplicate event identity in packet')
        identities.add(event['identity'])
        event['rank'] = rank
        events.append(event)
    return dict(path=str(path), version=version, precision=8*width, carrier=carrier,
                curve_id_bits=bits, kind=kind, generation=generation, time=time,
                previous_flush_generation=previous, previous_file_digest=previous_digest,
                prefix_digest=prefix, count=count, rank=rank, nranks=nranks,
                lineage=lineage, contract=contract, stem=stem, events=events,
                checksum=footer, file_digest=digest(data), rolling_digest=digest(data, prefix))


def _packets(directory, carrier):
    paths = binary_files(directory, carrier)
    _require(paths, 'No binary encounter packets')
    suffix = 'crv' if carrier else 'thp'
    packets, groups, metadata = [], {}, None
    _require(not list(Path(directory).glob('cycle*.rank*.csv')) and
             not list(Path(directory).glob('suspended*.rank*.csv')),
             'Mixed historical CSV and binary encounter streams')
    for path in paths:
        match = re.fullmatch(r'(.+)\.(cycle|suspended)([0-9]+)\.rank([0-9]+)of([0-9]+)\.'+suffix,
                            path.name)
        _require(match is not None, 'Malformed binary encounter filename')
        packet = read(path)
        expected_kind = 0 if match[2] == 'cycle' else 1
        _require(packet['carrier'] == carrier and packet['kind'] == expected_kind and
                 (packet['stem'], packet['generation'], packet['rank'], packet['nranks']) ==
                 (match[1], int(match[3]), int(match[4]), int(match[5])),
                 'Binary encounter filename/header mismatch')
        current = tuple(packet[k] for k in
                        ('lineage', 'contract', 'stem', 'precision', 'curve_id_bits', 'nranks'))
        _require(metadata is None or metadata == current, 'Binary encounter stream metadata changed')
        metadata = current
        group = groups.setdefault((packet['kind'], packet['generation']), {})
        _require(packet['rank'] not in group, 'Duplicate binary encounter rank alias')
        group[packet['rank']] = packet
        packets.append(packet)
    for group in groups.values():
        first = next(iter(group.values()))
        _require(len(group) == first['nranks'], 'Missing binary encounter rank packet')
        _require(all(p['time'] == first['time'] and
                     p['previous_flush_generation'] == first['previous_flush_generation']
                     for p in group.values()), 'Binary encounter rank clocks/coverage disagree')
    for rank in range(metadata[-1]):
        chain = {MASK: (0, OFFSET, -1)}
        previous, prior_time = MASK, -1
        identities = set()
        for packet in sorted((p for p in packets if p['rank'] == rank and p['kind'] == 0),
                             key=lambda p:p['generation']):
            _require(packet['previous_flush_generation'] == previous and
                     packet['previous_file_digest'] == chain[previous][0] and
                     packet['prefix_digest'] == chain[previous][1] and
                     packet['time'] > prior_time and
                     (previous == MASK or packet['generation'] > previous),
                     'Missing/modified binary encounter chain predecessor')
            for event in packet['events']:
                _require(event['identity'] not in identities, 'Duplicate committed event identity')
                identities.add(event['identity'])
            previous, prior_time = packet['generation'], packet['time']
            chain[previous] = (packet['file_digest'], packet['rolling_digest'], packet['time'])
        for packet in (p for p in packets if p['rank'] == rank and p['kind'] == 1):
            ancestor = chain.get(packet['previous_flush_generation'])
            _require(ancestor is not None and packet['previous_file_digest'] == ancestor[0] and
                     packet['prefix_digest'] == ancestor[1] and packet['time'] >= ancestor[2],
                     'Terminal view has missing/modified committed ancestor')
    # Optional qualification traces retain text manifests. They must describe this stream.
    manifest_ranks = set()
    for path in Path(directory).glob('manifest.rank*'):
        match = re.fullmatch(r'manifest\.rank([0-9]+)', path.name)
        _require(match is not None and int(match[1]) < metadata[-1],
                 'Invalid binary encounter trace manifest rank')
        rank = int(match[1])
        _require(rank not in manifest_ranks, 'Duplicate binary trace manifest rank alias')
        manifest_ranks.add(rank)
        try:
            lines = path.read_text().splitlines()
        except (OSError, UnicodeError) as exc:
            raise HistoryError('Cannot read encounter trace manifest') from exc
        _require(lines == [('AKCURVE-SHOCK-EVENTS 1' if carrier else 'AKSHOCK-EVENTS 1'),
                           metadata[0], metadata[1]], 'Binary encounter trace manifest differs')
    _require(not manifest_ranks or len(manifest_ranks) == metadata[-1],
             'Missing binary encounter trace manifest rank')
    return packets


def manifests(directory, carrier):
    packets = _packets(directory, carrier)
    first = packets[0]
    return {rank:(first['lineage'], first['contract']) for rank in range(first['nranks'])}


def read_history(directory, carrier, through_cycle=None):
    packets = [p for p in _packets(directory, carrier) if p['kind'] == 0]
    _require(packets, 'No committed binary encounter batches')
    latest = max(p['generation'] for p in packets)
    if through_cycle is None:
        through_cycle = latest
    _require(type(through_cycle) is int and 0 <= through_cycle <= latest,
             'Requested generation is not covered by committed batches')
    events = [e for p in packets for e in p['events']
              if e['emission_generation'] <= through_cycle]
    identities = [e['identity'] for e in events]
    _require(len(set(identities)) == len(identities), 'Duplicate committed event identity across ranks')
    return sorted(events, key=lambda e:e['identity'])


def read_suspended(directory, carrier, generation=None):
    packets = [p for p in _packets(directory, carrier) if p['kind'] == 1]
    if generation is None:
        if not packets:
            return []
        generation = max(p['generation'] for p in packets)
    _require(type(generation) is int and generation >= 0, 'Invalid terminal generation')
    selected = [p for p in packets if p['generation'] == generation]
    _require(selected, 'Missing binary terminal view')
    events = [e for p in selected for e in p['events']]
    identities = [e['identity'] for e in events]
    _require(len(set(identities)) == len(identities), 'Duplicate binary terminal identity')
    return sorted(events, key=lambda e:e['identity'])
