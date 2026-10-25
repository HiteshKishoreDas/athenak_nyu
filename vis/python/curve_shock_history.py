"""Strict curve-node shock histories; censored views remain separate from events.

Node IDs use a declared curve/sequence split; fixed 8/24 and historical 32/32 remain readable.
AMR level/spacing describe the sampled cells, not a fixed physical detector scale.
"""
import csv
from material_curves import node_curve_id
import math
from pathlib import Path
import re

from shock_history import HistoryError, integer, parse_event, parse_observation


def _real(text):
    try:
        value = float(text)
    except ValueError as exc:
        raise HistoryError('Invalid curve metadata real') from exc
    if not math.isfinite(value):
        raise HistoryError('Nonfinite curve metadata')
    return value


def id_layout(contract):
    """Contract prefixes are authoritative even when journals have no events."""
    if contract.startswith('curve-node-history-v3 '):
        fields=contract.split(maxsplit=2)
        if len(fields)!=3 or not re.fullmatch(r'[0-9]+',fields[1]) or int(fields[1])>30:
            raise HistoryError('Invalid curve_id_bits contract')
        return (32,int(fields[1]))
    return (32,8) if contract.startswith('curve-node-history-v2 ') else (64,32)


def _node_id(text, layout=None):
    node_id = integer(text)
    try:
        curve = (node_curve_id(node_id,curve_id_bits=layout[1]) if layout and layout[0]==32
                 else node_curve_id(node_id,32) if layout else node_curve_id(node_id))
    except ValueError as exc:
        raise HistoryError('Invalid curve node ID namespace') from exc
    return node_id,curve


def _spacing(fields):
    values = [_real(value) for value in fields]
    if any(value <= 0 for value in values):
        raise HistoryError('Curve sample spacing must be positive')
    return values


def manifests(directory):
    result = {}
    for path in sorted(Path(directory).glob('manifest.rank*')):
        match = re.fullmatch(r'manifest\.rank(\d+)', path.name)
        if not match:
            raise HistoryError('Malformed curve manifest filename')
        rank = integer(match[1], 32)
        if rank in result:
            raise HistoryError('Duplicate curve manifest rank alias')
        try:
            lines = path.read_text().splitlines()
        except (OSError, UnicodeError) as exc:
            raise HistoryError('Cannot read curve manifest') from exc
        if len(lines) != 3 or lines[0] != 'AKCURVE-SHOCK-EVENTS 1':
            raise HistoryError('Unknown curve history schema')
        if not re.fullmatch(r'[A-Za-z0-9_-]{1,128}', lines[1]) or not lines[2]:
            raise HistoryError('Invalid curve lineage/contract')
        id_layout(lines[2])
        result[rank] = (lines[1], lines[2])
    if not result or sorted(result) != list(range(len(result))):
        raise HistoryError('Missing/noncontiguous curve rank manifests')
    if len(set(result.values())) != 1:
        raise HistoryError('Curve rank contracts disagree')
    return result


def _files(directory, prefix):
    result = {}
    for path in Path(directory).glob(prefix+'*.rank*.csv'):
        match = re.fullmatch(prefix+r'(\d+)\.rank(\d+)\.csv', path.name)
        if not match:
            raise HistoryError('Malformed curve chunk filename')
        generation, rank = (integer(value) for value in match.groups())
        if (rank, generation) in result:
            raise HistoryError('Duplicate curve generation alias')
        result[rank, generation] = path
    return result


def _read_chunk(path, lineage, contract, generation, suspended=False):
    try:
        lines = Path(path).read_text().splitlines()
    except (OSError, UnicodeError) as exc:
        raise HistoryError('Cannot read curve chunk') from exc
    kind = 'AKCURVE-SHOCK-SUSPENDED' if suspended else 'AKCURVE-SHOCK-EVENTS'
    header = lines[0].split() if lines else []
    if len(lines) < 3 or len(header) != 6 or header[:3] != ['#', kind, '1']:
        raise HistoryError('Unknown/truncated curve chunk schema')
    count = integer(header[5])
    if header[3] != lineage or integer(header[4]) != generation:
        raise HistoryError('Curve lineage/generation mismatch')
    if lines[1] != '# '+contract or lines[-1] != '# END '+str(count):
        raise HistoryError('Curve contract/END mismatch')
    if len(lines) != count+3:
        raise HistoryError('Curve row count mismatch')
    layout=id_layout(contract)
    events = []
    for row in csv.reader(lines[2:-1]):
        if len(row) != 133:
            raise HistoryError('Curve event requires133fields')
        event = parse_event(row[:124], lineage, generation)
        node_id,curve = _node_id(row[0],layout)
        birth = _real(row[124])
        first_level = integer(row[125], 32, signed=True)
        last_level = integer(row[129], 32, signed=True)
        if birth < 0 or birth > event['first']['time'] or min(first_level, last_level) < 0:
            raise HistoryError('Curve birth/refinement metadata mismatch')
        event.update(carrier='curve_node',
                     identity=('curve_node',lineage,node_id,event['sequence']),
                     node_id=node_id, curve_id=curve, birth=birth,
                     first_level=first_level, first_spacing=_spacing(row[126:129]),
                     last_level=last_level, last_spacing=_spacing(row[130:133]))
        if suspended and (event['status'] != 1 or not event['flags'] & 64):
            raise HistoryError('Curve suspended view is not right censored')
        events.append(event)
    return events


def read_history(directory, through_cycle=None):
    """Return canonical events; require every per-rank generation including empty ones."""
    metadata = manifests(directory)
    paths = _files(directory, 'cycle')
    if any(rank not in metadata for rank, _ in paths):
        raise HistoryError('Curve chunk lacks rank manifest')
    if through_cycle is None:
        if not paths:
            raise HistoryError('No committed curve generations')
        through_cycle = max(generation for _, generation in paths)
    if not isinstance(through_cycle, int) or through_cycle < 0:
        raise HistoryError('through_cycle must be a nonnegative integer')
    result, identities = [], set()
    for rank, (lineage, contract) in metadata.items():
        for generation in range(through_cycle+1):
            path = paths.get((rank, generation))
            if path is None:
                raise HistoryError(f'Missing curve generation{generation} rank{rank}')
            for event in _read_chunk(path, lineage, contract, generation):
                if event['identity'] in identities:
                    raise HistoryError('Duplicate curve event identity')
                identities.add(event['identity'])
                event['rank'] = rank
                result.append(event)
    return sorted(result, key=lambda event: event['identity'])


def read_suspended(directory, generation=None):
    """Read a terminal censored view without adding it to committed events."""
    metadata = manifests(directory)
    paths = _files(directory, 'suspended')
    if any(rank not in metadata for rank, _ in paths):
        raise HistoryError('Curve suspended view lacks rank manifest')
    if generation is None:
        if not paths:
            return []
        generation = max(g for _, g in paths)
    if not isinstance(generation, int) or generation < 0:
        raise HistoryError('generation must be a nonnegative integer')
    result, identities = [], set()
    for rank, (lineage, contract) in metadata.items():
        path = paths.get((rank, generation))
        if path is None:
            raise HistoryError('Missing rank curve suspended view')
        for event in _read_chunk(path, lineage, contract, generation, suspended=True):
            if event['identity'] in identities:
                raise HistoryError('Duplicate curve suspended identity')
            identities.add(event['identity'])
            event['rank'] = rank
            result.append(event)
    return sorted(result, key=lambda event: event['identity'])


def read_trace(path, curve_id_bits=None):
    """Read one raw33-column trace/birthtrace file, keyed by exact node ID."""
    path = Path(path)
    match = re.fullmatch(r'(?:birthtrace|trace)(\d+)\.rank(\d+)\.csv', path.name)
    if not match:
        raise HistoryError('Malformed curve trace filename')
    generation, rank = (integer(value) for value in match.groups())
    result = {}
    if curve_id_bits is not None and (type(curve_id_bits) is not int or not 0<=curve_id_bits<=30):
        raise HistoryError('Invalid curve_id_bits')
    layout = None if curve_id_bits is None else (32,curve_id_bits)
    if (path.parent/f'manifest.rank{rank}').exists():
        contract = manifests(path.parent)[rank][1]
        declared = id_layout(contract)
        if layout is not None and layout!=declared:
            raise HistoryError('Curve trace layout differs from manifest')
        layout = declared
    try:
        lines = path.read_text().splitlines()
    except (OSError, UnicodeError) as exc:
        raise HistoryError('Cannot read curve trace') from exc
    for row in csv.reader(lines):
        if len(row) != 33:
            raise HistoryError('Curve trace requires33fields')
        if layout is None and integer(row[0]) < 1<<32:
            raise HistoryError('32-bit curve trace requires manifest or curve_id_bits')
        node_id,curve = _node_id(row[0],layout or (64,32))
        sequence, phase, flags = integer(row[1]), integer(row[2], 32), integer(row[3], 32)
        if node_id in result or phase > 2 or flags & ~1023:
            raise HistoryError('Duplicate/invalid curve trace identity/state')
        observation = parse_observation(row[4:28])
        birth, level = _real(row[28]), integer(row[29], 32, signed=True)
        if observation['generation'] != generation or birth < 0 or birth > observation['time'] or level < 0:
            raise HistoryError('Curve trace generation/birth/refinement mismatch')
        result[node_id] = dict(observation, node_id=node_id, curve_id=curve,
                               particle_id=node_id, sequence=sequence, phase=phase,
                               flags=flags, birth=birth, level=level,
                               spacing=_spacing(row[30:33]), rank=rank)
    return result
