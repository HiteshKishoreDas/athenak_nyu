"""Strict curve-node shock histories; censored views remain separate from events.

Node IDs retain all64bits. Curve ID is the high32bits of the native node ID.
AMR level/spacing describe the sampled cells, not a fixed physical detector scale.
"""
import csv
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


def _node_id(text):
    node_id = integer(text)
    if not node_id >> 32 or not node_id & ((1 << 32) - 1):
        raise HistoryError('Invalid curve node ID namespace')
    return node_id


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
    events = []
    for row in csv.reader(lines[2:-1]):
        if len(row) != 133:
            raise HistoryError('Curve event requires133fields')
        event = parse_event(row[:124], lineage, generation)
        node_id = _node_id(row[0])
        birth = _real(row[124])
        first_level = integer(row[125], 32, signed=True)
        last_level = integer(row[129], 32, signed=True)
        if birth < 0 or birth > event['first']['time'] or min(first_level, last_level) < 0:
            raise HistoryError('Curve birth/refinement metadata mismatch')
        event.update(carrier='curve_node',
                     identity=('curve_node',lineage,node_id,event['sequence']),
                     node_id=node_id, curve_id=node_id >> 32, birth=birth,
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


def read_trace(path):
    """Read one raw33-column trace/birthtrace file, keyed by exact node ID."""
    path = Path(path)
    match = re.fullmatch(r'(?:birthtrace|trace)(\d+)\.rank(\d+)\.csv', path.name)
    if not match:
        raise HistoryError('Malformed curve trace filename')
    generation, rank = (integer(value) for value in match.groups())
    result = {}
    try:
        lines = path.read_text().splitlines()
    except (OSError, UnicodeError) as exc:
        raise HistoryError('Cannot read curve trace') from exc
    for row in csv.reader(lines):
        if len(row) != 33:
            raise HistoryError('Curve trace requires33fields')
        node_id = _node_id(row[0])
        sequence, phase, flags = integer(row[1]), integer(row[2], 32), integer(row[3], 32)
        if node_id in result or phase > 2 or flags & ~1023:
            raise HistoryError('Duplicate/invalid curve trace identity/state')
        observation = parse_observation(row[4:28])
        birth, level = _real(row[28]), integer(row[29], 32, signed=True)
        if observation['generation'] != generation or birth < 0 or birth > observation['time'] or level < 0:
            raise HistoryError('Curve trace generation/birth/refinement mismatch')
        result[node_id] = dict(observation, node_id=node_id, curve_id=node_id >> 32,
                               particle_id=node_id, sequence=sequence, phase=phase,
                               flags=flags, birth=birth, level=level,
                               spacing=_spacing(row[30:33]), rank=rank)
    return result
