"""Strict reader for experimental AthenaK shock event schema 1.

Event identities are (lineage, particle_id, sequence), retaining integer widths.
Suspended terminal views are separate from committed events. Times describe mask
entry/exit evidence, never an inferred physical crossing time.
"""
import csv
import math
from pathlib import Path
import re


class HistoryError(ValueError):
    """Malformed, incomplete or inconsistent shock history."""


def integer(text, bits=64, signed=False):
    if not re.fullmatch(r'-?\d+' if signed else r'\d+', text):
        raise HistoryError('Invalid integer: ' + text)
    value = int(text)
    low = -(1 << (bits - 1)) if signed else 0
    high = (1 << (bits - int(signed))) - 1
    if not low <= value <= high:
        raise HistoryError('Integer outside schema range: ' + text)
    return value


def parse_observation(fields):
    if len(fields) != 24:
        raise HistoryError('Observation requires 24 fields')
    valid, inside = integer(fields[2], 1), integer(fields[3], 1)
    floats = {}
    names = {0: 'time', 10: 'rho', 11: 'pressure', 12: 'temperature',
             13: 'mach_pressure', 14: 'mach_temperature'}
    for n in (0, 4, 5, 6, *range(10, 24)):
        try:
            value = float(fields[n])
        except ValueError as exc:
            raise HistoryError('Invalid floating field') from exc
        if math.isinf(value) or (math.isnan(value) and valid and n not in (13, 14)):
            raise HistoryError('Nonfinite required observation field')
        floats[n] = value
    if not math.isfinite(floats[0]):
        raise HistoryError('Observation time must be finite')
    if valid and (floats[0] < 0 or any(floats[n] <= 0 for n in (10, 11, 12, *range(18, 24)))):
        raise HistoryError('Valid observation requires positive thermodynamics')
    result = {name: floats[n] for n, name in names.items()}
    result.update(generation=integer(fields[1]), valid=bool(valid), inside=bool(inside),
                  x=[floats[n] for n in range(4, 7)],
                  image=[integer(fields[n], signed=True) for n in range(7, 10)],
                  normal=[floats[n] for n in range(15, 18)],
                  upstream=[floats[n] for n in range(18, 21)],
                  downstream=[floats[n] for n in range(21, 24)])
    if inside and not valid:
        raise HistoryError('Inside observation must be valid')
    return result


def parse_event(row, lineage, generation):
    """Parse the unchanged 124-column event core shared by both carriers."""
    if len(row) != 124:
        raise HistoryError('Event requires 124 fields')
    event = dict(lineage=lineage, particle_id=integer(row[0]),
                 sequence=integer(row[1]), status=integer(row[2], 1),
                 flags=integer(row[3], 32), emission_generation=generation)
    event['identity'] = (lineage, event['particle_id'], event['sequence'])
    if event['flags'] & ~1023:
        raise HistoryError('Unknown schema-1 validity flag')
    if event['status'] == 0 and event['flags'] & (1 | 64 | 128):
        raise HistoryError('Completed event carries censoring flags')
    for n, name in enumerate(('pre', 'first', 'last', 'post', 'representative')):
        o = parse_observation(row[4 + 24*n:4 + 24*(n+1)])
        if o['valid'] and o['generation'] > generation:
            raise HistoryError('Observation after emission generation')
        event[name] = o
    if not event['first']['valid'] or not event['first']['inside']:
        raise HistoryError('Event has no valid first inside evidence')
    if event['last']['time'] < event['first']['time']:
        raise HistoryError('Reversed inside interval')
    if event['pre']['valid'] and event['pre']['time'] > event['first']['time']:
        raise HistoryError('Reversed entry bracket')
    if event['post']['valid'] and event['post']['time'] < event['last']['time']:
        raise HistoryError('Reversed exit bracket')
    if event['status'] == 0 and not event['post']['valid']:
        raise HistoryError('Completed event lacks outside evidence')
    return event


def read_chunk(path, lineage, contract, generation, suspended=False):
    try:
        lines = Path(path).read_text().splitlines()
    except (OSError, UnicodeError) as exc:
        raise HistoryError('Cannot read chunk: ' + str(path)) from exc
    kind = 'AKSHOCK-SUSPENDED' if suspended else 'AKSHOCK-EVENTS'
    if len(lines) < 3:
        raise HistoryError('Truncated chunk: ' + str(path))
    header = lines[0].split()
    if len(header) != 6 or header[:3] != ['#', kind, '1']:
        raise HistoryError('Unknown/malformed chunk schema: ' + str(path))
    count = integer(header[5])
    if header[3] != lineage or integer(header[4]) != generation:
        raise HistoryError('Chunk lineage/generation mismatch')
    if lines[1] != '# ' + contract or lines[-1] != '# END ' + str(count):
        raise HistoryError('Chunk contract or END mismatch')
    if len(lines) != count + 3:
        raise HistoryError('Chunk row count mismatch')
    events = []
    for row in csv.reader(lines[2:-1]):
        events.append(parse_event(row, lineage, generation))
    return events


def manifests(directory):
    root = Path(directory)
    paths = sorted(root.glob('manifest.rank*'))
    result = {}
    for path in paths:
        match = re.fullmatch(r'manifest\.rank(\d+)', path.name)
        if not match:
            raise HistoryError('Malformed manifest name')
        rank = int(match[1])
        if rank in result:
            raise HistoryError('Duplicate rank manifest alias')
        lines = path.read_text().splitlines()
        if len(lines) != 3 or lines[0] != 'AKSHOCK-EVENTS 1':
            raise HistoryError('Unknown/malformed manifest schema')
        if not re.fullmatch(r'[A-Za-z0-9_-]{1,128}', lines[1]) or not lines[2]:
            raise HistoryError('Invalid manifest lineage/contract')
        result[rank] = (lines[1], lines[2])
    if not result or sorted(result) != list(range(len(result))):
        raise HistoryError('Missing/ noncontiguous rank manifests')
    if len(set(result.values())) != 1:
        raise HistoryError('Rank manifest contracts disagree')
    return result


def read_history(directory, through_cycle=None):
    """Return canonical committed events; require contiguous per-rank chunks."""
    root = Path(directory)
    metadata = manifests(root)
    chunks = {}
    for path in root.glob('cycle*.rank*.csv'):
        match = re.fullmatch(r'cycle(\d+)\.rank(\d+)\.csv', path.name)
        if not match:
            raise HistoryError('Malformed chunk filename')
        generation, rank = map(int, match.groups())
        if rank not in metadata:
            raise HistoryError('Chunk has no rank manifest')
        if (rank, generation) in chunks:
            raise HistoryError('Duplicate chunk generation alias')
        chunks[rank, generation] = path
    if through_cycle is None:
        if not chunks:
            raise HistoryError('No committed generations')
        through_cycle = max(g for _, g in chunks)
    if not isinstance(through_cycle, int) or through_cycle < 0:
        raise HistoryError('through_cycle must be a nonnegative integer')
    events, identities = [], set()
    for rank, (lineage, contract) in metadata.items():
        for generation in range(through_cycle + 1):
            path = chunks.get((rank, generation))
            if path is None:
                raise HistoryError(f'Missing generation {generation} on rank {rank}')
            for event in read_chunk(path, lineage, contract, generation):
                if event['identity'] in identities:
                    raise HistoryError('Duplicate event identity')
                identities.add(event['identity'])
                event['rank'] = rank
                events.append(event)
    return sorted(events, key=lambda event: event['identity'])


def read_suspended(directory, generation=None):
    """Read a selected/latest terminal view, separately from committed history."""
    root = Path(directory)
    metadata = manifests(root)
    paths = {}
    for path in root.glob('suspended*.rank*.csv'):
        match = re.fullmatch(r'suspended(\d+)\.rank(\d+)\.csv', path.name)
        if not match:
            raise HistoryError('Malformed suspended filename')
        g, rank = map(int, match.groups())
        if rank not in metadata or (rank, g) in paths:
            raise HistoryError('Unknown rank or duplicate suspended alias')
        paths[rank, g] = path
    if generation is None:
        if not paths:
            return []
        generation = max(g for _, g in paths)
    result, identities = [], set()
    for rank, (lineage, contract) in metadata.items():
        path = paths.get((rank, generation))
        if path is None:
            raise HistoryError('Missing rank suspended view')
        for event in read_chunk(path, lineage, contract, generation, suspended=True):
            if event['identity'] in identities:
                raise HistoryError('Duplicate suspended identity')
            identities.add(event['identity'])
            event['rank'] = rank
            result.append(event)
    return sorted(result, key=lambda event: event['identity'])
