"""Native thermo readers for particles (.thp v2) and material curves (CSV v1/v2).

Both return {(cycle, exact_id): observation}, {cycle: time}; quantities live in
observation['fields']. Curve positions are continuous lifts, particle positions
are wrapped. Populations have separate identity namespaces: never join by ID
alone. No smoothing, resampling or reconstruction of pre-birth histories.
"""
import csv
import math
from pathlib import Path
import re
import struct

from material_curves import _integer, _real, _require, node_curve_id


class ThermoRows(dict):
    identical_replay_rows = 0


def read_particles(path, allow_identical_replay=False):
    """Read exact .thp v2 doubles, rejecting truncation and duplicate identities."""
    data = Path(path).read_bytes(); file_header = struct.Struct('=32s4i')
    if len(data) < file_header.size:
        raise ValueError('truncated thermo file header')
    magic, version, width, fields, namebytes = file_header.unpack_from(data)
    if magic.rstrip(b'\0') != b'ATHK_PRTCL_THERMO_HISTORY' or version != 2 or width != 8:
        raise ValueError('unsupported thermo schema or Real size')
    if fields < 0 or namebytes < 0 or file_header.size+namebytes > len(data):
        raise ValueError('invalid thermo name table')
    names = data[file_header.size:file_header.size+namebytes].decode().split('\n')
    if len(names) != fields:
        raise ValueError('thermo field count mismatch')
    offset = file_header.size+namebytes; result = ThermoRows(); times = {}
    block_header = struct.Struct('=16s5id')  # explicitly serialized, no native padding
    while offset < len(data):
        if offset+block_header.size > len(data):
            raise ValueError('truncated thermo block header')
        magic, version, count, ni, nr, cycle, time_value = block_header.unpack_from(data, offset)
        offset += block_header.size
        if magic.rstrip(b'\0') != b'ATHKTHPBLK' or version != 2 or count < 0 or ni != 4 or nr != 4+fields:
            raise ValueError('invalid thermo block schema')
        end = offset + count*(ni*4 + nr*8)
        if end > len(data):
            raise ValueError('truncated thermo records')
        ints = struct.unpack_from('='+str(count*ni)+'i', data, offset)
        real_offset = offset+count*ni*4
        reals = struct.unpack_from('='+str(count*nr)+'d', data, real_offset)
        for p in range(count):
            c, tag, seed_id, gid = ints[p*ni:(p+1)*ni]
            values = reals[p*nr:(p+1)*nr]
            if c != cycle or values[0] != time_value or not all(map(math.isfinite, values)):
                raise ValueError('invalid thermo observation time or finite value')
            key = (cycle, tag)
            observation = {'cycle':cycle, 'tag':tag, 'seed_id':seed_id, 'gid':gid,
                           'time':time_value, 'x':values[1:4],
                           'fields':dict(zip(names, values[4:]))}
            if key in result:
                if not allow_identical_replay or result[key] != observation:
                    raise ValueError('duplicate thermo observation identity')
                result.identical_replay_rows += 1
            else:
                result[key] = observation
        times[cycle] = time_value; offset = end
    return result, times


def read_curves(directory, stem=None):
    """Read one output stream, validating complete rank sets and exact identities.

    Set stem='basename.output_id' when the directory contains multiple streams.
    Field changes during continuation require a different output id/directory.
    """
    rows, times, ranks, expected_ranks = {}, {}, {}, {}
    identities = {}
    schema = None
    streams = set()
    encoding = None
    paths = sorted(Path(directory).glob('*.csv'))
    for path in paths:
        match = re.fullmatch(r'(.+)\.cycle([0-9]+)\.rank([0-9]+)\.csv', path.name)
        _require(match is not None, 'malformed curve thermo chunk filename')
        if stem is not None and match[1] != stem:
            continue
        streams.add(match[1])
        _require(len(streams) == 1, 'multiple curve thermo streams; select stem')
        lines = path.read_text().splitlines()
        head = lines[0].split() if lines else []
        _require(len(head) == 8 and head[:2] == ['#','ATHK_CURVE_THERMO'] and head[2] in ('1','2'),
                 'unknown/truncated curve thermo header')
        sequence_bits = 24 if head[2]=='2' else 32
        _require(encoding is None or encoding == sequence_bits,'curve thermo ID encoding changed')
        encoding = sequence_bits
        cycle, count, rank, nranks = (_integer(head[n], 32) for n in (3,5,6,7))
        time = _real(head[4])
        _require(time >= 0 and cycle == int(match[2]) and rank == int(match[3])
                 and rank < nranks, 'curve thermo clock/rank mismatch')
        _require(len(lines) == count+3 and lines[-1] == '# END '+str(count),
                 'truncated curve thermo chunk/count mismatch')
        columns = next(csv.reader([lines[1]]))
        fixed = ['node_id','curve_id','birth','a','x1','x2','x3','gid','level']
        _require(columns[:9] == fixed and len(columns)>9 and
                 len(set(columns)) == len(columns), 'curve thermo column mismatch')
        if schema is None:
            schema = columns
        _require(schema == columns, 'curve thermo field schema changed')
        _require(cycle not in times or times[cycle] == time, 'rank clocks disagree')
        _require(cycle not in expected_ranks or expected_ranks[cycle] == nranks,
                 'rank counts disagree')
        times[cycle], expected_ranks[cycle] = time, nranks
        _require(rank not in ranks.setdefault(cycle,set()), 'duplicate rank chunk')
        ranks[cycle].add(rank)
        for values in csv.reader(lines[2:-1]):
            _require(len(values) == len(columns), 'curve thermo row width mismatch')
            node, curve = _integer(values[0]), _integer(values[1],32)
            birth, a = _real(values[2]), _real(values[3])
            _require(curve > 0 and node_curve_id(node,sequence_bits) == curve
                     and 0 <= birth <= time and a >= 0, 'invalid node identity/birth/label')
            _require(curve % nranks == rank, 'curve thermo home rank mismatch')
            identity = curve, birth, a
            _require(node not in identities or identities[node] == identity,
                     'curve thermo immutable node metadata changed')
            identities[node] = identity
            key = cycle, node
            _require(key not in rows, 'duplicate curve thermo observation')
            rows[key] = dict(cycle=cycle,node_id=node,curve_id=curve,birth=birth,a=a,
                time=time,x=tuple(_real(v) for v in values[4:7]),
                gid=_integer(values[7],32),level=_integer(values[8],32),
                fields={k:_real(v) for k,v in zip(columns[9:],values[9:])})
    _require(times, 'no curve thermo chunks')
    _require(all(ranks[c] == set(range(expected_ranks[c])) for c in times),
             'missing curve thermo rank chunk')
    ordered = sorted(times)
    _require(all(times[a] < times[b] for a,b in zip(ordered,ordered[1:])),
             'curve thermo time did not increase')
    return rows, times
