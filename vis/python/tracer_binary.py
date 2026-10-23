"""Read native tracer tables: versioned header, typed binary rows, FNV footer.

Mesh .bin and legacy particle .thp remain their existing formats. This related
format preserves uint64 identities and declares its native endian/Real width.
No HDF5 or NumPy dependency. read() provides the common named-table interface.
"""
import math
from pathlib import Path
import re
import struct

MAGIC = b'Athena tracer binary output'
OFFSET, PRIME, MASK = 14695981039346656037, 1099511628211, (1 << 64)-1


def require(condition, message):
    if not condition:
        raise ValueError('tracer binary: '+message)


def read(path, max_records=10_000_000, max_bytes=512*1024*1024):
    path = Path(path)
    require(isinstance(max_records, int) and max_records >= 0 and
            isinstance(max_bytes, int) and max_bytes >= 0, 'invalid reader bounds')
    require(path.stat().st_size <= max_bytes, 'file exceeds reader byte cap')
    data = path.read_bytes()
    require(len(data) >= 84 and data[:32].rstrip(b'\0') == MAGIC, 'unknown/truncated magic')
    marker = data[36:40]
    require(marker in (b'\x78\x56\x34\x12', b'\x12\x34\x56\x78'), 'invalid endian marker')
    endian = '<' if marker[0] == 0x78 else '>'
    version, _, width, header_bytes, count, cycle = struct.unpack_from(endian+'4I2Q', data, 32)
    require(version == 1 and width in (4, 8) and header_bytes <= 65536 and
            count <= max_records, 'unsupported version/width or extent cap')
    real = 'd' if width == 8 else 'f'
    time, rank, ranks = struct.unpack_from(endian+real+'ii', data, 64)
    start = 64+width+8
    require(math.isfinite(time) and time >= 0 and 0 <= rank < ranks and
            start+header_bytes+8 <= len(data), 'invalid clock/rank/header extent')
    try:
        text = data[start:start+header_bytes].decode('utf-8')
    except UnicodeError as exc:
        raise ValueError('tracer binary: invalid header encoding') from exc
    metadata = {}
    for line in text.splitlines():
        key, sep, value = line.partition('=')
        require(sep and key and key not in metadata, 'invalid/duplicate metadata')
        metadata[key] = value
    require('kind' in metadata and 'columns' in metadata, 'missing table schema')
    schema = []
    codes = {'u64':'Q', 'i64':'q', 'u32':'I', 'i32':'i', 'real':real}
    for field in metadata['columns'].split(','):
        name, sep, dtype = field.partition(':')
        require(sep and re.fullmatch(r'[A-Za-z_][A-Za-z0-9_.]*', name) and
                dtype in codes, 'invalid column/type')
        schema.append((name, dtype))
    require(len(set(n for n, _ in schema)) == len(schema), 'duplicate column name')
    record = struct.Struct(endian+''.join(codes[t] for _, t in schema))
    start += header_bytes
    require(start+count*record.size+8 == len(data), 'truncated/trailing data or record extent')
    checksum = OFFSET
    for byte in memoryview(data)[:-8]:
        checksum = ((checksum ^ byte)*PRIME) & MASK
    require(checksum == struct.unpack_from(endian+'Q', data, len(data)-8)[0], 'checksum mismatch')
    names = [name for name, _ in schema]
    rows = [dict(zip(names, row)) for row in record.iter_unpack(memoryview(data)[start:-8])]
    return dict(kind=metadata['kind'], schema=schema, metadata=metadata, rows=rows,
                cycle=cycle, time=time, rank=rank, nranks=ranks, real_bytes=width)


OBSERVATION = [('time','real'),('generation','u64'),('valid','u32'),('inside','u32')]
OBSERVATION += [(f'x{d}','real') for d in (1,2,3)]
OBSERVATION += [(f'image{d}','i64') for d in (1,2,3)]
OBSERVATION += [(n,'real') for n in ('density','pressure','temperature','mach_pressure',
    'mach_temperature','normal1','normal2','normal3','upstream_density',
    'upstream_pressure','upstream_temperature','downstream_density',
    'downstream_pressure','downstream_temperature')]


def shock_rows(path, lineage=None, contract=None, generation=None, *, curve=False,
               suspended=False, trace=False):
    """Validate native schema, then reuse existing semantic shock parsers verbatim."""
    table = read(path)
    carrier = 'curve' if curve else 'particle'
    kind = carrier+('_trace' if trace else '_suspended' if suspended else '_events')
    allowed = (kind, 'curve_birthtrace') if trace and curve else (kind,)
    require(table['kind'] in allowed, 'wrong shock table kind')
    for key, expected in (('lineage',lineage),('contract',contract)):
        require(key in table['metadata'] and
                (expected is None or table['metadata'][key] == expected), key+' mismatch')
    require(generation is None or table['cycle'] == generation, 'generation mismatch')
    schema = [('particle_id','u64'),('sequence','u64'),
              ('phase' if trace else 'status','u32'),('flags','u32')]
    if trace:
        schema += OBSERVATION
    else:
        schema += [(prefix+'.'+n,t) for prefix in ('pre','first','last','post','representative')
                   for n,t in OBSERVATION]
    if curve:
        schema += [('birth','real')]
        if trace:
            schema += [('level','i32')]+[(f'dx{d}','real') for d in (1,2,3)]
        else:
            for p in ('first','last'):
                schema += [(p+'_level','i32')]+[(f'{p}_dx{d}','real') for d in (1,2,3)]
    require(table['schema'] == schema, 'unsupported shock column schema')
    return table, [[str(row[name]) for name,_ in schema] for row in table['rows']]


def event_lines(path, lineage, contract, generation, *, curve=False, suspended=False, expected_ranks=None):
    table, rows = shock_rows(path,lineage,contract,generation,curve=curve,suspended=suspended)
    require(expected_ranks is None or table['nranks']==expected_ranks, 'shock manifest rank count mismatch')
    match = re.fullmatch(r'(?:cycle|suspended)(\d+)\.rank(\d+)\.tbin',Path(path).name)
    require(match and int(match[1]) == table['cycle'] and int(match[2]) == table['rank']
            and table['metadata'].get('position') == 'wrapped', 'shock filename/position mismatch')
    kind = ('AKCURVE-SHOCK-' if curve else 'AKSHOCK-')+('SUSPENDED' if suspended else 'EVENTS')
    return [f'# {kind} 1 {lineage} {generation} {len(rows)}', '# '+contract,
            *(','.join(row) for row in rows), '# END '+str(len(rows))]


GEOMETRY_SCHEMA = [('curve_id','u64'),('closed','i32')]
GEOMETRY_SCHEMA += [(f'winding{d}','i64') for d in (1,2,3)]
GEOMETRY_SCHEMA += [('period','real'),('unresolved_since','real'),('reason','i32')]
GEOMETRY_SCHEMA += [(n,'u64') for n in ('next_sequence','inserted','migrations','samples','node_count')]
GEOMETRY_SCHEMA += [('length','real')]
GEOMETRY_SCHEMA += [(n,'u64') for n in ('node_id','prev','next','parent_left','parent_right')]
GEOMETRY_SCHEMA += [('a','real'),('birth','real'),('anchor','i32'),('placement','i32')]
GEOMETRY_SCHEMA += [(f'x{d}','real') for d in (1,2,3)]
GEOMETRY_SCHEMA += [(f'image{d}','i64') for d in (1,2,3)]
GEOMETRY_SCHEMA += [(n,'real') for n in ('edge_length','stretch','curvature')]+[('curvature_mask','i32')]


def read_geometry(path, max_nodes=1_000_000):
    """Normalize one rank's geometry; apply the existing topology/geometry audit.

    Native analysis snapshots are not checkpoint sidecars and carry no fluid
    fingerprint or AMR clock state. Never pass them to restart loading.
    """
    from material_curves import _audit, _integer, _real, geometry
    table = read(path,max_records=max_nodes)
    require(table['kind']=='curve_geometry' and table['schema']==GEOMETRY_SCHEMA and
            table['metadata'].get('position')=='lifted', 'unsupported geometry schema')
    domain = [_real(v) for v in table['metadata'].get('domain','').split()]
    dimensions = [_integer(v,32) for v in table['metadata'].get('dimensions','').split()]
    require(len(domain)==6 and len(dimensions)==3 and all(d>0 for d in dimensions)
            and all(domain[2*d+1]>domain[2*d] for d in range(3)), 'invalid geometry domain')
    curves, ids = {}, set()
    for row in table['rows']:
        require(all(math.isfinite(v) for k,v in row.items() if k!='curvature'),
                'nonfinite geometry value')
        cid = row['curve_id']
        require(0 < cid < 1<<32 and cid % table['nranks']==table['rank'], 'invalid curve/home rank')
        curve = {n:row[n] for n in ('closed','period','unresolved_since','reason',
                                   'next_sequence','inserted','migrations','samples','length')}
        curve.update(id=cid,winding=[row[f'winding{d}'] for d in (1,2,3)])
        if cid in curves:
            require(all(curves[cid][n]==v for n,v in curve.items()), 'curve metadata changed within snapshot')
        else:
            require(curve['closed'] in (0,1) and curve['period']>0 and
                    0<curve['next_sequence']<=1<<32 and 0<=curve['reason']<=4 and
                    curve['length']>=0, 'invalid curve metadata')
            require((curve['reason']==0 and curve['unresolved_since']==-1) or
                    (curve['reason']!=0 and 0<=curve['unresolved_since']<=table['time']), 'invalid unresolved state')
            require(curve['closed'] or curve['winding']==[0,0,0], 'open winding')
            curves[cid]=dict(curve,nodes=[],expected_count=row['node_count'])
        require(curves[cid]['expected_count']==row['node_count'] and row['node_id'] not in ids,
                'duplicate node or inconsistent node count')
        ids.add(row['node_id'])
        node = {n:row[n] for n in ('prev','next','parent_left','parent_right','a','birth',
                                 'anchor','placement','edge_length','stretch','curvature','curvature_mask')}
        node.update(id=row['node_id'],x=[row[f'x{d}'] for d in (1,2,3)],
                    image=[row[f'image{d}'] for d in (1,2,3)])
        mask, value = node['curvature_mask'], node['curvature']
        require(node['edge_length']>=0 and node['stretch']>=0 and
                ((mask==0 and math.isfinite(value) and value>=0) or
                 (mask in (1,2) and math.isnan(value))), 'invalid geometry diagnostics')
        curves[cid]['nodes'].append(node)
    for curve in curves.values():
        require(len(curve['nodes'])==curve.pop('expected_count') and
                len(curve['nodes']) >= (3 if curve['closed'] else 2), 'invalid curve node count')
        _audit(curve,table['time'],domain,dimensions)
        curve['reconstructed']=geometry(curve,domain)
    return dict(kind='snapshot',cycle=table['cycle'],time=table['time'],domain=domain,
                dimensions=dimensions,curves=sorted(curves.values(),key=lambda c:c['id']),
                curve_count=len(curves),node_count=len(ids),rank=table['rank'],nranks=table['nranks'])
