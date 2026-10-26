"""Typed, checksummed AthenaK .crv v1 snapshots; no particle THP changes.

Native scalar byte order is declared by the endian marker; structs have no padding.
Each rank owns whole curves. read_set validates complete rank sets before exposing
samples. Historical ASCII geometry/checkpoints remain material_curves.read().
"""
import math
from pathlib import Path
import re
import struct

from material_curves import OFFSET, PRIME, MASK64, _audit, _require, geometry


def checksum(data):
    value=OFFSET
    for byte in data:value=((value^byte)*PRIME)&MASK64
    return value


def read(path,max_nodes=1000000):
    _require(type(max_nodes) is int and max_nodes>=0,'invalid CRV reader cap')
    path=Path(path)
    # Bound allocation before reading an untrusted file; field count also bounded below.
    _require(path.stat().st_size<=20000000+max_nodes*8400,'CRV file exceeds reader cap')
    data=path.read_bytes()
    _require(len(data)>=48 and data[:16]==b'ATHKCRV'+b'\0'*9,'unknown/truncated CRV header')
    endian='<' if data[20:24]==b'\x04\x03\x02\x01' else '>' if data[20:24]==b'\x01\x02\x03\x04' else None
    _require(endian is not None,'invalid CRV endian marker')
    _require(checksum(data[:-8])==struct.unpack_from(endian+'Q',data,len(data)-8)[0],
             'CRV checksum mismatch')
    offset=16
    def take(fmt):
        nonlocal offset
        parser=struct.Struct(endian+fmt)
        _require(offset+parser.size<=len(data)-8,'truncated CRV records')
        values=parser.unpack_from(data,offset);offset+=parser.size
        return values[0] if len(values)==1 else list(values)
    version,marker,width,nfields=take('4I')
    _require(version==1 and width in (4,8) and 0<nfields<=1024,'unsupported CRV schema/fields/Real width')
    real='d' if width==8 else 'f'
    cycle,rank,nranks,bits=take('3iI')
    time=take(real);domain=take(6*real);dimensions=take('3i');bcs=take('6i')
    cap,ncurves,nnodes=take('3Q')
    config=dict(zip(('max_spacing','spacing_cells','max_angle','chord_error'),take(4*real)))
    config.update(max_nodes=cap,adaptive=take('i'),output_dcycle=take('i'))
    def text(cap):
        nonlocal offset
        count=take('I')
        _require(count<=cap and offset+count<=len(data)-8,'invalid CRV string extent')
        try:value=data[offset:offset+count].decode('utf-8')
        except UnicodeDecodeError as exc:raise ValueError('invalid CRV text encoding') from exc
        offset+=count
        return value
    basename=text(4096);names=text(65536).split('\n')
    _require(len(names)==nfields and len(set(names))==len(names) and
             all(re.fullmatch('[A-Za-z][A-Za-z0-9_]*',n) for n in names),'invalid CRV field names')
    _require(cycle>=0 and math.isfinite(time) and time>=0 and 0<=rank<nranks and
             0<=bits<=30 and 0<cap and nnodes<=min(cap,max_nodes) and
             ncurves<=min(100000,(1<<bits)-1 if bits else 1),'invalid CRV clock/rank/count/split')
    _require(all(math.isfinite(x) for x in domain) and
             all(domain[2*d]<domain[2*d+1] and dimensions[d]>0 for d in range(3)),
             'invalid CRV domain/dimensions')
    _require(all(math.isfinite(config[k]) for k in ('max_spacing','spacing_cells','max_angle','chord_error')) and
             config['max_spacing']>=0 and config['spacing_cells']>0 and 0<config['max_angle']<math.pi and
             config['chord_error']>=0 and config['adaptive'] in (0,1) and config['output_dcycle']>=0,
             'invalid CRV configuration')
    state=dict(kind='snapshot',version=5,binary_version=1,precision=8*width,
        cycle=cycle,time=time,rank=rank,nranks=nranks,curve_count=ncurves,node_count=nnodes,
        id_layout=(32,bits,True),domain=domain,dimensions=dimensions,bcs=bcs,
        config=config,basename=basename,field_names=names,curves=[])
    from tracer_thermo_history import _sample_fields
    ids=set();total=0
    for _ in range(ncurves):
        cid,closed,reason=take('3i');winding=take('3q');period,status=take(2*real)
        sequence,inserted,migrations,samples,count=take('5Q');length=take(real)
        _require(cid>0 and cid<=((1<<bits)-1 if bits else 1) and cid not in ids and cid%nranks==rank and
                 closed in (0,1) and 0<=reason<=4 and math.isfinite(period) and period>0 and
                 0<sequence<=1<<(32-bits) and math.isfinite(length) and length>=0 and
                 ((reason==0 and status==-1) or (reason!=0 and 0<=status<=time)) and
                 (closed or winding==[0,0,0]) and
                 (3 if closed else 2)<=count<=nnodes-total,'invalid CRV curve metadata/count')
        ids.add(cid)
        curve=dict(id=cid,closed=closed,reason=reason,winding=winding,period=period,
            unresolved_since=status,next_sequence=sequence,inserted=inserted,
            migrations=migrations,samples=samples,length=length,nodes=[])
        for _ in range(count):
            node=dict(zip(('id','prev','next','parent_left','parent_right','anchor','placement'),take('7i')))
            a,birth,*x=take(5*real);image=take('3q')
            edge,stretch,curvature=take(3*real);mask,gid,level=take('3i')
            _require(all(math.isfinite(v) for v in [a,birth,*x,edge,stretch]) and edge>=0 and stretch>=0 and
                     gid>=0 and level>=0 and all(abs(v)<=2**40 for v in image) and
                     ((mask==0 and math.isfinite(curvature) and curvature>=0) or
                      (mask in (1,2) and math.isnan(curvature))),'invalid CRV node diagnostics/location')
            raw=take(nfields*real);raw=[raw] if nfields==1 else raw
            node.update(a=a,birth=birth,x=x,image=image,edge_length=edge,stretch=stretch,
                        curvature=curvature,curvature_mask=mask,gid=gid,level=level,
                        fields=_sample_fields(names,raw))
            curve['nodes'].append(node)
        _audit(curve,time,domain,dimensions,32-bits,bits,True)
        curve['reconstructed']=geometry(curve,domain)
        state['curves'].append(curve);total+=count
    _require(total==nnodes and offset==len(data)-8,'CRV total count/trailing data mismatch')
    return state


def read_set(directory,stem=None,max_nodes=1000000):
    """Complete rank sets per stream/clock; no synthetic samples or times."""
    groups={};streams=set()
    for path in sorted(Path(directory).glob('*.crv')):
        match=re.fullmatch(r'(.+)\.cycle([0-9]+)\.rank([0-9]+)of([0-9]+)\.crv',path.name)
        _require(match is not None,'malformed CRV filename')
        if stem is not None and match[1]!=stem:continue
        state=read(path,max_nodes)
        streams.add(match[1]);_require(len(streams)==1,'multiple curve streams; select thermo_stem')
        _require((state['cycle'],state['rank'],state['nranks'])==tuple(map(int,match.group(2,3,4))),
                 'CRV filename clock/rank mismatch')
        group=groups.setdefault(state['cycle'],[])
        _require(all(s['rank']!=state['rank'] for s in group),'duplicate CRV rank chunk')
        group.append(state)
    frames=[];previous=-1;immutable={};schema=None;layout=None;fixed=None;curve_fixed={};counters={};seen=set()
    for cycle,group in sorted(groups.items()):
        first=group[0]
        _require(len(group)==first['nranks'],'missing CRV rank chunk')
        for state in group:
            for key in ('time','nranks','id_layout','precision','domain','dimensions','bcs','config','basename','field_names'):
                _require(state[key]==first[key],'CRV rank metadata differs: '+key)
        _require(first['time']>previous,'CRV time did not increase');previous=first['time']
        _require(schema is None or schema==first['field_names'],'CRV field schema changed')
        _require(layout is None or layout==first['id_layout'],'CRV ID layout changed')
        schema=first['field_names'];layout=first['id_layout']
        constants={k:first[k] for k in ('domain','dimensions','bcs','precision','basename')}
        constants['geometry_config']={k:first['config'][k] for k in
            ('max_spacing','spacing_cells','max_angle','chord_error','adaptive')}
        _require(fixed is None or fixed==constants,'CRV fixed stream metadata changed')
        fixed=constants

        curves=[c for s in group for c in s['curves']]
        total=sum(len(c['nodes']) for c in curves)
        _require(total<=min(max_nodes,first['config']['max_nodes']),'global CRV node cap exceeded')
        current=set()
        _require(len(curves)<=100000,'global CRV curve cap exceeded')
        _require(not curve_fixed or set(curve_fixed)=={c['id'] for c in curves},'CRV curve set changed')
        for curve in curves:
            cid=curve['id']
            constants=(curve['closed'],curve['period'],tuple(curve['winding']))
            _require(cid not in curve_fixed or curve_fixed[cid]==constants,'CRV fixed curve metadata changed')
            curve_fixed[cid]=constants
            values=tuple(curve[k] for k in ('next_sequence','inserted','migrations','samples'))
            _require(cid not in counters or all(a<=b for a,b in zip(counters[cid],values)),
                     'CRV counters decreased')
            counters[cid]=values
            for node in curve['nodes']:
                current.add(node['id'])
                identity=(curve['id'],node['birth'],node['a'],node['parent_left'],node['parent_right'],node['anchor'],node['placement'])
                _require(node['id'] not in immutable or immutable[node['id']]==identity,'CRV immutable node metadata changed')
                immutable[node['id']]=identity
        _require(seen<=current,'CRV existing nodes disappeared')
        seen=current
        frames.append(dict(first,curves=sorted(curves,key=lambda c:c['id']),node_count=total,
                           curve_count=len(curves),path=str(Path(directory)),stem=next(iter(streams))))
    return frames


def samples(frames):
    from tracer_thermo_history import ThermoRows
    rows=ThermoRows();times={};rows.id_layout=None;rows.field_names=None;rows.stem=None
    rows.geometry=frames
    for frame in frames:
        cycle=frame['cycle'];times[cycle]=frame['time']
        rows.id_layout=frame['id_layout'];rows.field_names=frame['field_names'];rows.stem=frame['stem']
        for curve in frame['curves']:
            for node in curve['nodes']:
                rows[cycle,node['id']]=dict(cycle=cycle,node_id=node['id'],curve_id=curve['id'],
                    birth=node['birth'],a=node['a'],time=frame['time'],x=tuple(node['x']),
                    gid=node['gid'],level=node['level'],fields=node['fields'])
    return rows,times
