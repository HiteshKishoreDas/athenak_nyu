"""Corrupt actual native CRV files; repair checksums to exercise semantic guards."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import struct
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'vis/python'))
import curve_binary as cb
from tracer_thermo_history import read_curves,read_particles


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--historical',type=Path,required=True)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    fixture=sorted((a.runs/'planar/curves').glob('*.crv'))[0]
    data=fixture.read_bytes();state=cb.read(fixture);nfields=len(state['field_names'])
    # Native fixtures are little-endian doubles; offsets below follow documented packed scalars.
    assert data[20:24]==b'\x04\x03\x02\x01' and struct.unpack_from('<I',data,24)[0]==8
    name_at=204+4+struct.unpack_from('<I',data,204)[0]
    name_size=struct.unpack_from('<I',data,name_at)[0]
    curve_at=name_at+4+name_size;node_at=curve_at+100;node_width=128+8*nfields
    results=[]
    def repair(blob):
        blob[-8:]=struct.pack('<Q',cb.checksum(blob[:-8]));return blob
    def mutate(offset,fmt,value):
        blob=bytearray(data);struct.pack_into('<'+fmt,blob,offset,value);return repair(blob)
    def reject(name,blob):
        path=out/(name+'.crv');path.write_bytes(blob)
        try:cb.read(path)
        except (ValueError,OverflowError) as exc:results.append(dict(id=name,status='pass',rejected=str(exc)));return
        raise AssertionError('accepted corrupt '+name)
    reject('truncated_header',data[:20])
    reject('truncated_records',data[:-10])
    blob=bytearray(data);blob[node_at]^=1;reject('checksum_corruption',blob)
    reject('unknown_version',mutate(16,'I',99))
    reject('unknown_real_size',mutate(24,'I',16))
    reject('unknown_endian',mutate(20,'I',0))
    reject('negative_cycle',mutate(32,'i',-1))
    reject('invalid_rank',mutate(36,'i',1))
    reject('invalid_split',mutate(44,'I',31))
    reject('nonfinite_time',mutate(48,'d',float('nan')))
    reject('invalid_domain',mutate(64,'d',0))
    reject('oversized_count',mutate(156,'Q',1000001))
    reject('invalid_link',mutate(node_at+8,'i',0))
    reject('invalid_ancestry',mutate(node_at+12,'i',123))
    reject('future_birth',mutate(node_at+36,'d',1))
    reject('negative_gid',mutate(node_at+120,'i',-1))
    reject('negative_level',mutate(node_at+124,'i',-1))
    reject('invalid_curvature_mask',mutate(node_at+116,'i',7))
    reject('duplicate_node',mutate(node_at+node_width,'i',state['curves'][0]['nodes'][0]['id']))
    index=state['field_names'].index('shock_valid')
    reject('nonboolean_shock_validity',mutate(node_at+128+8*index,'d',2))
    reject('infinite_thermo',mutate(node_at+128,'d',float('inf')))
    reject('trailing_bytes',repair(bytearray(data[:-8]+b'extra'+data[-8:])))
    # Complete-set validation, including fixed metadata across different clocks.
    frames=sorted((a.runs/'planar/curves').glob('*.crv'),key=lambda f:int(f.name.split('.cycle')[1].split('.')[0]))
    def set_reject(name,modify):
        folder=out/name;folder.mkdir()
        for path in frames[:2]:shutil.copy2(path,folder/path.name)
        modify(folder)
        try:cb.read_set(folder)
        except ValueError as exc:
            if name=='missing_rank_chunk':assert 'missing CRV rank chunk' in str(exc)
            results.append(dict(id=name,status='pass',rejected=str(exc)));return
        raise AssertionError('accepted corrupt set '+name)
    def domain_change(folder):
        path=folder/frames[1].name;blob=bytearray(path.read_bytes());struct.pack_into('<d',blob,64,2)
        path.write_bytes(repair(blob))
    set_reject('cross_clock_domain_changed',domain_change)
    def missing_rank(folder):
        for path in list(folder.glob('*.crv')):
            blob=bytearray(path.read_bytes());struct.pack_into('<i',blob,40,2);struct.pack_into('<i',blob,36,1)
            target=folder/path.name.replace('rank0of1','rank1of2');target.write_bytes(repair(blob));path.unlink()
    set_reject('missing_rank_chunk',missing_rank)
    # Exercise -1 and the exhausted 32-bit allocator without passing IDs through Real.
    blob=bytearray(data);struct.pack_into('<I',blob,44,0);struct.pack_into('<i',blob,curve_at,1)
    struct.pack_into('<Q',blob,curve_at+52,1<<32)
    nodes=state['curves'][0]['nodes'];mapping={n['id']:i+1 for i,n in enumerate(nodes)}
    mapping[nodes[-1]['id']]=-1;mapping[0]=0
    for i,node in enumerate(nodes):
        for j,key in enumerate(('id','prev','next','parent_left','parent_right')):
            struct.pack_into('<i',blob,node_at+i*node_width+j*4,mapping[node[key]])
    path=out/'full_width_minus_one.crv';path.write_bytes(repair(blob));full=cb.read(path)
    assert full['curves'][0]['nodes'][-1]['id']==-1 and full['curves'][0]['next_sequence']==1<<32
    results.append(dict(id='minus_one_and_exhausted_allocator_exact',status='pass'))
    # Explicit old/new stream selection must not shadow the requested historical stream.
    mixed=out/'mixed';mixed.mkdir()
    old=next(a.historical.glob('*.thp'));oldstem=old.name.split('.cycle')[0]
    for path in a.historical.glob('*.thp'):shutil.copy2(path,mixed/path.name)
    shutil.copy2(fixture,mixed/fixture.name.replace('.samples.cycle','.other.cycle'))
    before,_=read_curves(a.historical,oldstem);after,_=read_curves(mixed,oldstem)
    assert before.keys()==after.keys()
    results.append(dict(id='historical_stream_not_shadowed',status='pass',records=len(after)))
    # The unchanged particle reader rejects CRV, instead of misinterpreting integer slots.
    try:read_particles(fixture)
    except ValueError:results.append(dict(id='particle_reader_rejects_crv',status='pass'))
    else:raise AssertionError('particle reader accepted CRV')
    # Current material_curves dispatch and independent node audit both accept native file.
    import material_curves
    assert json.dumps(material_curves.read(fixture)['curves'],sort_keys=True)==json.dumps(state['curves'],sort_keys=True)
    results.append(dict(id='existing_geometry_reader_dispatch',status='pass'))
    (out/'report.json').write_text(json.dumps(dict(checks=results,fixture=str(fixture)),indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
