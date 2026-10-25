"""THP v2 clock/validity/corruption regressions; standard library only."""
from pathlib import Path
import struct
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'vis/python'))
from tracer_thermo_history import read_particles,read_curves


def header(names):
    blob='\n'.join(names).encode()
    return struct.pack('=32s4i',b'ATHK_PRTCL_THERMO_HISTORY',2,8,len(names),len(blob))+blob


def block(names,time=0.,cycle=0,values=None):
    count=int(values is not None)
    result=struct.pack('=16s5id',b'ATHKTHPBLK',2,count,4,4+len(names),cycle,time)
    if count:
        result+=struct.pack('=4i',cycle,-1,255,0)
        result+=struct.pack('='+str(4+len(names))+'d',time,.25,.5,.5,*values)
    return result


def reject(action):
    try:action()
    except ValueError:return
    raise AssertionError('corrupt THP accepted')


def main():
    checks=0
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);path=root/'run.thermo.cycle0.rank0of1.bits8.thp'
        names=['birth','a','level','density']
        valid=header(names)+block(names,values=[0,0,3,1])
        path.write_bytes(valid)
        rows,times=read_curves(root)
        assert rows[0,-1]['curve_id']==255 and times=={0:0.}
        checks+=1
        for time,cycle in [(float('inf'),0),(float('nan'),0),(-1.,0),(0.,-1)]:
            path.write_bytes(header(names)+block(names,time,cycle))
            reject(lambda:read_particles(path));reject(lambda:read_curves(root));checks+=1
        path.write_bytes(valid+block(names,time=.1))
        reject(lambda:read_particles(path));reject(lambda:read_curves(root));checks+=1
        path.write_bytes(valid+block(names))
        assert len(read_particles(path)[0])==1  # ordinary THP streams can have repeated clocks
        reject(lambda:read_curves(root))  # one immutable curve chunk holds exactly one block
        checks+=1
        for content in [valid[:-1],valid+b'junk',valid[:32]+struct.pack('=i',9)+valid[36:]]:
            path.write_bytes(content);reject(lambda:read_curves(root));checks+=1
        path.write_bytes(header(names)+block(names))
        wide=root/'run.thermo.cycle0.rank0of4294967295.bits8.thp'
        path.rename(wide);reject(lambda:read_curves(root));wide.rename(path);checks+=1
        invalid=root/'run.thermo.cycle0.rank0of1.bits31.thp'
        path.rename(invalid);reject(lambda:read_curves(root));invalid.rename(path);checks+=1
        names=['birth','a','level','density','shock_valid','shock_region',
            'shock_mach_pressure','shock_mach_temperature','shock_mach_pressure_valid','shock_mach_temperature_valid']
        values=[0,0,3,1,1,0,float('nan'),float('nan'),0,0]
        one=block(names,values=values);path.write_bytes(header(names)+one)
        assert read_curves(root)[0][0,-1]['fields']['shock_mach_pressure_valid']==0
        path.write_bytes(header(names)+one+one)
        rows,_=read_particles(path,allow_identical_replay=True)
        assert len(rows)==1 and rows.identical_replay_rows==1
        checks+=1
        for index,value in [(3,float('nan')),(6,float('inf')),(8,1),(4,0),(5,1)]:
            damaged=values.copy();damaged[index]=value
            # A valid region=1 sample with no estimate is permitted; make parent invalid too.
            if index==5:damaged[4]=0
            if index==4:damaged[8]=1
            path.write_bytes(header(names)+block(names,values=damaged))
            reject(lambda:read_particles(path));checks+=1
        for names,values in [
                (['shock_mach_pressure','shock_mach_pressure_valid'],[float('nan'),1]),
                (['shock_valid','shock_mach_pressure_valid'],[0,1])]:
            path.write_bytes(header(names)+block(names,values=values))
            reject(lambda:read_particles(path));checks+=1
        names=['shock_mach_pressure','shock_mach_pressure_valid']
        path.write_bytes(header(names)+block(names,values=[2,0]))
        assert read_particles(path)[0][0,-1]['fields']['shock_mach_pressure_valid']==0
        checks+=1  # absent parent validity makes a finite estimate with valid=0 legitimate
    print(f'{checks} THP v2 clock, signed-ID, NaN and corruption cases PASS')


if __name__=='__main__':main()
