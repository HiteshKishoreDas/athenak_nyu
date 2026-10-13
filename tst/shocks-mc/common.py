"""Small MC fixtures and strict existing thermo-output reader; standard library only."""
import hashlib
import math
from pathlib import Path
import re
import struct
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def edit_input(text, changes):
    remaining = dict(changes)
    block = None
    lines = []
    for line in text.splitlines():
        match = re.fullmatch(r'\s*<([^>]+)>\s*', line)
        if match:
            if block is not None:
                for key in list(remaining):
                    if key.startswith(block + '/'):
                        lines.append(key.split('/', 1)[1] + ' = ' + str(remaining.pop(key)))
            block = match[1]
        assignment = re.match(r'\s*([^#=]+?)\s*=', line)
        if assignment and block:
            key = block + '/' + assignment[1].strip()
            if key in remaining:
                line = assignment[1].strip() + ' = ' + str(remaining.pop(key))
        lines.append(line)
    if block:
        for key in list(remaining):
            if key.startswith(block + '/'):
                lines.append(key.split('/', 1)[1] + ' = ' + str(remaining.pop(key)))
    sections = {}
    for key, value in remaining.items():
        section, parameter = key.split('/', 1)
        sections.setdefault(section, []).append(parameter + ' = ' + str(value))
    for section, parameters in sections.items():
        lines += ['<' + section + '>'] + parameters
    return '\n'.join(lines) + '\n'


def uniform_changes(seed=11, count=4096, blocks=1):
    result = {'mesh/nx1':64, 'meshblock/nx1':64 // blocks,
              'mesh/ix1_bc':'periodic', 'mesh/ox1_bc':'periodic',
              'particles/random_seed':seed, 'tracer_seed_initial/count_per_event':count,
              'tracer_seed_initial/slab_min':.9, 'tracer_seed_initial/slab_max':.99,
              'time/tlim':.25, 'shock_tracking/enabled':'false'}
    for side in ('l', 'r'):
        result.update({'problem/d'+side:1, 'problem/p'+side:1, 'problem/u'+side:1})
    return result


def run(binary, directory, changes=None, restart=None, overrides=(), timeout=30):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=False)
    binary = Path(binary).resolve()
    if restart is None:
        path = directory / 'input.athinput'
        path.write_text(edit_input((ROOT/'inputs/tests/shock_mc.athinput').read_text(),
                                   changes or {}))
        command = [str(binary), '-i', str(path)]
    else:
        path = Path(restart).resolve()
        command = [str(binary), '-r', str(path)] + list(overrides)
    started = time.monotonic()
    with (directory/'run.log').open('w') as log:
        try:
            exit_status = subprocess.run(command, cwd=directory, stdout=log,
                stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            exit_status = 'timeout'
    return {'command':command, 'directory':str(directory), 'exit':exit_status,
            'elapsed_seconds':time.monotonic()-started, 'input_sha256':sha(path),
            'log_sha256':sha(directory/'run.log'), 'executable_sha256':sha(binary)}


class ThermoRows(dict):
    identical_replay_rows = 0


def read_thermo(path, allow_identical_replay=False):
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


def thermo(directory):
    files = list((Path(directory)/'prtcl_thermo_history').glob('*.thp'))
    if len(files) != 1:
        raise ValueError('expected exactly one thermo output file')
    return read_thermo(files[0], allow_identical_replay=True)


def moment_oracle(observations, times, velocity=1, dx=1/64, length=1):
    """Uniform advection: independent Bernoulli cell hops, p_n=u*dt_n/dx.

    The physical flux supplies the probability; this does not reproduce RNG or
    detector code. Six standard errors are predeclared statistical diagnostics.
    """
    cycles = sorted(times); ids = sorted(tag for cycle, tag in observations if cycle == cycles[0])
    probabilities = [velocity*(times[b]-times[a])/dx for a,b in zip(cycles,cycles[1:])]
    if not probabilities or any(p < 0 or p > 1 for p in probabilities):
        raise ValueError('oracle requires resolved positive single-axis advection')
    mu = sum(probabilities); variance = sum(p*(1-p) for p in probabilities)
    fourth = 3*variance*variance + sum(p*(1-p)*(1-6*p*(1-p)) for p in probabilities)
    hops = []
    for tag in ids:
        path = [observations[(c,tag)]['x'][0] for c in cycles]
        displacement = 0
        for a,b in zip(path,path[1:]):
            delta = b-a
            if delta < -length/2:delta += length
            if delta > length/2:delta -= length
            step = delta/dx
            if abs(step-round(step)) > 1e-12 or round(step) not in (0,1):
                raise ValueError('non-cell-hop displacement in uniform MC tracer')
            displacement += round(step)
        hops.append(displacement)
    n = len(hops); mean = sum(hops)/n
    sample_var = sum((value-mean)**2 for value in hops)/(n-1)
    mean_se = math.sqrt(variance/n)
    var_se = math.sqrt(max(0, (fourth - (n-3)/(n-1)*variance*variance)/n))
    return {'particles':n, 'steps':len(probabilities), 'expected_hops':mu,
            'expected_variance':variance, 'observed_hops_mean':mean,
            'observed_hops_variance':sample_var,
            'mean_standard_errors':(mean-mu)/mean_se,
            'variance_standard_errors':(sample_var-variance)/var_se,
            'mean_pass':abs(mean-mu) <= 6*mean_se,
            'variance_pass':abs(sample_var-variance) <= 6*var_se,
            'statistical_assumption':'fixed hash draws tested as independent Bernoulli ensemble',
            'scientific_signoff':'pending; no deterministic per-particle path oracle'}


def check_uniform_transfers(directory, observations, times, seed=11, dx=1/64):
    """TRC11: logged source, bin, physical probability and realized hop evidence."""
    import csv
    seen = set(); maximum_error = 0
    for path in (Path(directory)/'shock_events').glob('transfer*.rank*.csv'):
        match = re.fullmatch(r'transfer(\d+)\.rank(\d+)\.csv',path.name)
        generation,rank = map(int,match.groups())
        lines=path.read_text().splitlines();header=lines[0].split()
        if len(header)!=9 or header[:3]!=['#','AKSHOCK-MC-TRANSFERS','1']:
            raise ValueError('transfer schema/version header mismatch')
        g,count=int(header[4]),int(header[5]);start,end=map(float,header[6:8])
        if g!=generation or int(header[8])!=seed or lines[-1]!='# END '+str(count) or len(lines)!=count+2:
            raise ValueError('transfer count/generation/seed/END mismatch')
        if start!=times[g-1] or end!=times[g]:raise ValueError('transfer timestep bin mismatch')
        expected=(end-start)/dx
        for row in csv.reader(lines[1:-1]):
            if len(row)!=17:raise ValueError('transfer row schema mismatch')
            values=list(map(float,row))
            if not all(map(math.isfinite,values)):raise ValueError('nonfinite transfer evidence')
            if not all(values[i].is_integer() for i in (0,1,2,3,4,16)):
                raise ValueError('noninteger transfer identity/index/direction')
            tag,gid=int(values[0]),int(values[1]);key=(g,tag)
            if key in seen:raise ValueError('duplicate transfer identity')
            seen.add(key)
            before=observations[g-1,tag];after=observations[g,tag]
            if gid!=before['gid'] or tuple(values[12:15])!=before['x']:
                raise ValueError('transfer source parent/position mismatch')
            if values[5]!=1:raise ValueError('uniform start mass mismatch')
            probability=values[6:12];draw=values[15];direction=int(values[16])
            error=max(abs(p-e) for p,e in zip(probability,(0,expected,0,0,0,0)))
            maximum_error=max(maximum_error,error)
            if error>1e-12:raise ValueError('uniform physical probability mismatch')
            if not 0<=draw<1 or direction!=(2 if draw<probability[1] else 0):
                raise ValueError('transfer draw/direction mismatch')
            delta=after['x'][0]-before['x'][0]
            if delta<-.5:delta+=1
            if abs(delta-(dx if direction==2 else 0))>1e-12:
                raise ValueError('realized transfer cell hop mismatch')
    expected_count=sum(1 for g,tag in observations if g>min(times))
    if len(seen)!=expected_count:raise ValueError('missing transfer evidence')
    return {'test_id':'TRC11','records':len(seen),'source_parent_and_timestep_exact':True,
            'max_probability_error':maximum_error,'double_roundoff_tolerance':1e-12,
            'draw_direction_and_realized_cell_hop_consistent':True,
            'scope':'uniform fixedseed serial; conditional source evidence, no continuous physical path'}


def provenance(binary):
    binary=Path(binary).resolve()
    result={'executable':str(binary),'executable_sha256':sha(binary)}
    header=binary.parents[1]/'shock_version.hpp'
    if header.exists():
        match=re.search(r'SHOCK_IMPLEMENTATION_HASH "([0-9a-f]{64})"',header.read_text())
        if match and match[1].encode() in binary.read_bytes():
            result['compiled_source_fingerprint']=match[1]
        else:result['source_fingerprint_status']='build header does not match executable; unbound'
    cache=binary.parents[1]/'CMakeCache.txt'
    if cache.exists():result['build_cache_sha256']=sha(cache)
    return result
