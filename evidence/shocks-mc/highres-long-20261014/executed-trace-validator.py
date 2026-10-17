"""Stream-validate immutable full traces without keeping every particle in memory."""
import argparse,csv,hashlib,io,json,math,re,time
from pathlib import Path
import sys
root=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu');sys.path.insert(0,str(root/'vis/python'))
from shock_history import integer,parse_observation
p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run=a.run.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
selected={0,2487,4551,4575};last_sequence=[0]*8192;g=0;previous_time=-1;digests={};total=0
with (out/'frames.jsonl').open('w') as frames,(out/'selected-traces.csv').open('w',newline='') as selected_stream:
 writer=csv.writer(selected_stream);writer.writerow(['particle_id','sequence','phase','flags','time','generation','valid','inside','x1','x2','x3','image1','image2','image3','rho','pressure','temperature','mach_pressure','mach_temperature','normal1','normal2','normal3','rho_up','pressure_up','temperature_up','rho_down','pressure_down','temperature_down'])
 while True:
  path=run/'shock_events'/f'trace{g}.rank0.csv'
  if not path.exists():
   execution=json.loads((run/'execution.json').read_text())
   if execution['status']!='running':
    assert execution['status']=='completed' and execution['exit']==0,execution
    assert previous_time==1.0,previous_time
    break
   time.sleep(.25);continue
  data=path.read_bytes();rows=csv.reader(io.StringIO(data.decode()));seen=set();inside=valid=0;phases=[0,0,0];frame_time=None
  for row in rows:
   assert len(row)==28
   pid=integer(row[0]);sequence=integer(row[1]);phase=integer(row[2],32);flags=integer(row[3],32)
   assert pid<8192 and pid not in seen and phase in (0,1,2) and not flags&~1023
   assert last_sequence[pid]<=sequence<=last_sequence[pid]+1;last_sequence[pid]=sequence;seen.add(pid)
   observation=parse_observation(row[4:]);assert observation['generation']==g
   if frame_time is None:frame_time=observation['time']
   assert frame_time==observation['time']
   assert all(0<=v<1 for v in observation['x'])
   assert all(abs(v*128-(math.floor(v*128)+.5))<=1e-12 for v in observation['x'])
   inside+=observation['inside'];valid+=observation['valid'];phases[phase]+=1
   if pid in selected:writer.writerow(row)
  assert seen==set(range(8192)) and frame_time>previous_time
  previous_time=frame_time;total+=8192;digests[path.name]=hashlib.sha256(data).hexdigest()
  frames.write(json.dumps({'generation':g,'time':frame_time,'particles':8192,'valid':valid,'inside':inside,'phase_outside':phases[0],'phase_active':phases[1],'phase_pending':phases[2]})+'\n');frames.flush();selected_stream.flush()
  (out/'progress.json').write_text(json.dumps({'status':'collecting','last_generation':g,'last_time':frame_time,'validated_observations':total,'last_inside':inside,'last_active':phases[1],'last_pending':phases[2]},indent=2)+'\n')
  g+=1
(out/'trace-sha256.json').write_text(json.dumps(digests,indent=2)+'\n')
(out/'progress.json').write_text(json.dumps({'status':'completed','last_generation':g-1,'last_time':previous_time,'validated_observations':total,'all_8192_ids_every_frame':True,'continuous_generations_and_sequences':True,'full_precision_clocks':True,'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
print(json.dumps({'status':'completed','frames':g,'observations':total,'output':str(out)}))
