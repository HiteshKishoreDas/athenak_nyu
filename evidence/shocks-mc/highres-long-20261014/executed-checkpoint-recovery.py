"""Dump exact checkpoint fluid states with native restart; evolve zero steps."""
from pathlib import Path
import json,shutil,struct,subprocess,sys,time
root=Path('/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu');sys.path.insert(0,str(root/'tst/shocks-mc'));from common import sha,provenance
source=Path('/private/tmp/shock-highres128-long-20261014');cp=Path(sys.argv[1]).resolve();side=Path(str(cp)+'.shock.rank0')
# ponytail: current same-binary local ABI only; the native loader validates the
# complete paired checksums/config/state, rather than duplicating its parser.
with side.open('rb') as f:head=f.read(144)
assert head[:8]==b'AKSHOCK2' and struct.unpack_from('=5I',head,8)==(2,0x12345678,8,816,4)
generation,count=struct.unpack_from('=2Q',head,56);t=struct.unpack_from('=d',head,136)[0];assert count==8192
output=source/'checkpoint-views'/('generation-'+str(generation));output.mkdir(parents=True,exist_ok=False)
history=output/'shock_events';history.mkdir();shutil.copy2(source/'shock_events/manifest.rank0',history/'manifest.rank0')
for g in range(generation+1):shutil.copy2(source/'shock_events'/f'cycle{g}.rank0.csv',history/f'cycle{g}.rank0.csv')
binary=root/'build-release/src/athena';command=[str(binary),'-r',str(cp),f'time/nlim={generation}','output4/dt=-1','output3/dcycle=1','output5/dcycle=0','output6/dcycle=1','job/basename=checkpoint_view',f'shock_tracking/directory={history}','shock_tracking/trace=false']
start=time.monotonic()
with (output/'run.log').open('w') as f:code=subprocess.run(command,cwd=output,stdout=f,stderr=subprocess.STDOUT,timeout=180).returncode
assert code==0,(output/'run.log').read_text()
record={**provenance(binary),'checkpoint':str(cp),'checkpoint_sha256':sha(cp),'sidecar_sha256':sha(side),'generation':generation,'time':t,'evolved_steps':0,'command':command,'exit':code,'seconds':time.monotonic()-start,'log_sha256':sha(output/'run.log')}
(output/'recovery.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
