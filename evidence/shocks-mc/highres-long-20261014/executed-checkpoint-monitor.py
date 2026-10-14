from pathlib import Path
import json,subprocess,time
run=Path('/private/tmp/shock-highres128-long-20261014')
for index in (2,3):
 checkpoint=run/'rst'/f'shock_highres128_long.{index:05d}.rst'
 while not Path(str(checkpoint)+'.shock.rank0').exists():
  state=json.loads((run/'execution.json').read_text());assert state['status']=='running',state
  time.sleep(5)
 with (run/f'recovery-{index}.log').open('w') as log:
  code=subprocess.run(['/Users/hitesh/hitesh/git/shockit/.venv/bin/python','/private/tmp/recover-highres-checkpoint.py',str(checkpoint)],stdout=log,stderr=subprocess.STDOUT,timeout=180).returncode
 assert code==0,checkpoint
 print(json.dumps({'recovered':str(checkpoint)}),flush=True)
