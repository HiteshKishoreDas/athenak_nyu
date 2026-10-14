"""Final read-only manifest/source-script verification after replay."""
import hashlib
import json
from pathlib import Path
root=Path(__file__).resolve().parent
main=Path('/private/tmp/shock-highres128-long-20261014')
report=json.loads((root/'report.json').read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
assert (root/'shock_events/manifest.rank0').read_bytes()==(main/'shock_events/manifest.rank0').read_bytes()
assert sha(root/'replay_check.py')==report['script_sha256']
result={'manifest_byte_identical':True,'frozen_executed_script_sha256_verified':True,'manifest_sha256':sha(root/'shock_events/manifest.rank0'),'report_sha256':sha(root/'report.json'),'verification_script_sha256':sha(Path(__file__))}
(root/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
