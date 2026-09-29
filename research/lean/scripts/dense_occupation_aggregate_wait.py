"""Compile the occupation aggregate only after all indexed certificates pass."""
from pathlib import Path
import hashlib,json,re,subprocess,sys,time
r=Path(__file__).resolve().parents[1];data=r/'scripts/map_data'
out=data/'dense_occupation_aggregate_verification.json'
def save(x):
 tmp=out.with_suffix('.aggregate.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(out)
save({'status':'WAITING','required_indexed_boxes':433})
while True:
 try:report=json.loads((data/'dense_occupation_fixed_rates_verification.json').read_text())
 except json.JSONDecodeError:time.sleep(.1);continue
 if report['status']=='FAIL':save({'status':'BLOCKED','reason':'An indexed occupation rate failed.'});sys.exit(1)
 if report['status']=='PASS':
  assert len(report['boxes'])==433
  for tag,rec in report['boxes'].items():
   assert rec['status']=='PASS'
   name='DenseOccupationFixed'+tag+'Rate'
   assert hashlib.sha256((r/f'SpinCodes/Structured/{name}.lean').read_bytes()).hexdigest()==rec['source_sha256']
   assert hashlib.sha256((r/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean').read_bytes()).hexdigest()==rec['olean_sha256']
  break
 time.sleep(15)
save({'status':'RUNNING','required_indexed_boxes':433})
name='DenseOccupationAllCertified';path=f'SpinCodes/Structured/{name}.lean'
run=subprocess.run([sys.executable,str(r/'scripts/peach-lean.py'),path],cwd=r,text=True,encoding='utf-8',errors='replace',stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
(data/'dense_occupation_aggregate.log').write_text(run.stdout,encoding='utf-8')
if run.returncode:
 save({'status':'FAIL','tail':run.stdout[-4000:]});print('Aggregate FAIL',flush=True);sys.exit(1)
a=re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]",run.stdout,re.S)
assert len(a)==1 and all(set(x.strip() for x in ax.split(','))<={'propext','Classical.choice','Quot.sound'} for _,ax in a)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
save({'status':'PASS','indexed_rate_boxes_checked':433,'axioms':dict(a),'source_sha256':sha(r/path),'olean_sha256':sha(r/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean')})
print('PASS occupation aggregate: all433 indexed boxes',flush=True)
