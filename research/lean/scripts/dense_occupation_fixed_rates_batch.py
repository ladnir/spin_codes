"""Replay indexed occupation rate wrappers after their local boxes are checked."""
from pathlib import Path
import hashlib,json,re,subprocess,sys,time
r=Path(__file__).resolve().parents[1];data=r/'scripts/map_data'
entries=json.loads((data/'dense_occupation_fixed_rates_candidates.json').read_text())['entries']
out=data/'dense_occupation_fixed_rates_verification.json'
report={'status':'RUNNING','indexed_rate_boxes_checked':0,'uniform_constant':7000000,'eta':'4/10000000','boxes':{}}
def save():out.write_text(json.dumps(report,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
save()
for e in entries:
 tag=e['box'];record={**e,'status':'WAITING'};report['boxes'][tag]=record;save()
 dep=r/f'.lake/build/lib/lean/SpinCodes/Structured/DenseOccupationFixed{tag}.olean'
 while not dep.exists():
  batch=json.loads((data/'dense_occupation_fixed_all_boxes_verification.json').read_text())
  if batch['status']=='FAIL':raise RuntimeError('Local box replay failed; inspect box manifest')
  time.sleep(15)
 start=time.monotonic();record['status']='RUNNING';save();print(tag,'RATE RUNNING',flush=True)
 name=f'DenseOccupationFixed{tag}Rate';path=f'SpinCodes/Structured/{name}.lean'
 run=subprocess.run([sys.executable,str(r/'scripts/peach-lean.py'),path],cwd=r,text=True,encoding='utf-8',errors='replace',stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 (data/f'dense_occupation_rates_{tag}.log').write_text(run.stdout,encoding='utf-8')
 if run.returncode:
  record.update(status='FAIL',tail=run.stdout[-3000:]);report['status']='FAIL';save();print(tag,'RATE FAIL',flush=True);sys.exit(1)
 a=re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]",run.stdout,re.S)
 assert len(a)==2 and all(set(x.strip() for x in ax.split(','))<={'propext','Classical.choice','Quot.sound'} for _,ax in a)
 record.update(status='PASS',elapsed_seconds=time.monotonic()-start,source_sha256=sha(r/path),olean_sha256=sha(r/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean'),axioms=dict(a))
 report['indexed_rate_boxes_checked']+=1;save();print(tag,'RATE PASS',flush=True)
report['status']='PASS';save();print('PASS all433 indexed occupation rates',flush=True)
