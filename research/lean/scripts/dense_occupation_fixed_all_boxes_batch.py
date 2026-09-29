"""Kernel replay for the 433 generated occupation boxes, one lightweight worker."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,re,subprocess,sys,time
r=Path(__file__).resolve().parents[1]
data=r/'scripts/map_data'
out=data/'dense_occupation_fixed_all_boxes_verification.json'
candidates=json.loads((data/'dense_occupation_fixed_all_boxes_candidates.json').read_text())['boxes']
report={'status':'RUNNING','scope':'Local exponent rectangles with actual matrix association. Global cover is a separate theorem.','numerical_local_boxes_checked':0,'boxes':{}}
def save():out.write_text(json.dumps(report,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
save()
for b in candidates:
 tag=b['box'];record={'status':'RUNNING','witness':b['witness'],'segment':b['segment'],'started':datetime.now(timezone.utc).isoformat(),'modules':[]}
 report['boxes'][tag]=record;save();start=time.monotonic()
 print(tag,'RUNNING',flush=True)
 for suffix in ['Data','']:
  if not suffix:
   dep=r/f".lake/build/lib/lean/SpinCodes/Structured/DenseOccupationFixed{b['witness']}.olean"
   while not dep.exists():
    batch=json.loads((data/'dense_occupation_fixed_remaining_verification.json').read_text())
    state=batch['witnesses'].get(b['witness'],{}).get('status')
    if state=='FAIL':raise RuntimeError(f"Witness {b['witness']} failed; cannot certify {tag}")
    time.sleep(15)
  name=f'DenseOccupationFixed{tag}{suffix}';path=f'SpinCodes/Structured/{name}.lean'
  run=subprocess.run([sys.executable,str(r/'scripts/peach-lean.py'),path],cwd=r,text=True,encoding='utf-8',errors='replace',stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  (data/f'dense_occupation_all_boxes_{tag}{suffix}.log').write_text(run.stdout,encoding='utf-8')
  if run.returncode:
   record.update(status='FAIL',failed_module=name,tail=run.stdout[-3000:]);report['status']='FAIL';save();print(tag,'FAIL',flush=True);sys.exit(1)
  record['modules'].append({'name':name,'source_sha256':sha(r/path),'olean_sha256':sha(r/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean')})
  if not suffix:
   a=re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]",run.stdout,re.S)
   assert len(a)==2 and all(set(x.strip() for x in ax.split(','))<={'propext','Classical.choice','Quot.sound'} for _,ax in a)
   record['axioms']=dict(a)
 record.update(status='PASS',elapsed_seconds=time.monotonic()-start);report['numerical_local_boxes_checked']+=1;save();print(tag,'PASS',flush=True)
report['status']='PASS';save();print('PASS',report['numerical_local_boxes_checked'],'local boxes; no global coverage claim',flush=True)
