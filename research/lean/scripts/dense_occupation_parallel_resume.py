"""Retired occupation replay controller, retained for historical inspection.

This controller predates source/object provenance admission. Its existence-only
dependency checks and recovery from old logs must not publish new PASS records.
The completed closure used dense_occupation_admitted_resume.py. For a new replay,
follow scripts/README.md and FINAL_REPRODUCTION.md instead of resuming old lanes.
"""

raise SystemExit(
    "Retired controller: no compilation or report updates were performed. "
    "See scripts/README.md and FINAL_REPRODUCTION.md for current verification."
)

# Historical implementation below; execution is intentionally disabled above.
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import hashlib,json,re,subprocess,sys,time,threading
r=Path(__file__).resolve().parents[1];data=r/'scripts/map_data';mode=sys.argv[1];assert mode in ('boxes','rates')
file='dense_occupation_fixed_all_boxes_verification.json' if mode=='boxes' else 'dense_occupation_fixed_rates_verification.json'
count='numerical_local_boxes_checked' if mode=='boxes' else 'indexed_rate_boxes_checked'
items=json.loads((data/('dense_occupation_fixed_all_boxes_candidates.json' if mode=='boxes' else 'dense_occupation_fixed_rates_candidates.json')).read_text())['boxes' if mode=='boxes' else 'entries']
out=data/file;report=json.loads(out.read_text());lock=threading.Lock()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(path):
 while True:
  try:return json.loads(path.read_text())
  except json.JSONDecodeError:time.sleep(.05)
def audits(text):
 a=re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]",text,re.S)
 if len(a)!=2 or any(set(x.strip() for x in ax.split(','))-{'propext','Classical.choice','Quot.sound'} for _,ax in a):raise RuntimeError('Unexpected axiom audit')
 return dict(a)
def module_record(name):
 return {'name':name,'source_sha256':sha(r/f'SpinCodes/Structured/{name}.lean'),'olean_sha256':sha(r/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean')}
def save():
 report[count]=sum(v['status']=='PASS' for v in report['boxes'].values())
 tmp=out.with_suffix('.resume.tmp');tmp.write_text(json.dumps(report,indent=2)+'\n')
 for attempt in range(100):
  try:tmp.replace(out);return
  except PermissionError:
   if attempt==99:raise
   time.sleep(.02)
# The caller drains both old controllers before starting either new controller.
# A fetched object exists only after peach-lean.py's successful kernel compile and source hash check.
for tag,rec in report['boxes'].items():
 if rec['status']=='PASS':continue
 name='DenseOccupationFixed'+tag+('Rate' if mode=='rates' else '')
 obj=r/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean';log=data/f'peach_{name}.log'
 if obj.exists() and log.exists():
  txt=log.read_text(encoding='utf-8')
  if 'error:' not in txt:
   a=audits(txt)
   if mode=='boxes':rec.update(status='PASS',modules=[module_record(name+'Data'),module_record(name)],axioms=a,recovered_after_controller_drain=True)
   else:
    m=module_record(name);rec.update(status='PASS',source_sha256=m['source_sha256'],olean_sha256=m['olean_sha256'],axioms=a,recovered_after_controller_drain=True)
report.update(status='RUNNING',concurrency=2,resumed=True);save()
def run(item):
 tag=item['box'];start=time.monotonic();rec={**item,'status':'WAITING','modules':[]} if mode=='boxes' else {**item,'status':'WAITING'}
 with lock:report['boxes'][tag]=rec;save()
 try:
  dep=r/f".lake/build/lib/lean/SpinCodes/Structured/DenseOccupationFixed{item['witness'] if mode=='boxes' else tag}.olean"
  while not dep.exists():
   dependency=read(data/('dense_occupation_fixed_remaining_verification.json' if mode=='boxes' else 'dense_occupation_fixed_all_boxes_verification.json'))
   if dependency.get('status') in ('FAIL','FAILED'):raise RuntimeError('Dependency replay failed')
   time.sleep(10)
  with lock:rec.update(status='RUNNING',started=datetime.now(timezone.utc).isoformat());save()
  print(tag,mode,'RUNNING',flush=True)
  for suffix in (['Data',''] if mode=='boxes' else ['Rate']):
   name='DenseOccupationFixed'+tag+suffix;path=f'SpinCodes/Structured/{name}.lean'
   p=subprocess.run([sys.executable,str(r/'scripts/peach-lean.py'),path],cwd=r,text=True,encoding='utf-8',errors='replace',stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
   (data/f'dense_occupation_parallel_{name}.log').write_text(p.stdout,encoding='utf-8')
   if p.returncode:raise RuntimeError(p.stdout[-3500:])
   m=module_record(name)
   with lock:
    if mode=='boxes':rec['modules'].append(m)
    else:rec.update(source_sha256=m['source_sha256'],olean_sha256=m['olean_sha256'])
    if suffix!='Data':rec['axioms']=audits(p.stdout)
    save()
  with lock:rec.update(status='PASS',elapsed_seconds=time.monotonic()-start);save()
  print(tag,mode,'PASS',flush=True);return True
 except Exception as exc:
  with lock:rec.update(status='FAIL',error=str(exc));report['status']='FAIL';save()
  print(tag,mode,'FAIL',str(exc),flush=True);return False
pending=[i for i in items if report['boxes'].get(i['box'],{}).get('status')!='PASS']
with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,pending))
report['status']='PASS' if all(results) and report[count]==433 else 'FAIL';save()
print(report['status'],mode,report[count],flush=True)
