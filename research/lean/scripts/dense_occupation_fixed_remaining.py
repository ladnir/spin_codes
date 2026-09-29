"""Replay remaining 124 witnesses, two independent Lean processes at most."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime,timezone
import subprocess,sys,time,json,hashlib,re,threading
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'scripts/map_data'
OUT=DATA/'dense_occupation_fixed_remaining_verification.json'
lock=threading.Lock()
manifest={'status':'RUNNING','concurrency':2,'per_process_memory_cap_mb':20000,'numerical_exponent_boxes_checked':0,'scope':'Numerical occupation matrix witnesses and actual routed probability; no vertex exponent or cover claim. Checked cached dependencies reused, no full dependency replay.','witnesses':{f'W{i:03d}':{'status':'PENDING','index':i} for i in range(9,133)}}
def save():OUT.write_text(json.dumps(manifest,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(i):
 tag=f'W{i:03d}';start=time.monotonic();record=manifest['witnesses'][tag]
 with lock:record.update(status='RUNNING',started=datetime.now(timezone.utc).isoformat());save()
 print(tag,'RUNNING',flush=True)
 modules=[]
 for suffix in ['Data','']:
  name=f'DenseOccupationFixed{tag}{suffix}'
  path=f'SpinCodes/Structured/{name}.lean'
  res=subprocess.run([sys.executable,str(ROOT/'scripts/peach-lean.py'),path],cwd=ROOT,text=True,encoding='utf-8',errors='replace',stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  (DATA/f'dense_occupation_batch_{tag}{suffix}.log').write_text(res.stdout,encoding='utf-8')
  if res.returncode:
   with lock:record.update(status='FAIL',failed_module=name,elapsed_seconds=time.monotonic()-start,output_tail=res.stdout[-2500:]);save()
   print(tag,'FAIL',name,res.stdout[-700:],flush=True);return False
  modules.append({'name':name,'source_sha256':sha(ROOT/path),'olean_sha256':sha(ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean')})
  if not suffix:
   a=re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]",res.stdout,re.S)
   assert len(a)==3,(tag,a)
   assert all(set(x.strip() for x in ax.split(','))<={'propext','Classical.choice','Quot.sound'} for _,ax in a)
 with lock:
  record.update(status='PASS',elapsed_seconds=time.monotonic()-start,modules=modules,axioms=dict(a),candidate_sha256=sha(DATA/f'dense_occupation_fixed_{tag}_candidate.json'));save()
 print(tag,'PASS',round(record['elapsed_seconds'],1),'seconds',flush=True);return True
save()
with ThreadPoolExecutor(max_workers=2) as pool:
 results=list(pool.map(run,range(9,133)))
manifest['status']='PASS' if all(results) else 'FAIL'
manifest['numerical_witnesses_checked']=sum(r['status']=='PASS' for r in manifest['witnesses'].values())
save();print(manifest['status'],manifest['numerical_witnesses_checked'],'witnesses; 0 exponent boxes',flush=True)
