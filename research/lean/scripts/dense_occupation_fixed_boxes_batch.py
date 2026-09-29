from pathlib import Path
import json,subprocess,sys,time,re,hashlib
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'scripts/map_data';OUT=DATA/'dense_occupation_fixed_boxes_verification.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def entry(name):return {'name':name,'source_sha256':sha(ROOT/f'SpinCodes/Structured/{name}.lean'),'olean_sha256':sha(ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean')}
manifest={'status':'RUNNING','scope':'Local affine-segment exponent rectangles plus exact parameter-associated occupation Collatz witness. This does not prove the outer affine majorant or full-domain coverage.','first_box':{'status':'PASS','modules':[entry('DenseOccupationFixedFirstBoxData'),entry('DenseOccupationFixedFirstBox')]},'boxes':{f'W{i:03d}':{'status':'PENDING'} for i in range(1,9)}}
def save():OUT.write_text(json.dumps(manifest,indent=2)+'\n')
def run(i):
 tag=f'W{i:03d}';r=manifest['boxes'][tag];r['status']='RUNNING';save();start=time.monotonic();mods=[]
 print(tag,'BOX RUNNING',flush=True)
 for suffix in ['Data','']:
  if not suffix:
   deadline=time.monotonic()+1800
   while not (ROOT/f'.lake/build/lib/lean/SpinCodes/Structured/DenseOccupationFixed{tag}.olean').exists():
    upstream=json.loads((DATA/'dense_occupation_fixed_batch_verification.json').read_text())['witnesses'][tag]
    if upstream['status']=='FAIL':raise RuntimeError(f'{tag} upstream matrix failed')
    assert time.monotonic()<deadline
    time.sleep(5)
  name=f'DenseOccupationFixed{tag}Box{suffix}'
  p=subprocess.run([sys.executable,str(ROOT/'scripts/peach-lean.py'),f'SpinCodes/Structured/{name}.lean'],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace')
  (DATA/f'dense_occupation_box_{tag}{suffix}.log').write_text(p.stdout,encoding='utf-8')
  if p.returncode:
   r.update(status='FAIL',failed_module=name,output_tail=p.stdout[-2500:]);save();print(tag,'BOX FAIL',p.stdout[-900:],flush=True);return False
  mods.append(entry(name))
  if not suffix:
   a=re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]",p.stdout,re.S)
   assert len(a)==2
   assert all(set(x.strip() for x in ax.split(','))<={'propext','Classical.choice','Quot.sound'} for _,ax in a)
 r.update(status='PASS',elapsed_seconds=time.monotonic()-start,modules=mods,axioms=dict(a));save();print(tag,'BOX PASS',flush=True);return True
save();results=[run(i) for i in range(1,9)]
manifest['status']='PASS' if all(results) else 'FAIL';manifest['numerical_local_boxes_checked']=1+sum(results);save();print(manifest['status'],manifest['numerical_local_boxes_checked'],'local boxes; global coverage unchecked',flush=True)
