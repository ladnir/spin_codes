"""Assemble already checked scalar/Fourier families, then the common dense theorem."""
from pathlib import Path
import hashlib,json,re,subprocess,sys,time
r=Path(__file__).resolve().parents[1];data=r/'scripts/map_data'
out=data/'dense_mixed_aggregate_verification.json'
state=json.loads(out.read_text()) if out.exists() else {'status':'WAITING','eta':'4/10000000','common_constant':24000000000000,'groups':{}}
def save():
 tmp=out.with_suffix('.aggregate.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(out)
def read(n):
 try:return json.loads((data/n).read_text())
 except (FileNotFoundError,json.JSONDecodeError):return None
def compile(name,count):
 state['status']='RUNNING';save();print(name,'RUNNING',flush=True)
 path=f'SpinCodes/Structured/{name}.lean'
 run=subprocess.run([sys.executable,str(r/'scripts/peach-lean.py'),path],cwd=r,text=True,encoding='utf-8',errors='replace',stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 (data/(name+'_aggregate.log')).write_text(run.stdout,encoding='utf-8')
 if run.returncode:
  state.update(status='FAIL',failed_module=name,tail=run.stdout[-4000:]);save();print(name,'FAIL',flush=True);sys.exit(1)
 a=re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]",run.stdout,re.S)
 assert len(a)==1 and all(set(x.strip() for x in ax.split(','))<={'propext','Classical.choice','Quot.sound'} for _,ax in a),(name,a)
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 state['groups'][name]={'status':'PASS','indexed_boxes':count,'axioms':dict(a),'source_sha256':sha(r/path),'object_sha256':sha(r/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean')};save();print(name,'PASS',flush=True)
def valid_group(name):
 rec=state['groups'].get(name,{})
 src=r/f'SpinCodes/Structured/{name}.lean';obj=r/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean'
 return rec.get('status')=='PASS' and src.exists() and obj.exists() and hashlib.sha256(src.read_bytes()).hexdigest()==rec.get('source_sha256') and hashlib.sha256(obj.read_bytes()).hexdigest()==rec.get('object_sha256')
for name in list(state['groups']):
 if not valid_group(name):del state['groups'][name]
state.update(status='RUNNING');state.pop('failed_module',None);state.pop('tail',None)
save()
while 'DenseOccupationScalarCertified' not in state['groups']:
 s=read('dense_scalar_geometry_verification.json')
 if s and s.get('status')=='PASS':
  assert len(s['boxes'])==283 and all(x['status']=='PASS' for x in s['boxes'].values())
  break
 time.sleep(15)
if 'DenseOccupationScalarCertified' not in state['groups']:compile('DenseOccupationScalarCertified',283)
while 'DenseOccupationFourierCertified' not in state['groups']:
 reports=[read(n) for n in ['dense_fourier_box_replay.json','dense_fourier_box_second_replay.json','dense_fourier_box_prefix_tail_replay.json','dense_fourier_box_suffix_tail_replay.json']]
 if all(report and report.get('status')=='PASS' for report in reports) and (read('dense_fourier_complete_verification.json') or {}).get('status')=='PASS':
  rows=[x for report in reports for x in report.get('boxes',[])]
  if {x['index'] for x in rows}==set(range(307)):
   assert len(rows)==307
   break
 time.sleep(15)
if 'DenseOccupationFourierCertified' not in state['groups']:compile('DenseOccupationFourierCertified',307)
while True:
 a=read('dense_occupation_aggregate_verification.json')
 if a and a.get('status')=='PASS':
  name='DenseOccupationAllCertified'
  assert hashlib.sha256((r/f'SpinCodes/Structured/{name}.lean').read_bytes()).hexdigest()==a['source_sha256']
  assert hashlib.sha256((r/f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean').read_bytes()).hexdigest()==a['olean_sha256']
  break
 if a and a.get('status') in {'FAIL','BLOCKED'}:
  state.update(status='BLOCKED',reason='Occupation aggregate not checked');save();sys.exit(1)
 time.sleep(15)
if 'DenseOccupationMixedCertified' not in state['groups']:compile('DenseOccupationMixedCertified',1023)
state['status']='PASS';save();print('PASS all1023 indexed numerical boxes, mixed dense rate theorem',flush=True)
