"""Pipeline Fourier box checks behind checked matrix witnesses; no benchmarks."""
from pathlib import Path
import subprocess,sys,json,hashlib,re,time,os,argparse
os.environ["FOURIER_SECOND_HALF"]="1"
R=Path(__file__).resolve().parents[1];D=R/'scripts/map_data'
parser=argparse.ArgumentParser();parser.add_argument('lane');args=parser.parse_args()
config=json.loads((D/'dense_fourier_box_lanes.json').read_text())[args.lane]
os.environ['FOURIER_BOX_LANE']=args.lane
sys.stdout.reconfigure(encoding='utf-8')
source=json.loads((R.parent/'workstreams/inner_design/imt_asymptotic/d11/DENSE_REPLAY.json').read_text())
unique={};boxes=[]
for gi,b in enumerate(source['leaves']):
 w=b['rational_witness']
 if w['family']!='fourier':continue
 key=json.dumps(w,sort_keys=True)
 if key not in unique:unique[key]=len(unique)
 boxes.append((gi,unique[key]))
progress_path=D/config['manifest']
rows=json.loads(progress_path.read_text()).get('boxes',[]) if progress_path.exists() else []
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for i in config["indices"]:
 gi,wi=boxes[i]
 tag=f'B{i:03d}'
 done=next((r for r in rows if r['index']==i),None)
 if done is not None:
  for module in done['modules']:
   assert sha(R/module['source'])==module['source_sha256']
  continue
 while True:
  try:progress=json.loads((D/'dense_fourier_witness_replay.json').read_text())
  except (FileNotFoundError,json.JSONDecodeError):time.sleep(5);continue
  if progress['checked']>wi:break
  time.sleep(5)
 gen=subprocess.run([sys.executable,str(R/'scripts/route_fourier_box.py'),str(i)],cwd=R,capture_output=True,text=True,encoding='utf-8')
 if gen.returncode:print(gen.stdout+gen.stderr,flush=True);sys.exit(gen.returncode)
 for suffix in ['Data','']:
  stem=f'DenseFourierExact{tag}{suffix}'
  p=subprocess.run([sys.executable,str(R/'scripts/peach-lean.py'),f'SpinCodes/Structured/{stem}.lean'],cwd=R,capture_output=True,text=True,encoding='utf-8')
  (D/f'dense_fourier_batch_{tag}{suffix}.log').write_text(p.stdout+p.stderr,encoding='utf-8')
  if p.returncode:print('FAIL',tag,suffix,p.stdout+p.stderr,flush=True);sys.exit(p.returncode)
 audits=re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]",(D/f'peach_DenseFourierExact{tag}.log').read_text(encoding='utf-8'))
 assert len(audits)==4
 for t,a in audits:assert set(x.strip() for x in a.split(','))<={'propext','Classical.choice','Quot.sound'},t
 record={'index':i,'global_index':gi,'witness_index':wi,'audits':4,'modules':[]}
 for suffix in ['Data','']:
  stem=f'DenseFourierExact{tag}{suffix}';s=R/f'SpinCodes/Structured/{stem}.lean';o=R/f'.lake/build/lib/lean/SpinCodes/Structured/{stem}.olean'
  record['modules'].append({'source':str(s.relative_to(R)),'source_sha256':sha(s),'object_sha256':sha(o)})
 rows.append(record)
 (D/config['manifest']).write_text(json.dumps({'status':'PASS' if len(rows)==len(config['indices']) else 'RUNNING','checked':len(rows),'total':len(config['indices']),'owned_indices':config['indices'],'common_prefactor':24000000000000,'scope':'Actual Fourier CertifiedBox proofs; four endpoint exponent certificates, exact matrix Collatz witness, actual outer support and original global box identification. Checks reuse checked dependencies.','boxes':rows},indent=2)+'\n')
 print('PASS',tag,f'({len(rows)}/{len(config["indices"])} {args.lane})',flush=True)
