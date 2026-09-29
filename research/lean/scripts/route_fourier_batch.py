"""Sequential proof-check replay for 105 Fourier witnesses; no benchmarks."""
from pathlib import Path
import subprocess,sys,json,hashlib,re,time
R=Path(__file__).resolve().parents[1];D=R/'scripts/map_data'
sys.stdout.reconfigure(encoding='utf-8')
rows=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for i in range(105):
 tag=f'W{i:03d}'
 gen=subprocess.run([sys.executable,str(R/'scripts/route_fourier_exact.py'),str(i)],cwd=R,capture_output=True,text=True,encoding='utf-8')
 if gen.returncode: print(gen.stdout+gen.stderr,flush=True);sys.exit(gen.returncode)
 for suffix in ['Data','']:
  stem=f'DenseFourierExact{tag}{suffix}'
  p=subprocess.run([sys.executable,str(R/'scripts/peach-lean.py'),f'SpinCodes/Structured/{stem}.lean'],cwd=R,capture_output=True,text=True,encoding='utf-8')
  (D/f'dense_fourier_batch_{tag}{suffix}.log').write_text(p.stdout+p.stderr,encoding='utf-8')
  if p.returncode:print('FAIL',tag,suffix,p.stdout+p.stderr,flush=True);sys.exit(p.returncode)
 audits=re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]",(D/f'peach_DenseFourierExact{tag}.log').read_text(encoding='utf-8'))
 assert len(audits)==3
 for t,a in audits:assert set(x.strip() for x in a.split(','))<={'propext','Classical.choice','Quot.sound'},t
 record={'index':i,'audits':3,'modules':[]}
 for suffix in ['Data','']:
  stem=f'DenseFourierExact{tag}{suffix}';s=R/f'SpinCodes/Structured/{stem}.lean';o=R/f'.lake/build/lib/lean/SpinCodes/Structured/{stem}.olean'
  record['modules'].append({'source':str(s.relative_to(R)),'source_sha256':sha(s),'object_sha256':sha(o)})
 rows.append(record)
 (D/'dense_fourier_witness_replay.json').write_text(json.dumps({'status':'PASS' if len(rows)==105 else 'RUNNING','checked':len(rows),'total':105,'scope':'Exact integer Collatz witnesses; local boxes/global coverage tracked separately. Individual checks reuse checked dependencies.','witnesses':rows},indent=2)+'\n')
 print('PASS',tag,f'({len(rows)}/105)',flush=True)
