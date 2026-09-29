"""Read only: sample peak RSS of one existing occupation data-check wave."""
from pathlib import Path
import runpy,subprocess,json
r=Path(__file__).resolve().parents[1]
ssh=runpy.run_path(str(r/'scripts/peach-lean.py'))['SSH']
remote=r'''import pathlib,time,json
p=pathlib.Path('/proc');targets={f'DenseOccupationFixedW{i:03d}Data.lean' for i in range(63,67)}
rows={}
for q in p.iterdir():
 if not q.name.isdigit():continue
 try:a=(q/'cmdline').read_bytes().split(b'\0');name=pathlib.Path(a[0].decode()).name
 except (FileNotFoundError,PermissionError,IndexError):continue
 if name!='lean':continue
 matches=[s.decode() for s in a if any(s.endswith(t.encode()) for t in targets)]
 if matches:rows[int(q.name)]={'module':matches[0],'peak_rss_kib':0,'last_rss_kib':0,'finished':False}
print('MONITOR',list(rows),flush=True)
for _ in range(80):
 for pid,row in rows.items():
  if row['finished']:continue
  try:status=(p/str(pid)/'status').read_text()
  except FileNotFoundError:row['finished']=True;continue
  vals={line.split(':')[0]:int(line.split()[1]) for line in status.splitlines() if line.startswith(('VmHWM:','VmRSS:'))}
  row['peak_rss_kib']=max(row['peak_rss_kib'],vals.get('VmHWM',0));row['last_rss_kib']=vals.get('VmRSS',0)
 if all(v['finished'] for v in rows.values()):break
 time.sleep(10)
print(json.dumps({'status':'PASS','measurements':rows,'note':'Observed kernel proof-check processes only; no benchmark was run.'}),flush=True)
'''
p=subprocess.run(ssh+['python3 -'],input=remote,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=850)
print(p.stdout,flush=True)
if p.returncode:raise SystemExit(p.returncode)
d=json.loads(p.stdout.strip().splitlines()[-1]);(r/'scripts/map_data/dense_occupation_memory_highwater.json').write_text(json.dumps(d,indent=2)+'\n')
