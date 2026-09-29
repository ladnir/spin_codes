"""Read only: sample ordinary occupation semantic proof-check duration and RSS."""
from pathlib import Path
import runpy,subprocess,json
r=Path(__file__).resolve().parents[1];ssh=runpy.run_path(str(r/'scripts/peach-lean.py'))['SSH']
remote=r'''import pathlib,time,json,re,os
p=pathlib.Path('/proc');ticks=os.sysconf('SC_CLK_TCK');rows={};started=time.monotonic()
while time.monotonic()-started<500:
 current=set()
 for q in p.iterdir():
  if not q.name.isdigit():continue
  try:
   a=(q/'cmdline').read_bytes().split(b'\0');name=pathlib.Path(a[0].decode()).name
   if name!='lean':continue
   matches=[s.decode() for s in a if re.search(rb'DenseOccupationFixedB[0-9]{3}(Rate)?\.lean$',s)]
   if not matches:continue
   stat=(q/'stat').read_text().split();status=(q/'status').read_text()
  except (FileNotFoundError,PermissionError,ProcessLookupError,IndexError):continue
  pid=int(q.name);current.add(pid)
  if pid not in rows:rows[pid]={'module':matches[0],'start_boot_seconds':int(stat[21])/ticks,'peak_rss_kib':0,'finished':False}
  vals={line.split(':')[0]:int(line.split()[1]) for line in status.splitlines() if line.startswith(('VmHWM:','VmRSS:'))}
  rows[pid]['peak_rss_kib']=max(rows[pid]['peak_rss_kib'],vals.get('VmHWM',0))
 for pid,row in rows.items():
  if not row['finished'] and pid not in current:
   row.update(finished=True,elapsed_seconds_upper=time.monotonic()-row['start_boot_seconds'])
 if sum(v['finished'] for v in rows.values())>=12:break
 time.sleep(.5)
print(json.dumps({'status':'PASS','measurements':rows,'poll_resolution_seconds':.5,'note':'Passive observation of actual kernel proof checks; no benchmark.'}))
'''
p=subprocess.run(ssh+['python3 -'],input=remote,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=550)
print(p.stdout,flush=True)
if p.returncode:raise SystemExit(p.returncode)
d=json.loads(p.stdout.strip().splitlines()[-1]);(r/'scripts/map_data/dense_occupation_semantic_highwater.json').write_text(json.dumps(d,indent=2)+'\n')
