from pathlib import Path
import runpy,subprocess,json,time
r=Path('.');ssh=runpy.run_path(str(r/'scripts/peach-lean.py'))['SSH'];rows=[]
for i in range(24):
 p=subprocess.run(ssh+["awk '/MemAvailable:/{print $2}' /proc/meminfo; ps -C lean -o rss=,args="],text=True,capture_output=True,check=True)
 lines=p.stdout.splitlines();row={'time':time.time(),'available_kib':int(lines[0]),'lean':lines[1:]};rows.append(row)
 out=r/'scripts/map_data/dense_occupation_admission_memory.json';out.write_text(json.dumps({'status':'OBSERVING' if i<23 else 'PASS','minimum_available_kib':min(x['available_kib'] for x in rows),'samples':rows},indent=2)+'\n')
 if row['available_kib']<20*1024**2:print('MEMORY FLOOR BREACH',row['available_kib'],flush=True);break
 if i%4==0:print('Memory available GiB',round(row['available_kib']/1024**2,1),flush=True)
 time.sleep(10)
print('minimum available GiB',min(x['available_kib'] for x in rows)/1024**2,flush=True)
