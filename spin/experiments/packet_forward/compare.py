"""Serial paired comparison against the unchanged public packet transpose.

Run on Peach after run_remote.py build. Retains every forward timing sample.
"""
from pathlib import Path
import fcntl,subprocess,os,json,time,csv,io,hashlib,statistics
r=Path('/tmp/spin-forward-20261001');locks=[]
for name in ['prindal-addition-encoder-benchmark','bare-spin-benchmark','hypercat-benchmark','hypercat-global-benchmark']:
 f=open('/tmp/'+name+'.lock','a');print('Waiting for '+name,flush=True);fcntl.flock(f,fcntl.LOCK_EX);locks.append(f)
freq=[Path(p).read_text().strip() for p in ['/sys/devices/system/cpu/cpu15/cpufreq/scaling_governor','/sys/devices/system/cpu/cpu15/cpufreq/scaling_setspeed','/sys/devices/system/cpu/cpufreq/boost']]
assert freq==['userspace','4500000','0']
out=r/('paired-'+time.strftime('%Y%m%d-%H%M%S'));out.mkdir()
rec={'base_revision':'36394342820da99291cc736165ff13129a2ca359','frequency':freq,'cpu':15,'calls':501,'warmups':5,'runs':[],
 'source_sha256':{str(p.relative_to(r/'spin')):hashlib.sha256(p.read_bytes()).hexdigest() for p in (r/'spin').rglob('*') if p.is_file()},
 'method':'Two seeds (1,17), both forward/transpose execution orders, median of four process medians. Precomputed encoding only; setup/allocation excluded. One process at a time under four shared locks. Normal allocation at K16/K18; PreferHugePages at K20.'}
expected={}
for log in [16,18,20]:
 memory='huge' if log==20 else 'normal'
 for seed in [1,17]:
  for order in [0,1]:
   for direction in (['transpose','forward'] if order==0 else ['forward','transpose']):
    tag=f'{log}-{seed}-{order}-{direction}'
    binary=r/'build'/('spin_packet_bench' if direction=='transpose' else 'spin_packet_forward_bench')
    args=['taskset','-c','15',str(binary)]+([] if direction=='transpose' else ['selected'])+[str(2**log),str(seed),'501',memory]
    env=dict(os.environ,SPIN_FORWARD_RAW=str(out/(tag+'-trials.csv')))
    start=time.monotonic()
    with (out/(tag+'.stdout')).open('w') as stdout,(out/(tag+'.stderr')).open('w') as stderr:
     child=subprocess.Popen(args,env=env,stdout=stdout,stderr=stderr)
     _,status,usage=os.wait4(child.pid,0);child.returncode=os.waitstatus_to_exitcode(status)
     assert child.returncode==0,(tag,child.returncode)
    row=list(csv.DictReader(io.StringIO((out/(tag+'.stdout')).read_text())))[0]
    key=(log,seed,direction);assert row['checksum']==expected.setdefault(key,row['checksum'])
    record=dict(tag=tag,direction=direction,K=2**log,seed=seed,order=order,command=args,exit_code=child.returncode,
       elapsed_seconds=time.monotonic()-start,max_rss_kib=usage.ru_maxrss,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),result=row)
    rec['runs'].append(record);(out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(tag,row['median_ms'],flush=True)
 rec.setdefault('summary',{})[str(log)]={d:statistics.median(float(x['result']['median_ms']) for x in rec['runs'] if x['K']==2**log and x['direction']==d) for d in ['forward','transpose']}
(out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec['summary'],indent=2),flush=True);print(out,flush=True)
