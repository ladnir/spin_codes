from pathlib import Path
import fcntl,subprocess,os,sys,json,time
r=Path('/tmp/spin-forward-20261001');locks=[]
for name in ['prindal-addition-encoder-benchmark','bare-spin-benchmark','hypercat-benchmark','hypercat-global-benchmark']:
 f=open('/tmp/'+name+'.lock','a');print('Waiting for '+name,flush=True);fcntl.flock(f,fcntl.LOCK_EX);locks.append(f)
env=dict(os.environ);env.pop('CXXFLAGS',None)
build=r/'build';log=r/'build.log'
with log.open('w') as out:
 subprocess.run(['cmake','-S',str(r/'spin'),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DSPIN_BUILD_EXPERIMENTS=ON','-DSPIN_BUILD_BENCHMARKS=ON','-DSPIN_TUNE=znver4'],env=env,stdout=out,stderr=subprocess.STDOUT,check=True)
 subprocess.run(['cmake','--build',str(build),'--target','spin_packet_forward_bench','spin_packet_bench','-j','4'],env=env,stdout=out,stderr=subprocess.STDOUT,check=True)
subprocess.run([str(build/'spin_packet_forward_bench'),'verify'],check=True)
if len(sys.argv)>1 and sys.argv[1]=='build':sys.exit(0)
for logk in [16,18,20]:
 for mode in ['transpose','forward','gather']:
  args=['taskset','-c','15',str(build/'spin_packet_forward_bench'),mode,str(2**logk),'1','51','huge' if logk==20 else 'normal']
  print(subprocess.check_output(args,text=True),flush=True)

