from pathlib import Path
import subprocess,fcntl,os
r=Path('/tmp/spin-forward-20261001');locks=[]
for name in ['prindal-addition-encoder-benchmark','bare-spin-benchmark','hypercat-benchmark','hypercat-global-benchmark']:
 f=open('/tmp/'+name+'.lock','a');fcntl.flock(f,fcntl.LOCK_EX);locks.append(f)
b=r/'sanitized'
with (r/'sanitizer.log').open('w') as log:
 subprocess.run(['cmake','-S',str(r/'spin'),'-B',str(b),'-DCMAKE_BUILD_TYPE=Debug','-DCMAKE_CXX_FLAGS=-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer','-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined'],stdout=log,stderr=subprocess.STDOUT,check=True)
 subprocess.run(['cmake','--build',str(b),'--target','spin_packet_test','-j4'],stdout=log,stderr=subprocess.STDOUT,check=True)
 subprocess.run([str(b/'spin_packet_test')],env=dict(os.environ,ASAN_OPTIONS='detect_leaks=1',UBSAN_OPTIONS='halt_on_error=1'),stdout=log,stderr=subprocess.STDOUT,check=True)
print('Full-library packet ASan/UBSan checks passed.',flush=True)
with (r/'library-tests.log').open('w') as log:
 subprocess.run(['cmake','--build',str(r/'build'),'--target','spin_api_test','spin_packet_test','spin_known_answers','spin_wide_routing','spin_setup_test','spin_prepared_masks','-j','4'],stdout=log,stderr=subprocess.STDOUT,check=True)
 subprocess.run(['ctest','--test-dir',str(r/'build'),'-R','^(spin_packet|spin_api|spin_known_answers|spin_wide_routing|spin_setup|spin_prepared_masks)$','--output-on-failure','-j','1'],stdout=log,stderr=subprocess.STDOUT,check=True)
print('Six existing library tests passed serially.',flush=True)
