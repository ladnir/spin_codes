"""Run a named generated LP once; retain logs; no concurrent child solvers."""
import argparse
import subprocess
import time
from pathlib import Path
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('folder')
parser.add_argument('--seconds',type=int,default=90)
parser.add_argument('--precision',type=int,default=256)
args=parser.parse_args()
assert args.folder.isidentifier() and 1<=args.seconds<=600
folder=ROOT/'generated'/args.folder
assert folder.is_dir() and (folder/'h_38.lp').is_file()
assert not (folder/'attempt.json').exists() and not (folder/'h_38.sol').exists()
command=['wsl.exe','-d','Ubuntu-24.04','--cd',
    '/mnt/c/Users/peter/.codex/worktrees/ba80/permute_conv/bch_spectrum_work/bch_spectrum_codex_bundle/generated/'+args.folder,
    '--','env','LD_LIBRARY_PATH=../oa21_closure_probe/tools/runtime/usr/lib/x86_64-linux-gnu',
    '../oa21_closure_probe/tools/runtime/usr/bin/esolver','-S','-d','7','-P',str(args.precision),'-L','-R',str(args.seconds),
    '-b','h_38.bas','-O','h_38.sol','h_38.lp']
start=time.monotonic()
result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
elapsed=time.monotonic()-start
with (folder/'esolver.log').open('xb') as f:
    f.write(result.stdout)
write_new(folder/'attempt.json',dict(command=command,elapsed_seconds=elapsed,exit_code=result.returncode,
    solution_exists=(folder/'h_38.sol').exists(),source_sha256=sha(Path(__file__))))
print(result.stdout.decode(errors='replace')[-4000:])
print('Elapsed seconds:',elapsed,'Exit:',result.returncode)
