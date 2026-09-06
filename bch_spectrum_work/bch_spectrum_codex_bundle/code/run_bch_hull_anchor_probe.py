"""Bounded sequential solver run for the derived hull anchor."""
import subprocess
import time
from pathlib import Path
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_hull_anchor_probe'
assert not (FOLDER/'attempt.json').exists() and not (FOLDER/'h_38.sol').exists()
command=['wsl.exe','-d','Ubuntu-24.04','--cd',
    '/mnt/c/Users/peter/.codex/worktrees/ba80/permute_conv/bch_spectrum_work/bch_spectrum_codex_bundle/generated/oa21_hull_anchor_probe',
    '--','env','LD_LIBRARY_PATH=../oa21_closure_probe/tools/runtime/usr/lib/x86_64-linux-gnu',
    '../oa21_closure_probe/tools/runtime/usr/bin/esolver','-S','-d','7','-P','256','-L','-R','90',
    '-b','h_38.bas','-O','h_38.sol','h_38.lp']
start=time.monotonic()
result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
elapsed=time.monotonic()-start
with (FOLDER/'esolver.log').open('xb') as f:
    f.write(result.stdout)
write_new(FOLDER/'attempt.json',dict(command=command,elapsed_seconds=elapsed,exit_code=result.returncode,
    solution_exists=(FOLDER/'h_38.sol').exists(),source_sha256=sha(Path(__file__))))
print(result.stdout.decode(errors='replace')[-4000:])
print('Elapsed seconds:',elapsed,'Exit:',result.returncode)
