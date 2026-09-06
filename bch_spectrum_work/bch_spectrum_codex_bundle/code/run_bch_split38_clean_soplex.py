"""Bounded multiprecision/exact attempt on the cleaned split relaxation."""
import subprocess
from pathlib import Path
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/split38_clean_probe'
OUTPUT=FOLDER/'soplex_mp_exact.log'
linux='/mnt/c/Users/peter/.codex/worktrees/ba80/permute_conv/bch_spectrum_work/bch_spectrum_codex_bundle/generated/split38_clean_probe'
command=['wsl.exe','-d','Ubuntu-24.04','--cd',linux,'--','env',
         'LD_LIBRARY_PATH=../split38_soplex_tools/runtime/usr/lib/x86_64-linux-gnu',
         '../split38_soplex_tools/runtime/usr/bin/soplex','--arithmetic=2','--precision=100',
         '--readmode=1','--solvemode=2','--int:syncmode=1','--int:checkmode=2',
         '-f0','-o0','-t120','-s0','-X','-Y','-q','h38.lp']
print('Starting cleaned exact LP with 100-digit arithmetic and 120-second limit',flush=True)
with OUTPUT.open('x',encoding='utf-8') as stream:
    result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
write_new(FOLDER/'soplex_mp_attempt.json',dict(classification='Process receipt, not a proof certificate',
          command=command,exit_code=result.returncode,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
          (Path(__file__),FOLDER/'model.json',FOLDER/'h38.lp',OUTPUT)}))
print('Process exit',result.returncode,flush=True)
print('\n'.join(OUTPUT.read_text().splitlines()[-18:]))
