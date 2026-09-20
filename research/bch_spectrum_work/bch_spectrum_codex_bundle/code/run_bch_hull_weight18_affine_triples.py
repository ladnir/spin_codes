"""Compile and run the expanded candidate search sequentially, once."""
import json
import subprocess
import time
from pathlib import Path
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/hull_weight18_affine_triples'
FOLDER.mkdir(exist_ok=False)
wroot='/mnt/c/'+str(ROOT.resolve())[3:].replace('\\','/')
compile_command=['wsl.exe','-d','Ubuntu-24.04','--cd',wroot,'--','g++','-O3','-march=native','-std=c++20',
    'code/bch_hull_weight18_affine_triples.cpp','-o','generated/hull_weight18_affine_triples/search']
subprocess.run(compile_command,check=True)
command=['wsl.exe','-d','Ubuntu-24.04','--cd',wroot,'--','generated/hull_weight18_affine_triples/search',
    'generated/hull_weight18_triples/input.bin','generated/hull_weight18_affine_triples/seeds.bin']
start=time.monotonic()
result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
elapsed=time.monotonic()-start
stats=json.loads(result.stdout)
write_new(FOLDER/'attempt.json',dict(compile_command=compile_command,command=command,elapsed_seconds=elapsed,
    exit_code=result.returncode,statistics=stats,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
        (Path(__file__),ROOT/'code/bch_hull_weight18_affine_triples.cpp',
         ROOT/'generated/hull_weight18_triples/input.bin',FOLDER/'search',FOLDER/'seeds.bin')}))
print(json.dumps(dict(elapsed_seconds=elapsed,**stats),indent=2))
