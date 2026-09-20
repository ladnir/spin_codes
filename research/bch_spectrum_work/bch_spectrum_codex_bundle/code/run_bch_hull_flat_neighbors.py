"""Compile and run a prepared affine-flat neighbor search, once."""
import argparse
import json
import subprocess
import time
from pathlib import Path
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('folder_name')
    args=parser.parse_args()
    assert args.folder_name.isidentifier()
    folder=ROOT/'generated'/args.folder_name
    assert (folder/'input.bin').is_file() and not (folder/'attempt.json').exists() and not (folder/'seeds.bin').exists()
    wroot='/mnt/c/'+str(ROOT.resolve())[3:].replace('\\','/')
    stem='generated/'+args.folder_name+'/'
    compile_command=['wsl.exe','-d','Ubuntu-24.04','--cd',wroot,'--','g++','-O3','-march=native','-std=c++20',
        'code/bch_hull_flat_neighbors.cpp','-o',stem+'search']
    subprocess.run(compile_command,check=True)
    command=['wsl.exe','-d','Ubuntu-24.04','--cd',wroot,'--',stem+'search',stem+'input.bin',stem+'seeds.bin']
    start=time.monotonic()
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True,timeout=180)
    elapsed=time.monotonic()-start
    stats=json.loads(result.stdout)
    write_new(folder/'attempt.json',dict(compile_command=compile_command,command=command,elapsed_seconds=elapsed,
        exit_code=result.returncode,statistics=stats,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/bch_hull_flat_neighbors.cpp',folder/'input.bin',folder/'search',folder/'seeds.bin')}))
    print(json.dumps(dict(elapsed_seconds=elapsed,**stats),indent=2))
