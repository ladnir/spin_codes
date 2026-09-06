"""Compile and run the fixed-width orbit kernel once in a new folder."""
import argparse
import json
import subprocess
import time
from pathlib import Path
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('weight',type=int,choices=[14,16,18,20,22])
    parser.add_argument('seed_folder')
    parser.add_argument('output_folder')
    args=parser.parse_args()
    assert args.seed_folder.isidentifier() and args.output_folder.isidentifier()
    folder=ROOT/'generated'/args.output_folder
    source=ROOT/'code/bch_affine_canonical.cpp'
    seed=ROOT/'generated'/args.seed_folder/'seeds.bin'
    assert seed.is_file()
    folder.mkdir(exist_ok=False)
    wroot='/mnt/c/'+str(ROOT.resolve())[3:].replace('\\','/')
    binary=str((folder/'canonical').relative_to(ROOT)).replace('\\','/')
    compile_command=['wsl.exe','-d','Ubuntu-24.04','--cd',wroot,'--','g++','-O3','-march=native','-std=c++20',
        'code/bch_affine_canonical.cpp','-o',binary]
    subprocess.run(compile_command,check=True)
    command=['wsl.exe','-d','Ubuntu-24.04','--cd',wroot,'--',binary,str(args.weight),
        str(seed.relative_to(ROOT)).replace('\\','/'),binary+'.bin']
    start=time.monotonic()
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True,timeout=180)
    stats=json.loads(result.stdout)
    value=dict(weight=args.weight,seed_folder=args.seed_folder,compile_command=compile_command,command=command,
        elapsed_seconds=time.monotonic()-start,statistics=stats,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),source,seed,folder/'canonical',folder/'canonical.bin')})
    write_new(folder/'attempt.json',value)
    print(json.dumps({k:v for k,v in value.items() if k not in ('source_sha256','compile_command','command')},indent=2))
