"""Two checked candidate leaves for a single unresolved Fourier case."""
import argparse
import json
from pathlib import Path
from bch_quotient import cyclotomic_coset
from probe_bch_shift_rank import search,verify_steps
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('source_folder')
    parser.add_argument('pivot',type=int)
    parser.add_argument('split',type=int)
    parser.add_argument('target',type=int)
    parser.add_argument('folder_name')
    parser.add_argument('--width',type=int,default=500)
    parser.add_argument('--seconds',type=int,default=30)
    args=parser.parse_args()
    assert args.folder_name.isidentifier() and args.source_folder.isidentifier()
    source=ROOT/'generated'/args.source_folder/f'case_{args.pivot:03d}.json'
    old=json.loads(source.read_text())
    zeros=set(old['zero_indices']); orbit=set(cyclotomic_coset(args.split))
    assert orbit.isdisjoint(zeros) and args.pivot not in orbit
    folder=ROOT/'generated'/args.folder_name
    folder.mkdir(exist_ok=False)
    results=[]
    for zero_branch in (False,True):
        childzeros=zeros|orbit if zero_branch else zeros
        pivot=args.pivot if zero_branch else args.split
        value=search(childzeros,pivot,args.target,args.width,args.seconds)
        verify_steps(childzeros,value)
        value.update(zero_indices=sorted(childzeros),split_coset_is_zero=zero_branch)
        write_new(folder/('zero.json' if zero_branch else 'nonzero.json'),value)
        brief={k:value[k] for k in ('pivot','rank','target_reached','stop_reason','elapsed_seconds','split_coset_is_zero')}
        results.append(brief)
        print(json.dumps(brief),flush=True)
    write_new(folder/'summary.json',dict(classification='Candidate binary refinement of one Fourier case; full cover audit required',
        source_folder=args.source_folder,source_case_pivot=args.pivot,split_coset_representative=args.split,
        target_rank=args.target,results=results,all_targets_reached=all(x['target_reached'] for x in results),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/probe_bch_shift_rank.py',ROOT/'code/bch_quotient.py',source)}))
