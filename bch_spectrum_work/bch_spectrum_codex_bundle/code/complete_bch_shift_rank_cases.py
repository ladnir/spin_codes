"""Complete a Fourier-case cover, preferring explicit root runs and retained witnesses."""
import argparse
import json
import math
from pathlib import Path
from bch_quotient import defining_cosets,cyclotomic_coset
from probe_bch_shift_rank import search,verify_steps
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]


def best_run(zeros):
    best=(0,0,1)
    for step in range(1,255):
        if math.gcd(step,255)!=1:
            continue
        for start in zeros:
            if (start-step)%255 in zeros:
                continue
            count=0
            while (start+count*step)%255 in zeros:
                count+=1
            best=max(best,(count,start,step))
    return best


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('distance',type=int,choices=[37,39])
    parser.add_argument('source_folder')
    parser.add_argument('folder_name')
    parser.add_argument('target',type=int)
    parser.add_argument('--width',type=int,default=300)
    parser.add_argument('--seconds',type=int,default=25)
    args=parser.parse_args()
    assert args.folder_name.isidentifier() and (args.source_folder=='-' or args.source_folder.isidentifier())
    folder=ROOT/'generated'/args.folder_name
    folder.mkdir(exist_ok=False)
    original=set().union(*defining_cosets(args.distance))
    base={-i%255 for i in set(range(255))-original}
    reps=sorted({min(cyclotomic_coset(i)) for i in set(range(255))-base})
    zeros=set(base); records=[]; paths=[]
    for pivot in reps:
        count,start,step=best_run(zeros)
        old=ROOT/'generated'/args.source_folder/f'case_{pivot:03d}.json'
        if count+1>=args.target:
            value=dict(pivot=pivot,rank=count+1,proof_type='BCH_root_run',
                start=start,step=step,root_count=count,target_reached=True)
        elif args.source_folder!='-' and old.is_file():
            value=json.loads(old.read_text())
            assert value['zero_indices']==sorted(zeros)
            verify_steps(zeros,value)
            if value['rank']<args.target:
                value=search(zeros,pivot,args.target,args.width,args.seconds)
            else:
                paths.append(old)
        else:
            value=search(zeros,pivot,args.target,args.width,args.seconds)
        if value.get('proof_type')!='BCH_root_run':
            verify_steps(zeros,value)
        value['zero_indices']=sorted(zeros)
        value['target_reached']=value['rank']>=args.target
        write_new(folder/f'case_{pivot:03d}.json',value)
        records.append(dict(pivot=pivot,rank=value['rank'],proof_type=value.get('proof_type','independent_set'),
            target_reached=value['target_reached']))
        print(json.dumps(records[-1]),flush=True)
        zeros|=set(cyclotomic_coset(pivot))
    assert zeros==set(range(255))
    summary=dict(classification='Candidate full Fourier-case cover; independent audit required',
        original_BCH_designed_distance=args.distance,target_rank=args.target,case_representatives=reps,cases=records,
        all_targets_reached=all(r['target_reached'] for r in records),independently_certified=False,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            [Path(__file__),ROOT/'code/probe_bch_shift_rank.py',ROOT/'code/bch_quotient.py']+paths})
    write_new(folder/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('source_sha256','cases')},indent=2))
