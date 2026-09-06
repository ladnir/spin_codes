"""Sequential independent-set searches covering first-nonzero Fourier cosets."""
import argparse
import json
from pathlib import Path
from bch_quotient import defining_cosets,cyclotomic_coset
from probe_bch_shift_rank import search,verify_steps
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('distance',type=int,choices=[37,39])
    parser.add_argument('folder_name')
    parser.add_argument('--width',type=int,default=300)
    parser.add_argument('--seconds',type=int,default=20)
    args=parser.parse_args()
    assert args.folder_name.isidentifier()
    folder=ROOT/'generated'/args.folder_name
    folder.mkdir(exist_ok=False)
    original=set().union(*defining_cosets(args.distance))
    base={-i%255 for i in set(range(255))-original}
    remaining=set(range(255))-base
    reps=sorted({min(cyclotomic_coset(i)) for i in remaining})
    zeros=set(base)
    target=25 if args.distance==37 else 23
    records=[]
    for pivot in reps:
        orbit=set(cyclotomic_coset(pivot))
        assert orbit<=remaining and not orbit&zeros
        value=search(zeros,pivot,target,args.width,args.seconds)
        verify_steps(zeros,value)
        value['zero_indices']=sorted(zeros)
        write_new(folder/f'case_{pivot:03d}.json',value)
        records.append({k:value[k] for k in ('pivot','rank','target_reached','elapsed_seconds','stop_reason')})
        print(json.dumps(records[-1]),flush=True)
        zeros|=orbit
        remaining-=orbit
    assert not remaining and zeros==set(range(255))
    summary=dict(classification='Candidate complete case cover; requires independent mathematical audit',
        original_BCH_designed_distance=args.distance,target_rank=target,case_representatives=reps,cases=records,
        all_targets_reached=all(r['target_reached'] for r in records),
        independently_certified=False,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/probe_bch_shift_rank.py',ROOT/'code/bch_quotient.py')})
    write_new(folder/'summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('source_sha256','cases')},indent=2),flush=True)
