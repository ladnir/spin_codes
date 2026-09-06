"""Bounded search for explicit independent-set rank witnesses; no implied exhaustiveness."""
import argparse
import json
import time
from pathlib import Path
from bch_quotient import defining_cosets,cyclotomic_coset
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
MASK=(1<<255)-1


def rotate(word,shift):
    shift%=255
    return ((word<<shift)|(word>>(255-shift)))&MASK


def indices(word):
    while word:
        bit=word&-word
        yield bit.bit_length()-1
        word^=bit


def search(zeros,pivot,target,width,seconds):
    start=time.monotonic()
    # Work in relative exponents: the known nonzero Fourier index is 0.
    relative={(z-pivot)%255 for z in zeros}
    assert 0 not in relative
    zword=sum(1<<i for i in relative)
    allowed=[rotate(zword,-i) for i in range(255)]
    powers=[[i*(1<<k)%255 for i in range(255)] for k in range(8)]
    def choices(values):
        out=MASK
        for i in values:
            out&=allowed[i]
            if not out:
                break
        return out
    beam={1:[]}
    best=(1,[])
    levels=[]
    stopped='no further state'
    for depth in range(2,target+1):
        candidates={}
        for state,path in beam.items():
            support=list(indices(state))
            for k in range(8):
                moved=[powers[k][i] for i in support]
                word=sum(1<<i for i in moved)
                for shift in indices(choices(moved)):
                    new=rotate(word,shift)|1
                    if new not in candidates:
                        candidates[new]=path+[(k,shift)]
            if time.monotonic()-start>seconds:
                stopped='time limit'
                break
        if not candidates:
            break
        levels.append(dict(rank=depth,candidates=len(candidates)))
        if depth==target:
            best=min(candidates.items())
            stopped='target reached'
            break
        # Prefer states with more legal next shifts; hash tie-breaking retains diversity.
        ranked=[]
        for word,path in candidates.items():
            support=list(indices(word))
            score=sum(choices([powers[k][i] for i in support]).bit_count() for k in range(8))
            if score:
                tie=(word^(word>>64)^(word>>128)^(word>>192))*0x9e3779b97f4a7c15&((1<<64)-1)
                ranked.append((score,tie,word,path))
        best=min(candidates.items())
        if not ranked or stopped=='time limit':
            break
        ranked.sort(reverse=True)
        beam={w:p for _,_,w,p in ranked[:width]}
    word,path=best
    return dict(pivot=pivot,rank=len(path)+1,relative_independent_set=list(indices(word)),
        steps=path,levels=levels,stop_reason=stopped,elapsed_seconds=time.monotonic()-start,
        target_reached=len(path)+1>=target)


def verify_steps(zeros,record):
    pivot=record['pivot']; state={pivot}
    assert pivot not in zeros
    for k,shift in record['steps']:
        # Relative normalization corresponds to this actual exponent shift.
        actual_shift=(shift+(1-(1<<k))*pivot)%255
        moved={((1<<k)*i+actual_shift)%255 for i in state}
        assert len(moved)==len(state) and moved<=zeros
        state=moved|{pivot}
    assert len(state)==record['rank']
    assert sorted((i-pivot)%255 for i in state)==record['relative_independent_set']


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('distance',type=int,choices=[37,39])
    parser.add_argument('pivot',type=int)
    parser.add_argument('folder_name')
    parser.add_argument('--width',type=int,default=1000)
    parser.add_argument('--seconds',type=int,default=30)
    args=parser.parse_args()
    assert args.folder_name.isidentifier()
    original=set().union(*defining_cosets(args.distance))
    zeros={-i%255 for i in set(range(255))-original}
    target=25 if args.distance==37 else 23
    value=search(zeros,args.pivot,target,args.width,args.seconds)
    verify_steps(zeros,value)
    value.update(classification='Verified independent-set witness for ONE nonzero Fourier case; NOT a full distance certificate',
        original_BCH_designed_distance=args.distance,dual_zero_indices=sorted(zeros),
        full_dual_distance_certified=False,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/bch_quotient.py')})
    folder=ROOT/'generated'/args.folder_name
    folder.mkdir(exist_ok=False)
    write_new(folder/'receipt.json',value)
    print(json.dumps({k:v for k,v in value.items() if k not in ('source_sha256','dual_zero_indices','steps','relative_independent_set')},indent=2))
