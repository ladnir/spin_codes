"""Bounded discovery of P-dual rank witnesses under explicit Fourier splits.

Each output proves only its saved conditional case. A separate tree audit
must verify the complete partition before any global distance improvement.
"""
import argparse
from pathlib import Path
import sys
import bridge as base

sys.path.insert(0,str(base.BCH/'code'))
from certify_bch_shift_rank import code_data, orbit, check_path
from probe_bch_shift_rank import search
from complete_bch_shift_rank_cases import best_run


def assumptions(parent, zero, nonzero):
    data=code_data(37)
    zeros=set(data['dual_zero_indices'])
    reps=sorted({min(orbit(e)) for e in set(range(255))-zeros})
    assert parent in reps
    for r in reps:
        if r==parent:break
        zeros |= orbit(r)
    known=orbit(parent)
    for r in zero:
        assert r==min(orbit(r)) and orbit(r).isdisjoint(known)
        zeros |= orbit(r)
    for r in nonzero:
        assert r==min(orbit(r)) and orbit(r).isdisjoint(zeros)
        known |= orbit(r)
    assert zeros.isdisjoint(known)
    return zeros,known


def verify(record):
    zeros,known=assumptions(record['parent'],record['zero_cosets'],record['nonzero_cosets'])
    assert record['zero_indices']==sorted(zeros) and record['known_nonzero_indices']==sorted(known)
    witness=record['witness']
    if witness.get('proof_type')=='BCH_root_run':
        import math
        start,step,count=(witness[k] for k in ('start','step','root_count'))
        assert 0<=start<255 and 0<step<255 and math.gcd(step,255)==1
        assert 0<count<255 and all((start+j*step)%255 in zeros for j in range(count))
        rank=count+1
    else:
        assert witness['pivot'] in known
        state=check_path(zeros,witness['pivot'],witness['steps'])
        assert sorted((e-witness['pivot'])%255 for e in state)==witness['relative_independent_set']
        rank=len(state)
    assert rank==witness['rank']
    return rank


def run(args):
    if args.verify:
        record=base.read(args.output)
        for name,digest in record['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
        print('Independently verified conditional rank',verify(record),flush=True)
        return
    assert not args.output.exists()
    zeros,known=assumptions(args.parent,args.zero,args.nonzero)
    assert args.pivot in known
    count,start,step=best_run(zeros)
    if count+1>=31:
        witness=dict(proof_type='BCH_root_run',start=start,step=step,root_count=count,rank=count+1)
    else:
        witness=search(zeros,args.pivot,31,args.width,args.seconds)
    record=dict(status='CONDITIONAL_CASE_WITNESS_NOT_A_COMPLETE_DISTANCE_PROOF',parent=args.parent,
        zero_cosets=args.zero,nonzero_cosets=args.nonzero,zero_indices=sorted(zeros),
        known_nonzero_indices=sorted(known),witness=witness,required_rank=31,
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in
            (Path(__file__),base.BCH/'code/certify_bch_shift_rank.py',base.BCH/'code/probe_bch_shift_rank.py',
             base.BCH/'code/complete_bch_shift_rank_cases.py',base.BCH/'code/bch_quotient.py',
             base.BCH/'code/affine_wambach.py')})
    rank=verify(record)
    base.write_new(args.output,record)
    print(dict(parent=args.parent,zero=args.zero,nonzero=args.nonzero,pivot=args.pivot,rank=rank,
               stop=witness.get('stop_reason','root run'),elapsed=witness.get('elapsed_seconds')),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parent',type=int,choices=[7,15],default=7)
    p.add_argument('--zero',type=int,nargs='*',default=[])
    p.add_argument('--nonzero',type=int,nargs='*',default=[])
    p.add_argument('--pivot',type=int,default=7)
    p.add_argument('--width',type=int,default=600)
    p.add_argument('--seconds',type=int,default=20)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
