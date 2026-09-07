"""Rank discovery using every specified nonzero Fourier coset as an append pivot.

Search is heuristic; the independent set-based verifier checks each transition.
No completeness of the beam is asserted or needed for a successful witness.
"""
import argparse
from pathlib import Path
import time
import bridge as base
import pdual_rank_refinement as single

MASK=(1<<255)-1


def positions(word):
    while word:
        bit=word&-word
        yield bit.bit_length()-1
        word ^= bit


def rotate(word,shift):
    return ((word<<shift)|(word>>(255-shift)))&MASK


def check(zeros,known,witness):
    assert witness['initial_pivot'] in known
    state={witness['initial_pivot']}
    for k,shift,pivot in witness['steps']:
        assert type(k) is int and 0<=k<8
        assert type(shift) is int and 0<=shift<255 and pivot in known
        moved={(e*(1<<k)+shift)%255 for e in state}
        assert len(moved)==len(state) and moved<=zeros and pivot not in moved
        state=moved|{pivot}
    assert witness['rank']==len(state) and witness['independent_set']==sorted(state)
    return len(state)


def discover(zeros,known,width,seconds,seed_records,start_rank,tie_seed):
    begun=time.monotonic();deadline=begun+seconds
    zword=sum(1<<e for e in zeros)
    allowed=[rotate(zword,(255-i)%255) for i in range(255)]
    powers=[[i*(1<<k)%255 for i in range(255)] for k in range(8)]
    pivots=sorted({min(single.orbit(e)) for e in known})
    def moves(state):
        indices=list(positions(state))
        result=[]
        for k,power in enumerate(powers):
            shifts=MASK;word=0
            for i in indices:
                moved=power[i];word |= 1<<moved;shifts &= allowed[moved]
                if not shifts:break
            if shifts:result.append((k,word,shifts))
        return result
    seeds={}
    for record in seed_records:
        witness=record['witness'] if 'witness' in record else record
        if 'steps' not in witness:continue
        if 'initial_pivot' in witness:
            initial=witness['initial_pivot'];steps=witness['steps']
        else:
            initial=witness['pivot']
            steps=[(k,(s+(1-(1<<k))*initial)%255,initial) for k,s in witness['steps']]
        if initial not in known:continue
        state={initial};path=[]
        seeds.setdefault(1,{})[1<<initial]=(initial,[])
        for k,s,pivot in steps:
            moved={(e*(1<<k)+s)%255 for e in state}
            if not (moved<=zeros and pivot in known and pivot not in moved):break
            state=moved|{pivot};path=path+[(k,s,pivot)]
            seeds.setdefault(len(state),{})[sum(1<<e for e in state)]=(initial,path)
    depth=min(start_rank,max(seeds,default=1))
    beam=seeds.get(depth,{1<<p:(p,[]) for p in pivots})
    best_word,(best_initial,best_path)=next(iter(beam.items()))
    levels=[];stop='no further state'
    for rank in range(depth+1,32):
        candidates=dict(seeds.get(rank,{}))
        for index,(state,(initial,path)) in enumerate(beam.items()):
            for k,word,shifts in moves(state):
                for shift in positions(shifts):
                    moved=rotate(word,shift)
                    for pivot in pivots:
                        new=moved|(1<<pivot)
                        if new not in candidates:candidates[new]=(initial,path+[(k,shift,pivot)])
            if time.monotonic()>=deadline:
                stop='time limit';break
        if not candidates:break
        best_word,(best_initial,best_path)=next(iter(candidates.items()))
        levels.append(dict(rank=rank,candidates=len(candidates)))
        if rank>=31:
            stop='target reached';break
        ranked=[]
        for index,(word,(initial,path)) in enumerate(candidates.items()):
            score=sum(shifts.bit_count() for _,_,shifts in moves(word))
            if score:
                tie=((word^(word>>64)^(word>>128)^(word>>192)^tie_seed)*0x9e3779b97f4a7c15)&((1<<64)-1)
                ranked.append((score,tie,word,initial,path))
            if index%32==0 and time.monotonic()>=deadline:
                stop='time limit';break
        if not ranked or stop=='time limit':break
        ranked.sort(reverse=True)
        beam={w:(initial,path) for _,_,w,initial,path in ranked[:width]}
        beam.update(seeds.get(rank,{}))
    witness=dict(initial_pivot=best_initial,steps=best_path,rank=len(best_path)+1,
                 independent_set=sorted(positions(best_word)))
    check(zeros,known,witness)
    return witness,dict(levels=levels,stop_reason=stop,elapsed_seconds=time.monotonic()-begun,
                        width=width,start_rank=depth,tie_seed=tie_seed)


def verify(record):
    zeros,known=single.assumptions(record['parent'],record['zero_cosets'],record['nonzero_cosets'])
    assert record['zero_indices']==sorted(zeros) and record['known_nonzero_indices']==sorted(known)
    return check(zeros,known,record['witness'])


def run(args):
    if args.verify:
        record=base.read(args.output)
        for name,digest in record['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
        print('Independent multi-pivot conditional rank',verify(record),flush=True);return
    assert not args.output.exists()
    zeros,known=single.assumptions(args.parent,args.zero,args.nonzero)
    witness,search=discover(zeros,known,args.width,args.seconds,[base.read(p) for p in args.seeds],args.start_rank,args.tie_seed)
    record=dict(status='CONDITIONAL_MULTI_PIVOT_WITNESS_NOT_FULL_DISTANCE_PROOF',parent=args.parent,
        zero_cosets=args.zero,nonzero_cosets=args.nonzero,zero_indices=sorted(zeros),
        known_nonzero_indices=sorted(known),witness=witness,search=search,
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in (Path(__file__),Path(single.__file__))},
        seed_sha256={str(p.resolve()):base.sha(p) for p in args.seeds})
    rank=verify(record);base.write_new(args.output,record)
    print(dict(parent=args.parent,zero=args.zero,nonzero=args.nonzero,rank=rank,**search),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parent',type=int,choices=[7,15],required=True)
    p.add_argument('--zero',type=int,nargs='*',default=[])
    p.add_argument('--nonzero',type=int,nargs='*',default=[])
    p.add_argument('--seeds',type=Path,nargs='*',default=[])
    p.add_argument('--start-rank',type=int,default=24)
    p.add_argument('--tie-seed',type=int,default=1)
    p.add_argument('--width',type=int,default=800)
    p.add_argument('--seconds',type=int,default=30)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
