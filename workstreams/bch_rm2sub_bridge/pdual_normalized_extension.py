"""Search a two-term final append using the normalization F7=1.

Unlike the monomial rank search, the modified append can be v7+vz, z in the
zero set. A following transform may cancel two known unit Fourier values.
"""
import argparse
from pathlib import Path
import time
import bridge as base
import pdual_rank_refinement as rank
import pdual_multi_pivot as multi


def check(witness):
    zeros,ones=rank.assumptions(7,[],[])
    state=rank.check_path(zeros,7,witness['prefix_steps'])
    assert len(state)==29
    k,s=witness['penultimate_transform']
    actual=(s+(1-(1<<k))*7)%255
    moved={(e*(1<<k)+actual)%255 for e in state}
    assert moved<=zeros and len(moved)==29
    z=witness['zero_append_exponent'];assert z in zeros
    rows=[{e} for e in sorted(moved)]+[{7,z}]
    k,s=witness['final_transform']
    assert 0<=k<8 and 0<=s<255
    rows=[{(e*(1<<k)+s)%255 for e in row} for row in rows]
    assert all(row<=zeros|ones and len(row&ones)%2==0 for row in rows)
    # All rows have known sum zero; v7 has sum one and raises rank to 31.
    return 31


def discover(seed,width,seconds):
    zeros,ones=rank.assumptions(7,[],[])
    mask=(1<<255)-1;zword=sum(1<<i for i in zeros);oword=sum(1<<i for i in ones)
    allowed=[multi.rotate(zword,(-i)%255) for i in range(255)]
    powers=[[i*(1<<k)%255 for i in range(255)] for k in range(8)]
    def moves(state):
        out=[];inds=list(multi.positions(state))
        for k,power in enumerate(powers):
            choices=mask;word=0
            for e in inds:
                j=power[e];choices &= allowed[j];word |= 1<<j
                if not choices:break
            if choices:out.append((k,word,choices))
        return out
    state={7};path=[];seeds={1:{1<<7:[]}}
    for k,s in seed['steps']:
        a=(s+(1-(1<<k))*7)%255
        state={(e*(1<<k)+a)%255 for e in state}|{7};path=path+[(k,s)]
        assert len(state)==len(path)+1
        seeds[len(state)]={sum(1<<e for e in state):path}
    beam=seeds[24];begun=time.monotonic();levels=[];tried=0
    for size in range(25,31):
        candidates=dict(seeds.get(size,{}))
        for word,path in beam.items():
            for k,moved,choices in moves(word):
                for s in multi.positions(choices):
                    new=multi.rotate(moved,s)|(1<<7)
                    relative=(s-(1-(1<<k))*7)%255
                    if new not in candidates:candidates[new]=path+[(k,relative)]
            if time.monotonic()-begun>seconds:break
        levels.append(dict(rank=size,candidates=len(candidates)))
        if size==30:
            for word,path in candidates.items():
                tried+=1
                for k,moved,choices in moves(word^(1<<7)):
                    possible=choices&multi.rotate(oword,(-powers[k][7])%255)
                    for s in multi.positions(possible):
                        inverse=pow(1<<k,-1,255)
                        choices_z=sorted({((o-s)*inverse)%255 for o in ones}&zeros)
                        if choices_z:
                            witness=dict(prefix_steps=path[:-1],penultimate_transform=path[-1],
                                zero_append_exponent=choices_z[0],final_transform=[k,s])
                            check(witness)
                            return dict(found=True,witness=witness,rank=31,tried=tried,levels=levels,
                                        elapsed_seconds=time.monotonic()-begun)
                if time.monotonic()-begun>seconds:break
            break
        ranked=[]
        for word,path in candidates.items():
            score=sum(s.bit_count() for _,_,s in moves(word))
            if score:
                tie=((word^(word>>64)^(word>>128))*0x9e3779b97f4a7c15)&((1<<64)-1)
                ranked.append((score,tie,word,path))
            if time.monotonic()-begun>seconds:break
        if not ranked or time.monotonic()-begun>seconds:break
        ranked.sort(reverse=True);beam={word:path for _,_,word,path in ranked[:width]}
        beam.update(seeds.get(size,{}))
    return dict(found=False,tried=tried,levels=levels,elapsed_seconds=time.monotonic()-begun,
                proof_of_impossibility=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--width',type=int,default=1600)
    p.add_argument('--seconds',type=int,default=45)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists()
    result=discover(base.read(a.seed),a.width,a.seconds)
    result['source_sha256']={x.relative_to(base.ROOT).as_posix():base.sha(x) for x in
        (Path(__file__),Path(rank.__file__),Path(multi.__file__))}
    result['seed_sha256']=base.sha(a.seed)
    base.write_new(a.output,result)
    print({k:v for k,v in result.items() if k not in ('source_sha256','witness')},flush=True)
