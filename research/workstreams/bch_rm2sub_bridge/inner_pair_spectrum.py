"""Exact pair spectrum of the fixed A code via XOR convolution.

P[a,b,c] counts ordered (u,v) with wt(Au)=a, wt(Av)=b,
wt(A(u+v))=c. Integer Walsh transforms avoid 2^(2s) enumeration.
"""
import argparse
from collections import Counter
from pathlib import Path
import numpy as np
import bridge as base


def walsh(values):
    result=values.copy();n=result.shape[-1];assert n and n&(n-1)==0
    for h in (1<<k for k in range(n.bit_length()-1)):
        blocks=result.reshape(*result.shape[:-1],-1,2*h)
        left=blocks[...,:h].copy();right=blocks[...,h:]
        blocks[...,:h]=left+right
        blocks[...,h:]=left-right
    return result


def pair_counts(generators):
    n=1<<len(generators);words=[0]*n
    for u in range(1,n):
        bit=u&-u;words[u]=words[u^bit]^generators[bit.bit_length()-1]
    weights=np.array([word.bit_count() for word in words],dtype=np.int64)
    shells=sorted(set(map(int,weights)));indicators=np.array([weights==w for w in shells],dtype=np.int64)
    transformed=walsh(indicators);counts={}
    # Worst absolute accumulation <= n^4. For s<=15 this fits signed int64.
    assert n**4<1<<63
    for i,a in enumerate(shells):
        for j,b in enumerate(shells):
            product=transformed[i]*transformed[j]
            for k,c in enumerate(shells):
                numerator=int(np.sum(product*transformed[k],dtype=np.int64))
                value,remainder=divmod(numerator,n);assert remainder==0 and value>=0
                if value:counts[a,b,c]=value
    return counts,Counter(map(int,weights))


def run(verify=False):
    name='t128_s15';t,s,spectrum=base.load_map(name)
    selection=base.HERE/'inputs'/f'{name}_selection.json'
    generators=[int(x,16) for x in base.read(selection)['selected']['A_generator_words_hex']]
    counts,marginal=pair_counts(generators);n=1<<s
    assert marginal==Counter({0:1,**spectrum}) and sum(counts.values())==n*n
    for (a,b,c),value in counts.items():
        assert (a+b+c)%2==0
        cells=((2*t-a-b-c)//2,(b+c-a)//2,(a+c-b)//2,(a+b-c)//2)
        assert min(cells)>=0 and sum(cells)==t
        assert counts[a,c,b]==counts[b,a,c]==value
    for a,count in marginal.items():
        assert sum(v for (x,y,z),v in counts.items() if x==a)==n*count
        assert counts[0,a,a]==counts[a,0,a]==counts[a,a,0]==count
    path=base.HERE/'generated/t128_s15_pair_spectrum.json'
    result=dict(status='EXACT_FIXED_A_PAIR_SPECTRUM',configuration=name,length=t,dimension=s,
        rows=[dict(first_weight=a,second_weight=b,sum_weight=c,ordered_pairs=v) for (a,b,c),v in sorted(counts.items())],
        total_ordered_pairs=n*n,marginal={str(w):v for w,v in sorted(marginal.items())},
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),selection,
            base.HERE/'inputs/t128_s15_a_spectrum.json',base.HERE/'inputs/manifest.json',Path(base.__file__)]})
    if verify:assert result==base.read(path)
    else:base.write_new(path,result)
    print('Exact A pair spectrum:',len(counts),'nonzero triples;',n*n,'ordered pairs; verified' if verify else 'ordered pairs',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
