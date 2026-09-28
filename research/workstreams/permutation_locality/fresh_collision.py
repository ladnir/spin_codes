"""Exact fresh-state cancellation from distinct-window kernel counts."""
from collections import Counter
from fractions import Fraction
from itertools import permutations,product
from math import comb,prod
from flint import arb
from occupancy_memory import F,Z
from group_rank_one_verify import up


def pair_completions(a,b,d,width=4):
    # For a fixed D of weight d, count A,B with A xor B=D and weights a,b.
    if (a-b+d)%2:
        return 0
    inside=(a-b+d)//2;outside=(a+b-d)//2
    return comb(d,inside)*comb(width-d,outside) if 0<=inside<=d and 0<=outside<=width-d else 0


def probability(zeros,a,weights,windows=32,width=4):
    """Requires injectivity on every pair of windows; j is two or three."""
    j=len(weights);assert j in (2,3)
    numerator=zeros.get(j+1,{}).get((a,)+tuple(weights),0)
    for i,b in enumerate(weights):
        for d in range(1,width+1):
            shape=list(weights);shape[i]=d
            numerator+=pair_completions(a,b,d,width)*zeros.get(j,{}).get(tuple(shape),0)
    denominator=windows*comb(width,a)*prod(range(windows-j+1,windows+1))*prod(comb(width,b) for b in weights)
    return Fraction(numerator,denominator)


def refine(base,prepared,tilt,penalty,input_penalty=1,odd_penalty=1):
    single,pair,zeros=prepared[0]
    assert all(row[1]==0 for row in pair[1].values())
    tables={1:{},2:{},**zeros}
    states=(1<<19)-1
    for j in (2,3):
        assert j+1 in tables
        candidates=[]
        for weights in product(range(1,5),repeat=j):
            cancellation=max(probability(tables,a,weights) for a in range(1,5))
            factor=(-arb(tilt)*max(0,48-sum(weights))).exp()*arb(penalty)**weights.count(4)*arb(input_penalty)**sum(weights)*arb(odd_penalty)**sum(w%2 for w in weights)
            value=factor*(arb(cancellation.numerator)/cancellation.denominator+arb(1)/states)/2
            candidates.append(up(value))
        base[j][F,Z]=min(base[j][F,Z],max(candidates))
    return base


def self_test():
    for a,b,d in product(range(1,5),range(1,5),range(5)):
        target=(1<<d)-1
        observed=sum(mask.bit_count()==a and (mask^target).bit_count()==b for mask in range(16))
        assert observed==pair_completions(a,b,d)
    width=2;windows=3
    columns=((1,2),(4,8),(5,10))
    images={}
    for i,row in enumerate(columns):
        for mask in range(1,1<<width):
            value=0
            for bit,column in enumerate(row):
                if mask>>bit&1:value^=column
            images[i,mask]=value
    zeros={j:Counter() for j in range(1,5)}
    for j in range(1,4):
        for positions in permutations(range(windows),j):
            for masks in product(range(1,1<<width),repeat=j):
                value=0
                for i,m in zip(positions,masks):value^=images[i,m]
                if not value:zeros[j][tuple(m.bit_count() for m in masks)]+=1
    assert not zeros[1] and not zeros[2]
    checks=0
    for j in (2,3):
        actual=Counter()
        for old in images:
            for positions in permutations(range(windows),j):
                for masks in product(range(1,1<<width),repeat=j):
                    value=images[old]
                    for i,m in zip(positions,masks):value^=images[i,m]
                    if not value:actual[old[1].bit_count(),tuple(m.bit_count() for m in masks)]+=1
        for a in range(1,width+1):
            for weights in product(range(1,width+1),repeat=j):
                denominator=windows*comb(width,a)*prod(range(windows-j+1,windows+1))*prod(comb(width,b) for b in weights)
                assert probability(zeros,a,weights,windows,width)==Fraction(actual[a,weights],denominator)
                checks+=1
    print('Fresh overlap/distinct cancellation checks:',checks,'exact shape identities',flush=True)


if __name__=='__main__':self_test()
