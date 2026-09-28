"""Exact histogram of feedback dual characters across the 32 windows."""
from collections import Counter
from math import comb,factorial,prod
from fractions import Fraction
from itertools import product
import numpy as np
from group_moment import maps


def census(with_character_index=False):
    _,columns,_=maps()
    characters=np.arange(1<<19,dtype=np.uint32)
    counts=np.zeros((len(characters),5),dtype=np.uint8)
    for window in range(32):
        weight=np.zeros(len(characters),dtype=np.uint8)
        for bit in range(4):
            weight+=np.bitwise_count(characters&columns[4*window+bit])&1
        for r in range(5):counts[:,r]+=weight==r
    assert (counts.sum(axis=1)==32).all()
    keys=sum(counts[:,r].astype(np.uint64)*33**r for r in range(5))
    values,inverse,multiplicities=np.unique(keys,return_inverse=True,return_counts=True)
    result={tuple(int(key//33**r%33) for r in range(5)):int(n)
            for key,n in zip(values,multiplicities)}
    assert sum(result.values())==1<<19
    # Independent scalar checks on selected characters.
    for character in (0,1,17,131071,262144,524287):
        direct=Counter(sum((character&columns[4*w+b]).bit_count()%2 for b in range(4)) for w in range(32))
        assert tuple(direct[r] for r in range(5))==tuple(counts[character])
    print('Feedback character window histograms:',len(result),'cover',sum(result.values()),'characters',flush=True)
    print('All-four feedback parity weight distribution:',sorted(Counter({w:sum(n for h,n in result.items() if h[1]+h[3]==w) for w in range(33)}).items()),flush=True)
    return (result,inverse) if with_character_index else result


def character_polynomials(maximum=6):
    """Compressed character polynomials and their exact shape denominators.

    Average each dual character polynomial over all 2^19 characters. The
    coefficient is an integer count by character orthogonality. All integer
    bounds below are checked before vectorized int64 arithmetic.
    """
    assert 1<=maximum<=10
    records,inverse=census(with_character_index=True)
    histograms=np.array(list(records),dtype=np.int64)
    multiplicities=np.array(list(records.values()),dtype=np.int64)
    shapes=sorted((s for s in product(range(maximum+1),repeat=4) if sum(s)<=maximum),key=lambda s:(sum(s),s))
    index={s:i for i,s in enumerate(shapes)}
    # Krawtchouk sums for a four-bit window containing r negative signs.
    kraw=np.array([[sum((-1)**v*comb(r,v)*comb(4-r,b-v)
                         for v in range(max(0,b-(4-r)),min(b,r)+1))
                   for b in range(1,5)] for r in range(5)],dtype=np.int64)
    denoms={s:comb(32,sum(s))*factorial(sum(s))//prod(factorial(x) for x in s)
               *prod(comb(4,b)**s[b-1] for b in range(1,5)) for s in shapes}
    # A partial polynomial coefficient is a signed count of a subset of
    # these assignments. This bounds every partial update, not just the
    # final coefficient. Walsh inversion needs its own L1 bound below.
    assert max(denoms.values())<1<<63
    coefficients=np.zeros((len(shapes),len(records)),dtype=np.int64)
    coefficients[index[(0,0,0,0)]]=1
    cumulative=np.cumsum(histograms,axis=1)
    for window in range(32):
        types=np.sum(cumulative<=window,axis=1)
        for s in reversed(shapes):
            if not 1<=sum(s)<=window+1:continue
            for b in range(4):
                if s[b]:
                    prior=list(s);prior[b]-=1
                    coefficients[index[s]]+=coefficients[index[tuple(prior)]]*kraw[types,b]
    return shapes,coefficients,denoms,multiplicities,inverse


def zero_counts(maximum=6):
    shapes,coefficients,denoms,multiplicities,_=character_polynomials(maximum)
    result={}
    for s,row in zip(shapes,coefficients):
        numerator=sum(int(x)*int(m) for x,m in zip(row,multiplicities))
        assert numerator>=0 and numerator%(1<<19)==0
        numerator//=1<<19
        assert numerator<=denoms[s]
        weights=tuple(b for b in range(1,5) for _ in range(s[b-1]))
        ordered=numerator*prod(factorial(x) for x in s)
        result[weights]=(ordered,Fraction(numerator,denoms[s]))
    for j in range(1,maximum+1):
        print('Exact character zero feedback',j,'worst',max((p,w) for w,(_,p) in result.items() if len(w)==j),flush=True)
    return result


def self_test(result):
    from occupancy_model import local_data
    data=local_data(4);checked=0
    for j,rows in data[2].items():
        for weights in product(range(1,5),repeat=j):
            assert result[tuple(sorted(weights))][0]==rows.get(weights,0)
            checked+=1
    print('Character vs independent direct-rank zero counts:',checked,'ordered shapes',flush=True)


if __name__=='__main__':self_test(zero_counts())
