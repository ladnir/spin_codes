"""Proof screen helpers for shifted and coordinate-scattered feedback.

The permutation of sixteen feedback coordinates replaces a native layout
transpose. This is a construction change, not an implementation-only win.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb
import numpy as np
from hd_bucket_feedback import maps as block_maps
from hd_systematic_mds import novel
from hd_common_refresh_enumerator import extend,field
from hd_blockwise_shared_transfer import mm


@lru_cache(None)
def maps():
    _,c0=block_maps(local_rank=4);a=tuple(zip(*c0))
    c=tuple(tuple(novel((i%16)^16,r//4) if i//16==r%4 else 0 for i in range(64)) for r in range(16))
    if any(x for row in mm(c,a,65536) for x in row):raise ArithmeticError('Nonzero C A')
    return a,c


@lru_cache(None)
def check_ranks(mask):
    indices=[i for i in range(4) if mask>>i&1]
    columns=tuple(tuple(novel(j^16,i) for i in indices) for j in range(16))
    spaces={()}
    for col in columns:spaces.update(extend(basis,col) for basis in tuple(spaces))
    maximum=[0]*(len(indices)+1)
    for basis in spaces:
        maximum[len(basis)]=max(maximum[len(basis)],sum(extend(basis,col)==basis for col in columns))
    return tuple(next(r for r,n in enumerate(maximum) if n>=l) for l in range(17))


@lru_cache(None)
def small_return_counts():
    return count_small_returns(*maps())


@lru_cache(None)
def count_small_returns(a,c):
    ac=mm(a,c,65536);p=tuple(tuple(x^int(i==j) for j,x in enumerate(row)) for i,row in enumerate(ac))
    columns=tuple(zip(*p));table,inv=field();q=65536
    result={};totals=Counter()
    for j in (1,2):
        for u in combinations(range(64),j):
            occ=tuple(sorted(sum(i//16==b for i in u) for b in range(4)))
            weights=result.setdefault(occ,Counter());totals[occ]+=1
            if j==1:weights[sum(x!=0 for x in columns[u[0]])]+=q-1;continue
            aa,bb=(columns[i] for i in u);ratios=Counter(table[x][inv[y]] for x,y in zip(aa,bb) if y)
            generic=sum(bool(x or y) for x,y in zip(aa,bb))
            weights[generic]+=(q-len(ratios))*(q-1)
            weights[sum(y!=0 for y in bb)]+=q-1
            for zeros in ratios.values():weights[generic-zeros]+=q-1
    for occ,counts in result.items():
        if sum(counts.values())!=totals[occ]*(q**sum(occ)-1):raise ArithmeticError('Incomplete projective line count')
    return {occ:(dict(counts),totals[occ]*q**sum(occ)) for occ,counts in result.items()}


def return_caps(occ,z,q=65536.):
    rho=np.minimum(4,occ).sum(axis=1);blocks=4;ret=np.zeros(len(occ))
    # The labelled occupancy is uniform over permutations of its profile.
    # Average the output floors before the paired Jensen comparison.
    elementary=np.zeros((len(occ),5));elementary[:,0]=1
    for b in range(4):
        v=z**np.maximum(0,13-occ[:,b])
        for w in range(b+1,0,-1):elementary[:,w]+=elementary[:,w-1]*v
    for active in range(1,16):
        ranks=np.array(check_ranks(active^15));rank=ranks[occ].sum(axis=1)
        probability=np.maximum(0,q**(-rank.astype(float))-q**(-rho.astype(float)))
        ret+=probability*elementary[:,active.bit_count()]/comb(4,active.bit_count())
    for pattern,(counts,denominator) in small_return_counts().items():
        match=np.all(occ==pattern,axis=1)
        # For <=2 inputs, C_U is injective, so all nonzero words return.
        value=sum(count*z**weight for weight,count in counts.items())/denominator
        ret[match]=np.minimum(ret[match],value)
    return ret


if __name__=='__main__':
    import json
    for mask in range(16):print(json.dumps(dict(check_mask=mask,minimum_ranks=check_ranks(mask))))
    print(json.dumps({str(k):v for k,v in small_return_counts().items()}))
