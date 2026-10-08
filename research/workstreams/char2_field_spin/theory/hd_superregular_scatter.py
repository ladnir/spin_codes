"""PROPOSED Lagrange-basis variant; not the timed mode36 construction.

Each local map evaluates the four Lagrange polynomials for nodes16..19
at nodes0..15. Every square minor is nonzero. Scattering the feedback
coordinates preserves C A=0. This module only supplies FLOAT proof screens.
"""
from functools import lru_cache
from itertools import combinations_with_replacement,permutations,product
from math import comb
import numpy as np
from hd_common_refresh_enumerator import field
from hd_bucket_feedback import profiles
from hd_blockwise_shared_transfer import mm
from hd_scattered_feedback import count_small_returns


@lru_cache(None)
def local_map():
    table,inv=field();rows=[]
    for x in range(16):
        row=[]
        for j in range(4):
            value=1
            for k in range(4):
                if j!=k:value=table[value][table[x^(16+k)][inv[j^k]]]
            row.append(value)
        rows.append(tuple(row))
    return tuple(rows)


@lru_cache(None)
def maps():
    local=local_map()
    a=tuple(tuple(local[i%16][j%4] if i//16==j//4 else 0 for j in range(16)) for i in range(64))
    c=tuple(tuple(local[i%16][r//4] if i//16==r%4 else 0 for i in range(64)) for r in range(16))
    if any(x for row in mm(c,a,65536) for x in row):raise ArithmeticError('Nonzero C A')
    return a,c


@lru_cache(None)
def rank_data():
    occ=np.array([row for row,count in profiles(16,64,4)],dtype=np.uint8)
    representatives=tuple(combinations_with_replacement(range(5),4))
    lookup={row:i for i,row in enumerate(representatives)}
    zeros=np.array(tuple(product(range(5),repeat=4)),dtype=np.uint8)
    indices=np.array([lookup[tuple(sorted(row))] for row in zeros])
    bounds=np.zeros((len(occ),len(representatives)),dtype=np.uint8)
    for j,row in enumerate(representatives):
        heights=[sum(r>p for r in row) for p in range(4)]
        for order in permutations(range(4)):
            candidate=sum(np.minimum(heights[p],occ[:,order[p]]) for p in range(4))
            bounds[:,j]=np.maximum(bounds[:,j],candidate)
    rho=np.minimum(4,occ).sum(axis=1)
    if np.any(bounds>rho[:,None]):raise ArithmeticError('Rank exceeds syndrome dimension')
    return occ,zeros,bounds[:,indices],rho


@lru_cache(1)
def return_caps(z,q=65536.):
    occ,zeros,ranks,rho=rank_data()
    marks=np.zeros((17,5))
    for n in range(17):
        for e in range(n+1):marks[n,min(e,4)]+=comb(n,e)*z**(n-e)*(1-z)**e
    term=np.maximum(0,q**(-ranks.astype(float))-q**(-rho.astype(float))[:,None])
    for block in range(4):term*=marks[(16-occ[:,block])[:,None],zeros[None,:,block]]
    result=term.sum(axis=1)*(z+(1-z)/q)**(occ.sum(axis=1)-rho)

    elementary=np.zeros((len(occ),5));elementary[:,0]=1
    for b in range(4):
        v=z**np.maximum(0,13-occ[:,b])
        for w in range(b+1,0,-1):elementary[:,w]+=elementary[:,w-1]*v
    floor_cap=np.zeros(len(occ))
    for w in range(1,5):
        rank=np.minimum(4-w,occ).sum(axis=1)
        floor_cap+=(q**(-rank.astype(float))-q**(-rho.astype(float)))*elementary[:,w]
    result=np.minimum(result,floor_cap)
    for pattern,(counts,denominator) in count_small_returns(*maps()).items():
        exact=sum(count*z**weight for weight,count in counts.items())/denominator
        match=np.all(occ==pattern,axis=1)
        result[match]=np.minimum(result[match],exact)
    return result


if __name__=='__main__':
    import json
    print(json.dumps({str(k):v for k,v in count_small_returns(*maps()).items()}))
