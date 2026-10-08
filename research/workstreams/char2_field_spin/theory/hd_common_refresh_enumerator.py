"""Exact averaged local return enumerators, using the lattice of zero sets.

For each potential support U of size <=4, enumerate the row-space lattice
of (I+A A^T)[:,U] over GF256. Ranks remain unchanged over GF65536. Count
vectors in each exact zero set by subtraction over this finite lattice.
There are no randomized tests or FLOAT operations in the enumerator.
"""
from functools import lru_cache
from itertools import combinations
from math import comb
from hd_bucket_feedback import maps
from hd_blockwise_shared_transfer import multiply,mm


# Exact proof data from the complete support/lattice enumeration. Rebuild
# with this module's CLI; these are integer counts, not timing samples.
PINNED_COUNTS={
    1:(16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,196605,851955),
    2:(120,0,0,0,0,0,0,0,0,0,0,0,0,0,1834980,122157240,515272083180),
    3:(560,0,0,0,0,0,0,0,0,0,0,131070,1310700,72285105,4177921785,38474328366165,157587508377951975),
    4:(1820,0,0,0,0,0,0,0,0,0,131070,5701545,17790589845,124418328570,937540799721720,8194675795768190925,33564878600672601275625),
}


def pinned_counts(dimension):
    values=PINNED_COUNTS[dimension];denominator=comb(16,dimension)*65536**dimension
    if sum(values)!=denominator:raise ArithmeticError('Invalid pinned enumerator total')
    return values,denominator


@lru_cache(1)
def field():
    table=tuple(tuple(multiply(i,j,65536) for j in range(256)) for i in range(256))
    inv=(0,)+tuple(table[i].index(1) for i in range(1,256))
    return table,inv


def extend(basis,vector):
    table,inv=field();v=list(vector)
    for pivot,row in basis:
        factor=v[pivot]
        if factor:v=[x^table[factor][y] for x,y in zip(v,row)]
    if not any(v):return basis
    pivot=next(i for i,x in enumerate(v) if x);scale=inv[v[pivot]]
    v=tuple(table[scale][x] for x in v)
    rows=[]
    for p,row in basis:
        factor=row[pivot]
        rows.append((p,tuple(x^table[factor][y] for x,y in zip(row,v))))
    rows.append((pivot,v));return tuple(sorted(rows))


def enumerator(rows,dimension,q=65536):
    spaces={()}
    for row in rows:
        spaces.update(extend(basis,row) for basis in tuple(spaces) if len(basis)<dimension)
    flats={}
    for basis in spaces:
        mask=sum(1<<i for i,row in enumerate(rows) if extend(basis,row)==basis)
        if mask in flats and flats[mask]!=len(basis):raise ArithmeticError('Inconsistent flat rank')
        flats[mask]=len(basis)
    exact={};weights=[0]*(len(rows)+1)
    for mask,rank in sorted(flats.items(),key=lambda item:-item[1]):
        count=q**(dimension-rank)-sum(n for sup,n in exact.items() if sup&mask==mask)
        if count<0:raise ArithmeticError('Negative lattice count')
        exact[mask]=count;weights[len(rows)-mask.bit_count()]+=count
    if sum(weights)!=q**dimension:raise ArithmeticError('Incomplete lattice')
    return tuple(weights)


@lru_cache(None)
def average_counts(dimension,q=65536):
    if not 0<=dimension<=4:raise ValueError('Only the small-support local return is enumerated')
    _,c=maps(local_rank=4);local_c=tuple(row[:16] for row in c[:4]);a=tuple(zip(*local_c))
    ac=mm(a,local_c,65536)
    p=tuple(tuple(x^int(i==j) for j,x in enumerate(row)) for i,row in enumerate(ac))
    total=[0]*17
    for u in combinations(range(16),dimension):
        row=tuple(tuple(p[i][j] for j in u) for i in range(16))
        for w,n in enumerate(enumerator(row,dimension,q)):total[w]+=n
    denominator=comb(16,dimension)*q**dimension
    if sum(total)!=denominator:raise ArithmeticError('Incomplete support average')
    return tuple(total),denominator


if __name__=='__main__':
    import argparse,json,time
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dimensions',type=int,nargs='+',default=(1,2,3,4))
    args=p.parse_args()
    for l in args.dimensions:
        start=time.perf_counter();values,denominator=average_counts(l)
        if l in PINNED_COUNTS and (values,denominator)!=pinned_counts(l):raise ArithmeticError('Pinned enumerator differs from replay')
        print(json.dumps(dict(dimension=l,counts=values,denominator=denominator,seconds=time.perf_counter()-start)),flush=True)
