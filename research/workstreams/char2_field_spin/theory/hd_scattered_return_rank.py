"""Joint rank bounds for inactive output zeros in the scattered return map.

The input syndromes are four rows of four coordinates. Output equations
constrain columns. MDS row elimination gives a block-triangular collection
of independent constraints. This module supplies FLOAT screening caps,
not a whole-code certificate.
"""
from functools import lru_cache
from itertools import combinations_with_replacement,permutations,product
from math import comb
import numpy as np
from hd_bucket_feedback import profiles
from hd_scattered_feedback import check_ranks


@lru_cache(None)
def minimum_check_ranks():
    return np.array([[min(check_ranks(mask)[l] for mask in range(16) if mask.bit_count()==k)
                      for l in range(17)] for k in range(5)],dtype=np.uint8)


def grid_rank_lower(occupancies,output_ranks):
    """Rank lower bound uniform over labellings of four input occupancies.

    RS guarantees MDS prefixes only in degree order, not arbitrary column
    orders. Therefore take the WORST labelling, never the best assignment.
    """
    local=minimum_check_ranks()
    heights=[sum(r>p for r in output_ranks) for p in range(4)]
    return min(sum(int(local[heights[p],order[p]]) for p in range(4))
               for order in permutations(occupancies))


@lru_cache(None)
def rank_data():
    occ=np.array([row for row,count in profiles(16,64,4)],dtype=np.uint8)
    representatives=tuple(combinations_with_replacement(range(5),4))
    lookup={row:i for i,row in enumerate(representatives)}
    zeros=np.array(tuple(product(range(5),repeat=4)),dtype=np.uint8)
    indices=np.array([lookup[tuple(sorted(row))] for row in zeros])
    local=minimum_check_ranks();bounds=np.full((len(occ),len(representatives)),16,dtype=np.uint8)
    for j,row in enumerate(representatives):
        heights=[sum(r>p for r in row) for p in range(4)]
        for order in permutations(range(4)):
            candidate=sum(local[heights[p],occ[:,order[p]]] for p in range(4))
            bounds[:,j]=np.minimum(bounds[:,j],candidate)
    rho=np.minimum(4,occ).sum(axis=1)
    if np.any(bounds>rho[:,None]):raise ArithmeticError('Rank exceeds syndrome dimension')
    return occ,zeros,bounds[:,indices],rho


@lru_cache(1)
def inactive_caps(z,q=65536.):
    occ,zeros,ranks,rho=rank_data()
    # Group zero-set sizes >=4: they impose the same four-dimensional rank.
    marks=np.zeros((17,5))
    for n in range(17):
        for e in range(n+1):marks[n,min(e,4)]+=comb(n,e)*z**(n-e)*(1-z)**e
    term=np.maximum(0,q**(-ranks.astype(float))-q**(-rho.astype(float))[:,None])
    for block in range(4):term*=marks[(16-occ[:,block])[:,None],zeros[None,:,block]]
    free=occ.sum(axis=1)-rho
    return term.sum(axis=1)*(z+(1-z)/q)**free


@lru_cache(None)
def labelled_data():
    """All labelled occupancies, with their sorted-profile aggregation."""
    labelled=np.array(tuple(product(range(17),repeat=4)),dtype=np.uint8)
    lookup={tuple(row):i for i,(row,_) in enumerate(profiles(16,64,4))}
    indices=np.array([lookup[tuple(sorted(row))] for row in labelled])
    multiplicities=np.bincount(indices)
    checks=np.array([check_ranks(mask) for mask in range(16)],dtype=np.uint8)
    rho=np.minimum(labelled,4).sum(axis=1)
    return labelled,indices,multiplicities,checks,rho


@lru_cache(1)
def averaged_inactive_caps(z,q=65536.):
    """Average natural-degree pivot bound over labelled occupancies.

    This retains which input block receives each occupancy and the exact
    set of check indices at each pivot. It does not reorder RS columns.
    The expectation is taken before the paired Jensen comparison.
    """
    occ,indices,multiplicities,checks,rho=labelled_data()
    marks=np.zeros((17,5))
    for n in range(17):
        for e in range(n+1):marks[n,min(e,4)]+=comb(n,e)*z**(n-e)*(1-z)**e
    powers=q**(-np.arange(17,dtype=float));quiet=powers[rho]
    total=np.zeros(len(occ))
    for zeros in product(range(5),repeat=4):
        rank=np.zeros(len(occ),dtype=np.uint8);weight=np.ones(len(occ))
        for p in range(4):
            mask=sum(1<<j for j in range(4) if zeros[j]>p)
            rank+=checks[mask,occ[:,p]]
            weight*=marks[16-occ[:,p],zeros[p]]
        total+=np.maximum(0,powers[rank]-quiet)*weight
    total*=(z+(1-z)/q)**(occ.sum(axis=1)-rho)
    return np.bincount(indices,weights=total)/multiplicities


if __name__=='__main__':
    import json
    occ,_,_,_=rank_data()
    for z in (.7,.5,.07):
        values=inactive_caps(z)
        print(json.dumps(dict(mark=z,selected={str(tuple(int(i) for i in row)):float(v)
              for row,v in zip(occ,values) if tuple(row) in ((0,0,0,1),(0,0,1,1),(2,2,2,2),(4,4,4,4),(6,6,6,6))})),flush=True)
