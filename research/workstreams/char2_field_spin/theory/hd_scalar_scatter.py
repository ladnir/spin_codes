"""PROPOSED fresh-scalar emission: Y=X+gamma A Q, Q'=T(Q+C_sc X).

The common-refresh density lemma is unchanged. At a return the output is
(I+gamma A C_sc)X. A determinant-polynomial bound retains active output
constraints outside at most min(c,16) exceptional scalars. No timing or
whole-code certificate is claimed. FLOAT screening only.
"""
from functools import lru_cache
from itertools import product
from math import comb
import numpy as np
from hd_scattered_return_rank import labelled_data


@lru_cache(1)
def components(z,q=65536.):
    occ,indices,multiplicities,checks,rho=labelled_data()
    marks=np.zeros((17,5));active=np.zeros((16,17))
    quiet_local=np.zeros(17)
    for n in range(17):
        for e in range(n+1):marks[n,min(e,4)]+=comb(n,e)*z**(n-e)*(1-z)**e
        for c in range(n+1):
            weight=comb(n,c)*z**(n-c)*(1-z)**c*q**(-c)
            quiet_local[n]+=weight*q**(-min(4,n-c))
            for mask in range(16):active[mask,n]+=weight*q**(-float(checks[mask,n-c]))
    powers=q**(-np.arange(17,dtype=float));quiet=powers[rho]
    quiet_mark=np.prod(quiet_local[occ],axis=1)
    inactive=np.zeros(len(occ));generic=np.zeros(len(occ))
    for zeros in product(range(5),repeat=4):
        rank=np.zeros(len(occ),dtype=np.uint8);weight=np.ones(len(occ));full=np.ones(len(occ))
        for p in range(4):
            mask=sum(1<<j for j in range(4) if zeros[j]>p)
            rank+=checks[mask,occ[:,p]]
            weight*=marks[16-occ[:,p],zeros[p]]
            full*=active[mask,occ[:,p]]
        inactive+=np.maximum(0,powers[rank]-quiet)*weight
        generic+=np.maximum(0,full-quiet_mark)*weight
    return (np.bincount(indices,weights=inactive)/multiplicities,
            np.bincount(indices,weights=generic)/multiplicities)


def return_caps(z,scalar_order=65536,q=65536.):
    from hd_scattered_return_rank import rank_data
    occ,_,_,rho=rank_data();degrees=occ.sum(axis=1)
    inactive,generic=components(z,q)
    exceptions=np.array([sum(comb(j,c)*z**(j-c)*(1-z)**c*min(c,16)
                             for c in range(j+1))/(scalar_order-1) for j in range(65)])
    # The fiber information-set cap is valid for every nonzero gamma.
    fiber=inactive*(z+(1-z)/q)**(degrees-rho)
    return np.minimum(fiber,generic+exceptions[degrees]*inactive)
