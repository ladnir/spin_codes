"""Exact support-CDF caps from all hyperplanes and short coset representatives.

Counts subspaces through their codimension-one subspaces. This is an outer
counting bound only, not a full SPIN failure certificate.
"""
from fractions import Fraction
from itertools import accumulate,product
from math import prod,log2

from joint_support import span,gf2_rank


def extension_bound(spectrum,h,u,v):
    """Upper extension count for a fixed (h-1)-space of support v.

    Every h-space extending it with support <=u has 2^(h-1) coset words.
    Their average weight is <=u-v/2. Each short word determines its extension.
    """
    possible=[w for w in range(1,u+1) if spectrum[w]]
    if not possible:
        return Fraction(0)
    d=possible[0]
    mean=Fraction(2*u-v,2)
    if mean<d:
        return Fraction(0)
    size=1<<(h-1)
    cumulative=list(accumulate(spectrum))
    best=None
    for i,t in enumerate(possible):
        if i+1==len(possible):
            good=size
        else:
            nxt=possible[i+1]
            if nxt<=mean:
                continue
            lower=size*(nxt-mean)/(nxt-d)
            good=-(-lower.numerator//lower.denominator)
        assert 1<=good<=size
        # The zero word is never a coset representative.
        candidate=Fraction(cumulative[t]-spectrum[0],good)
        best=candidate if best is None else min(best,candidate)
    assert best is not None
    return best


def refine(caps,spectrum,verbose=False):
    g=len(caps)
    result=[row[:] for row in caps]
    embeddings=[prod((1<<g)-(1<<j) for j in range(h)) for h in range(g+1)]
    for h in range(2,g+1):
        lower=result[h-2]
        for u in range(len(spectrum)):
            weights=[extension_bound(spectrum,h,u,v) for v in range(u+1)]
            assert all(a>=b for a,b in zip(weights,weights[1:]))
            # CDF domination, not a claim that increments bound shell counts.
            weighted=sum(((lower[v]-(lower[v-1] if v else 0))*weights[v]
                          for v in range(u+1)),Fraction(0))
            bound=weighted*Fraction(embeddings[h],((1<<h)-1)*embeddings[h-1])
            result[h-1][u]=min(result[h-1][u],bound.numerator//bound.denominator)
        for u in range(len(spectrum)-2,-1,-1):
            result[h-1][u]=min(result[h-1][u],result[h-1][u+1])
        if verbose:
            print('Hyperplane rank',h,'support, old bits, new bits',
                  [(u,round(log2(caps[h-1][u]),4),round(log2(result[h-1][u]),4))
                   for u in (72,76,80,83,84,96,128,192,256) if result[h-1][u]],flush=True)
    return result


def self_test():
    from basis_lattice import bound as basis_bound
    checks=0
    for rows,n,g in (([1,2,4],3,3),([15,51,85],7,3),([15,3],4,3),([0x97,0x4b,0x2d,0x1e],8,4)):
        words=span(rows)
        spectrum=[sum(w.bit_count()==u for w in words) for u in range(n+1)]
        shells=[[0]*(n+1) for _ in range(g+1)]
        for values in product(words,repeat=g):
            union=0
            for value in values:
                union|=value
            shells[gf2_rank(values)][union.bit_count()]+=1
        caps=[[basis_bound(spectrum,g,h,u) for u in range(n+1)] for h in range(1,g+1)]
        tightened=refine(caps,spectrum)
        for h in range(1,g+1):
            actual=list(accumulate(shells[h]))
            assert all(a<=b for a,b in zip(actual,tightened[h-1]))
            checks+=len(actual)
    print('Hyperplane-CDF exact small-code inequalities passed:',checks,flush=True)


if __name__=='__main__':
    from bch_joint_support import authenticated_caps,support_caps
    from shortened_bound import dimension_caps
    from basis_lattice import improve_caps
    self_test()
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    refine(caps,spectrum,True)
