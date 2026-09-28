"""Exact rank-two composition/counting tools; no full SPIN certificate."""
from collections import Counter
from itertools import combinations,combinations_with_replacement,permutations,product
from math import comb,factorial,prod,log2

from bch_joint_support import authenticated_caps
from joint_support import span,joint_bound,cumulative_basis_counts


def composition(weights,n):
    a,b,c=weights
    if (a+b+c)%2 or max(weights)*2>sum(weights) or sum(weights)>2*n:
        return None
    u=(a+b+c)//2
    return n-u,u-b,u-a,u-c  # zero, first only, second only, both


def subspace_cap(weights,spectrum):
    """Upper count of 2-spaces with the given sorted nonzero weight triple."""
    n=len(spectrum)-1
    if composition(weights,n) is None:
        return 0
    bounds=[]
    for a,b,c in set(permutations(weights)):
        # Count ordered distinct pairs (x,y) with weights a,b and wt(x+y)=c.
        pairs=spectrum[a]*(spectrum[b]-int(a==b))
        intersection=(a+b-c)//2
        patterns=comb(a,intersection)*comb(n-a,b-intersection)
        bounds.append(min(pairs,spectrum[a]*patterns))
    # Every 2-space contributes this many pairs to each fixed weight ordering.
    multiplicity=prod(factorial(k) for k in Counter(weights).values())
    return min(bounds)//multiplicity


def embeddings():
    """Counts of labelled output-weight triples over 210 rank-two row maps."""
    result=Counter()
    for rows in product(range(4),repeat=4):
        if len(set(rows)-{0})<2:
            continue
        weights=tuple(sum((row&v).bit_count()%2 for row in rows) for v in (1,2,3))
        result[weights]+=1
    assert sum(result.values())==210
    return result


def self_test():
    for rows,n in (([1,2,4],3),([15,51,85],7),([15,3],4),([0x97,0x4b,0x2d,0x1e],8)):
        words=span(rows)
        spectrum=[0]*(n+1)
        for w in words:
            spectrum[w.bit_count()]+=1
        spaces={tuple(sorted((x,y,x^y))) for x,y in combinations(words[1:],2)}
        counts=Counter(tuple(sorted(w.bit_count() for w in space)) for space in spaces)
        for weights in combinations_with_replacement(range(1,n+1),3):
            assert counts[weights]<=subspace_cap(weights,spectrum)
        for space in spaces:
            x,y=space[:2]
            c=composition((x.bit_count(),y.bit_count(),(x^y).bit_count()),n)
            assert c==(n-(x|y).bit_count(),(x&~y).bit_count(),(y&~x).bit_count(),(x&y).bit_count())
    assert {tuple(sorted(v)) for v in embeddings()}=={(1,1,2),(1,2,3),(1,3,4),(2,2,2),(2,2,4),(2,3,3)}
    print('Rank-two composition, subspace caps, and all 210 row maps checked',flush=True)


def main():
    self_test()
    spectrum=authenticated_caps()
    weights=[w for w in range(1,257) if spectrum[w]]
    shells=[0]*257
    types=0
    for triple in combinations_with_replacement(weights,3):
        cap=subspace_cap(triple,spectrum)
        if cap:
            shells[sum(triple)//2]+=210*cap
            types+=1
    cumulative=cumulative_basis_counts(spectrum,2)
    print('Feasible types',types,flush=True)
    for u in (57,58,60,64,72,80,96,128,192,256):
        typed=sum(shells[:u+1])
        old=joint_bound(spectrum,4,2,u,cumulative=cumulative)[0]
        print('support',u,'typed log2',log2(typed) if typed else '-inf',
              'previous basis log2',log2(old) if old else '-inf',flush=True)
    print('Row-map weight triples',sorted(embeddings().items()))
    print('No output-weight bound or complete distance certificate.')


if __name__=='__main__':
    main()
