"""An explicit iid envelope for uniformly shuffled independent categories.

For m categories and n slots, the universal density loss is the inverse
multinomial probability at its mean, maximized at balanced integer counts.
Positive category tilts preserve pointwise domination and can make the
envelope much sharper than the un-tilted mean distribution.
"""
from fractions import Fraction as Q
from math import factorial
from functools import lru_cache


def density_loss(n,categories=3):
    if type(n) is not int or n<1 or type(categories) is not int or categories<1:
        raise ValueError('positive integer slot and category counts required')
    counts=[n//categories+int(i<n%categories) for i in range(categories)]
    value=Q(n**n,factorial(n))
    for a in counts:value*=Q(factorial(a),a**a)
    return value


@lru_cache(maxsize=4096)
def capped_density_loss(n,caps):
    """Universal density loss when each category count has a known cap.

    The factor a!/a^a is log-concave on integers. Allocating units to the
    currently smallest uncapped count maximizes the product exactly.
    """
    if type(n) is not int or n<1 or not caps or any(type(c) is not int or c<0 for c in caps) or sum(caps)<n:
        raise ValueError('feasible integer category caps required')
    counts=[0]*len(caps)
    for _ in range(n):
        j=min((j for j,c in enumerate(caps) if counts[j]<c),key=lambda j:counts[j])
        counts[j]+=1
    value=Q(n**n,factorial(n))
    for a in counts:value*=Q(factorial(a),a**a)
    return value


def tilted_reference(probabilities,counts,tilt):
    """Return (reference law, pointwise multiplicative density upper).

    Each row of probabilities is a categorical law with the indicated
    count of independent slots. The slots are then uniformly permuted.
    All computations are exact rationals.
    """
    probabilities=[tuple(map(Q,row)) for row in probabilities];tilt=tuple(map(Q,tilt))
    if (not probabilities or len(probabilities)!=len(counts) or not tilt or min(tilt)<=0
            or any(type(c) is not int or c<0 for c in counts) or sum(counts)<1
            or any(len(row)!=len(tilt) or min(row)<0 or sum(row)!=1 for row in probabilities)):
        raise ValueError('categorical laws, integer slot counts and positive tilt required')
    n=sum(counts);normalizers=[sum(p*t for p,t in zip(row,tilt)) for row in probabilities]
    weights=[sum(Q(c,n)*row[j]/z for c,row,z in zip(counts,probabilities,normalizers)) for j in range(len(tilt))]
    scale=sum(weights);factor=density_loss(n,len(tilt))*scale**n
    for c,z in zip(counts,normalizers):factor*=z**c
    return tuple(w/scale for w in weights),factor
