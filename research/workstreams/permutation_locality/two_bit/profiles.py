"""Positive exact (union, intersection) enumerator for independently shuffled pairs.

Input spectrum caps give individual profile caps, not only cumulative caps.
This retains the two-bit profile needed by a future profile-resolved bound.
"""
from fractions import Fraction as Q
from math import comb,lcm


def pair_profiles(spectrum):
    """Include the zero pair. Key (u,b) has a=u-b weight-one packets."""
    if not spectrum or any(type(x) is not int or x<0 for x in spectrum):
        raise ValueError('a nonempty nonnegative integer spectrum is required')
    n=len(spectrum)-1
    denominator=lcm(*(comb(n,w) for w,a in enumerate(spectrum) if a))
    scaled=[a*(denominator//comb(n,w)) if a else 0 for w,a in enumerate(spectrum)]
    binomials=[[comb(a,v) for v in range(a+1)] for a in range(n+1)]
    result={}
    for u in range(n+1):
        for b in range(u+1):
            a=u-b
            count=sum(c*scaled[b+v]*scaled[u-v] for v,c in enumerate(binomials[a]))
            if count:result[u,b]=Q(comb(n,u)*comb(u,b)*count,denominator**2)
    assert sum(result.values(),Q(0))==sum(spectrum)**2
    return result


def weighted_shells(profiles,length,full_weight=Q(1)):
    full_weight=Q(full_weight)
    if full_weight<=0:raise ValueError('positive packet weight is required')
    result=[Q(0) for _ in range(length+1)]
    for (u,b),count in profiles.items():
        if not 0<=b<=u<=length or count<0:raise ValueError('invalid pair profile')
        result[u]+=count*full_weight**b
    return result


def conditioned_shells(profiles,length,double_probability):
    """Profile caps divided by a binomial double-packet conditioning mass.

    Given union size u, the reference distribution has b double packets
    with probability binomial(u,b) r^b (1-r)^(u-b). The returned shell cap
    sums each profile cap divided by that positive probability. It includes
    the zero pair; callers remove it before a nonzero-message union bound.
    This is not a claim that real pair profiles follow a binomial law.
    """
    r=Q(double_probability)
    if not 0<r<1:raise ValueError('an interior double-packet probability is required')
    result=[Q(0) for _ in range(length+1)]
    for (u,b),count in profiles.items():
        if not 0<=b<=u<=length or count<0:raise ValueError('invalid pair profile')
        result[u]+=count/(comb(u,b)*r**b*(1-r)**(u-b))
    return result


if __name__=='__main__':
    import model  # Establish the authenticated research imports.
    from bch_joint_support import authenticated_caps
    from support import union_shells
    caps=authenticated_caps();profile=pair_profiles(caps)
    assert weighted_shells(profile,256)==union_shells([caps,caps])
    print('Exact BCH pair-profile caps:',len(profile),'nonzero entries; all 257 union shells agree.')
