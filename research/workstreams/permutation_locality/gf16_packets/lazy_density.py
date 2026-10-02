"""Outward local envelopes preserving density through selected lazy updates.

These operators bound one conditional packet step. Complete outer-domain
coverage and aggregation remain obligations of the calling verifier.
"""
from fractions import Fraction as Q
from math import comb
from flint import arb
import occupancy_birth_classes as base


def density_caps(data,z,single=None):
    """Bound E[z^wt(A(b+CX)+X) 1{b+CX!=0} | J=j], for b!=0."""
    if not 0<z<=1:raise ValueError('positive output weight at most one required')
    W=data['windows']
    distance=min(sum(w*int(n) for w,n in enumerate(row)) for row in data['histograms'])
    weights=[Q(1)];caps=[]
    for j in range(W+1):
        bound=base.up(sum((base.aq(p)*z**max(0,distance-w) for w,p in enumerate(weights)),arb(0)))
        if single is not None and j:
            if j==1:exact=arb(single['density'])/single['denominator']
            else:
                factor=arb(W)/(W-j+1)*(((1+1/z)**4-1)/15)**(j-1)
                exact=factor*single['global_density']/single['denominator']
            bound=min(bound,base.up(exact))
        caps.append(min(arb(1),bound))
        next_weights=[Q(0)]*(len(weights)+4)
        for w,p in enumerate(weights):
            for b in range(1,5):next_weights[w+b]+=p*Q(comb(4,b),15)
        weights=next_weights
    return caps


def candidate(data,local,z,caps,through,quiet_classes=True):
    """Replace only the U-row lazy branch; retain its existing zero bound."""
    W=data['windows'];L=(1<<data['bits'])-1;alpha=arb(2)**-data['updates']
    if type(through) is not int or not 0<=through<=W or len(local)!=W+1 or len(caps)!=W+1:
        raise ValueError('matching geometry and valid occupancy cutoff required')
    result=[]
    for j,matrix in enumerate(local):
        value=matrix*arb(1)
        if j==0 and quiet_classes:
            value[2,1]=arb(0)
            for k,(level,n) in enumerate(zip(data['birth_class_levels'],data['birth_class_counts']),3):
                value[2,k]=base.up(value[2,k]+alpha*int(n)*z**int(level)/L)
        elif 1<=j<=through:
            value[2,1]=arb(0)
            value[2,2]=base.up(value[2,2]+alpha*caps[j])
        result.append(value)
    return result


def actual_caps(data,z,through):
    """Regenerate the checked one-packet census; never accept a saved cap."""
    from group_moment import maps
    import single_packet
    if type(through) is not int or not 0<=through<=data['windows']:
        raise ValueError('valid lazy-density occupancy cutoff required')
    images,columns,_=maps()
    if columns!=data['columns'] or len(images)!=(1<<data['bits']):
        raise ArithmeticError('lazy-density census map mismatch')
    single=single_packet.census(images,columns,z) if through else None
    return density_caps(data,z,single)


def refine_actual(data,local,z,through):
    return candidate(data,local,z,actual_caps(data,z,through),through)
