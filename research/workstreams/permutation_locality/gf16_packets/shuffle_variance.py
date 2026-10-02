"""Local shuffle/iid density bound retaining Bernoulli-count variance.

The scalar cover can supply a mean interval and an exactly checked
affine lower bound on the composition's variance.
See SHUFFLE_VARIANCE.md for the finite-binomial ratio argument.
"""
from fractions import Fraction as Q
from math import comb
from flint import arb


def aq(value):
    value=Q(value)
    return arb(value.numerator)/value.denominator


def upper(value):
    m,e=value.upper().man_exp()
    return arb(int(m))*arb(2)**int(e)


def density_upper(slots,mean,variance_lower):
    mean=Q(mean);variance_lower=Q(variance_lower)
    if (type(slots) is not int or slots<1 or not 0<=mean<=slots
            or not 0<=variance_lower<=mean*(1-mean/slots)):
        raise ValueError('valid slot count, mean, and variance lower bound required')
    if mean in (0,slots):return arb(1)
    # Fourier inversion: every atom is <= exp(-V) I_0(V).
    # This expression decreases with V, so a variance lower bound suffices.
    peak=arb(1) if not variance_lower else (-aq(variance_lower)).exp()*aq(variance_lower).bessel_i(0)
    p=aq(mean/slots);lo=mean.numerator//mean.denominator
    indices={lo,-(-mean.numerator//mean.denominator)}
    return max(upper(peak/(comb(slots,k)*p**k*(1-p)**(slots-k))) for k in indices)


def density_interval_upper(slots,mean_interval,variance_lower):
    """One density factor valid for every mean in a closed interior interval."""
    lo,hi=map(Q,mean_interval);variance_lower=Q(variance_lower)
    if (type(slots) is not int or slots<1 or not 0<lo<=hi<slots
            or not 0<=variance_lower<=Q(slots,4)):
        raise ValueError('interior mean interval and nonnegative variance lower bound required')
    peak=arb(1) if not variance_lower else (-aq(variance_lower)).exp()*aq(variance_lower).bessel_i(0)
    result=arb(0)
    for k in range(lo.numerator//lo.denominator,-(-hi.numerator//hi.denominator)+1):
        # For fixed k, k=floor(mu) or ceil(mu) requires |mu-k|<=1.
        # The binomial mass is unimodal in mu, so its minimum over this
        # (slightly enlarged, closed) interval is attained at an endpoint.
        left=max(lo,Q(k-1));right=min(hi,Q(k+1))
        if left>right:continue
        for mean in {left,right}:
            p=aq(mean/slots)
            result=max(result,upper(peak/(comb(slots,k)*p**k*(1-p)**(slots-k))))
    return result


def variance_dual(features,active,cell,active_min,dual=None):
    """Exact affine minorant; floating LP only chooses its two slopes."""
    features=list(map(Q,features));active=list(map(Q,active));lo,hi=map(Q,cell);active_min=Q(active_min)
    if (not features or len(features)!=len(active) or any(not 0<=f<=1 for f in features)
            or any(a not in (0,1) for a in active) or not 0<=lo<=hi<=1 or not 0<=active_min<=1):
        raise ValueError('valid component features and composition constraints required')
    values=[f*(1-f) for f in features]
    if dual is None:
        import numpy as np
        from scipy.optimize import linprog
        fit=linprog(np.array(list(map(float,values))),
            A_ub=np.array([list(map(float,features)),list(map(float,(-f for f in features))),list(map(float,(-a for a in active)))]),
            b_ub=[float(hi),-float(lo),-float(active_min)],A_eq=np.ones((1,len(features))),b_eq=[1],bounds=(0,None),method='highs')
        slopes=(fit.ineqlin.marginals[0]-fit.ineqlin.marginals[1],max(0.,-fit.ineqlin.marginals[2])) if fit.success else (0.,0.)
        dual=[Q(round(float(x)*10**10),10**10) for x in slopes]
    if len(dual)!=2:raise ValueError('two variance dual slopes required')
    b,c=map(Q,dual)
    if c<0:raise ValueError('nonnegative active-count dual required')
    a=min(v-b*f-c*t for v,f,t in zip(values,features,active))
    lower=max(Q(0),a+min(b*lo,b*hi)+c*active_min)
    return lower,(b,c)
