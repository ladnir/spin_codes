"""Reduce a positive shell majorant without increasing any component.

Floating LP solutions are proposals only. Exact convex repair against the
original valid majorant restores every shell inequality after rounding.
"""
from fractions import Fraction as Q
from math import comb,ceil,isfinite,log

import numpy as np
from scipy.optimize import linprog
from shared_mixture import verify


def limiting_shells(caps,mixture,count=12):
    """Smallest majorant/cap ratios, for diagnostics only (not proof inputs)."""
    if type(count) is not int or count<1:raise ValueError('positive diagnostic count required')
    verify(caps,mixture)
    n=len(caps)-1;rows=[]
    for u,cap in enumerate(caps):
        if not cap:continue
        upper=comb(n,u)*sum((Q(c)*Q(p)**u*(1-Q(p))**(n-u) for c,p in mixture),Q(0))
        rows.append((upper/cap,u))
    return [dict(support=u,ratio=float(ratio)) for ratio,u in sorted(rows)[:count]]


def pool(caps,mixtures):
    """Coordinatewise maximum of valid majorants on their union of centers."""
    result={}
    for mixture in mixtures:
        verify(caps,mixture)
        merged={}
        for c,p in mixture:merged[Q(p)]=merged.get(Q(p),Q(0))+Q(c)
        for p,c in merged.items():result[p]=max(result.get(p,Q(0)),c)
    output=[(c,p) for p,c in sorted(result.items())]
    verify(caps,output)
    return output


def repair(caps,mixture,multipliers,bits=48):
    verify(caps,mixture)
    if (len(multipliers)!=len(mixture) or any(not 0<=Q(x)<=1 for x in multipliers)
            or type(bits) is not int or not 16<=bits<=128):
        raise ValueError('one bounded rational multiplier per component required')
    n=len(caps)-1;delta=Q(0)
    for u,cap in enumerate(caps):
        if not cap:continue
        row=[Q(c)*comb(n,u)*Q(p)**u*(1-Q(p))**(n-u) for c,p in mixture]
        old=sum(row,Q(0));new=sum((v*Q(x) for v,x in zip(row,multipliers)),Q(0))
        if new<cap:
            delta=max(delta,(cap-new)/(old-new))
    # Round toward the original valid bound, never away from it.
    grid=1<<bits;delta=Q(ceil(delta*grid),grid)
    factors=[delta+(1-delta)*Q(x) for x in multipliers]
    result=[(Q(c)*x,Q(p)) for (c,p),x in zip(mixture,factors) if x]
    verify(caps,result)
    assert all(0<=x<=1 for x in factors)
    return result,delta


def union_gradient_logs(mixture,n,tilt):
    """Log gradient in component multipliers of the union's tilted mass.

    For P=sum_i c_i Bernoulli(p_i)^n, the nonzero-group comparison is
    2P+P union P. This is a proposal objective, not a distance bound.
    The common derivative factor two is omitted.
    """
    tilt=Q(tilt)
    if (type(n) is not int or n<1 or not 0<tilt<=1 or not mixture
            or any(Q(c)<=0 or not 0<Q(p)<=1 for c,p in mixture)):
        raise ValueError('positive mixture, dimension, and bounded union tilt required')
    rows=[(log(Q(c).numerator)-log(Q(c).denominator),Q(p)) for c,p in mixture]
    result=[]
    for lc,p in rows:
        terms=[n*log(float(1-p+p*tilt))]
        for ld,q in rows:
            activity=p+q-p*q
            terms.append(ld+n*log(float(1-activity+activity*tilt)))
        result.append(lc+float(np.logaddexp.reduce(terms)))
    return result


def prune(caps,mixture,cost_tilt=Q(1,4),seconds=10.,*,union_tilt=None):
    verify(caps,mixture)
    cost_tilt=Q(cost_tilt)
    if not 0<cost_tilt<=1 or not isfinite(seconds) or not 0<seconds<=60:
        raise ValueError('positive bounded cost tilt and LP time required')
    if not mixture or not any(caps):return [],dict(status=0,repair='0')
    n=len(caps)-1
    matrix=[]
    for u,cap in enumerate(caps):
        if cap:
            # Clipping only strengthens the covering constraint. The
            # all-ones proposal stays feasible because its row sum >=1.
            matrix.append([float(min(Q(10**6),Q(c)*comb(n,u)*Q(p)**u*
                (1-Q(p))**(n-u)/cap)) for c,p in mixture])
    logs=[log(c.numerator)-log(c.denominator)+n*log(float(1-p+p*cost_tilt))
          for c,p in ((Q(c),Q(p)) for c,p in mixture)]
    if union_tilt is not None:logs=union_gradient_logs(mixture,n,union_tilt)
    costs=np.maximum(np.exp(np.array(logs)-max(logs)),1e-16)
    fit=linprog(costs,A_ub=-np.array(matrix),b_ub=-np.ones(len(matrix)),
                bounds=(0,1),method='highs',options={'time_limit':seconds})
    info=dict(status=int(fit.status),message=fit.message)
    if union_tilt is not None:info['union_cost_tilt']=str(Q(union_tilt))
    if not fit.success or not np.isfinite(fit.x).all():
        return mixture,dict(info,repair='1')
    grid=1<<48
    proposed=[Q(min(grid,max(0,ceil(float(x)*grid))),grid) for x in fit.x]
    result,delta=repair(caps,mixture,proposed)
    return result,dict(info,repair=str(delta),limiting_shells=limiting_shells(caps,result))
