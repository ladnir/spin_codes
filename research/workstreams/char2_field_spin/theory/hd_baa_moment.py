"""Finite positive-transfer first-moment screen for striped MDS/acc/perm/acc.

First and final scans use independent uniform nonzero full-field feedback
multipliers. Source rows are coordinate-major before the first scan. There
is exactly one uniform global permutation between the scans. The row transfer
is exact; its positive-vector/boundary bound and numerical optimization are
floating screens, not directed certificates or a complete class cover.
"""
import argparse
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp

from gf8_exact_spectrum_dense import support_counts


@lru_cache(maxsize=8)
def row_log_transfer(b,k,q):
    """Log row transfer, divided by q^(k/2), before output-weight marking.

    Entry(r,s) sums over all MDS row words and outgoing binary states of
    weight s for any specified incoming state of weight r. The summand with
    a births, c returns, d live-to-live inputs has input weight a+c+d.
    """
    spectrum=support_counts(b,k,q)
    scale=k*math.log(q)/2
    loga=np.array([math.log(n)-scale if n else -math.inf for n in spectrum])
    logsurvive=math.log1p(-1/(q-1)) if q>2 else -math.inf
    # B[l,t] = sum_d C(l,d) survive^d A[t+d] / q^(k/2).
    weighted=np.full((b+1,b+1),-math.inf)
    weighted[0]=loga
    for size in range(1,b+1):
        count=b-size+1
        weighted[size,:count]=np.logaddexp(weighted[size-1,:count],
                                          logsurvive+weighted[size-1,1:count+1])
    combinations=[np.array([math.log(math.comb(n,j)) for j in range(n+1)])
                  for n in range(b+1)]
    out=np.full((b+1,b+1),-math.inf)
    for incoming in range(b+1):
        births=np.arange(b-incoming+1)
        for returns in range(incoming+1):
            following=incoming-returns+births
            terms=(combinations[b-incoming]+combinations[incoming][returns]
                   -returns*math.log(q-1)+weighted[incoming-returns,births+returns])
            out[incoming,following]=np.logaddexp(out[incoming,following],terms)
    return out


def first_log_moment(logy,b=128,k=64,q=65536,k_log=20,steps=48):
    """Bound log E[sum_message y^middle_weight], including the zero message.

    Every positive vector is valid; the iterations merely propose one.
    Boundary gluing costs at most 2^b max(v)/min(v), independent of R.
    """
    rows=2**k_log//k
    matrix=row_log_transfer(b,k,q)+np.arange(b+1)[None,:]*logy
    vector=np.zeros(b+1)
    for _ in range(steps):
        following=logsumexp(matrix+vector[None,:],axis=1)
        vector=following-float(following.max())
    loglambda=float((logsumexp(matrix+vector[None,:],axis=1)-vector).max())
    boundary=b*math.log(2)+float(vector.max()-vector.min())
    return rows*(loglambda+k*math.log(q)/2)+boundary


def final_log_moment(logx,logy,n,q):
    """Positive-vector bound for e0(T0(y)+x T1(y))^n 1."""
    x,y=math.exp(logx),math.exp(logy)
    off01=x*y
    off10=x/(q-1)
    diag11=y*(1+x*(q-2)/(q-1))
    eigen=(1+diag11+math.hypot(1-diag11,2*math.sqrt(off01*off10)))/2
    # off10/(lambda-diag11) avoids cancellation in (lambda-1)/off01.
    v=off10/(eigen-diag11)
    return n*math.log(eigen)+max(0.,-math.log(v))


def middle_class(alpha,b=128,k=64,q=65536,k_log=20,distance=.4):
    n=(2**k_log//k)*b
    h=round(alpha*n)
    e=math.ceil(distance*n)-1
    if not 0<h<=min(n,2*e):
        raise ValueError('middle class outside inverse-support range')
    first=lambda z:first_log_moment(z,b,k,q,k_log)-h*z
    optimum=minimize_scalar(first,bounds=(-math.log(q)-3,0),method='bounded',
                            options={'maxiter':50,'xatol':1e-6})
    logchoose=math.lgamma(n+1)-math.lgamma(h+1)-math.lgamma(n-h+1)
    second=lambda p:final_log_moment(p[0],p[1],n,q)-h*p[0]-e*p[1]-logchoose
    trials=[minimize(second,start,method='L-BFGS-B',bounds=[(-30,30),(-30,0)],
                     options={'maxiter':40,'ftol':1e-12})
            for start in ((0.,-1.),(math.log(q)/2,-math.log(q)/2))]
    best=min(trials,key=lambda r:r.fun)
    logbound=float(optimum.fun)+float(best.fun)-math.log(q-1)
    return dict(outer=[b,k],q=q,k_log=k_log,n=n,intermediate_weight=h,
                final_cutoff=e,margin_bits=-logbound/math.log(2),
                first_log_mark=float(optimum.x),final_log_marks=best.x.tolist(),
                first_component_bits=float(optimum.fun)/math.log(2),
                final_component_bits=float(best.fun)/math.log(2),
                rigorous=False,whole_code=False,
                scope='one intermediate-weight class, union over all messages; no first permutation')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--q',type=int,default=65536)
    p.add_argument('--outer',type=int,default=128)
    p.add_argument('--density',type=float,nargs='+',default=[.2,.4,.6,.75])
    a=p.parse_args()
    for alpha in a.density:
        print(json.dumps(middle_class(alpha,a.outer,a.outer//2,a.q)),flush=True)
