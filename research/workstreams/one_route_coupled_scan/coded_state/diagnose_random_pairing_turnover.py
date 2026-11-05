"""Changing masks at constant state weight: FLOAT aggregate screen.

Actual construction; no whole-code certificate. Serial <=30s, no files.
"""
import math
import time
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from diagnose_random_pairing_confined import choose2,crossing_logs
from diagnose_random_pairing_projection import sum2,softplus2


@lru_cache(maxsize=513)
def changes(w,n=512):
    b=np.arange(min(w,n-w)+1)
    return b,np.array([choose2(w,int(j))+choose2(n-w,int(j)) for j in b])


def history_count(T,w,B,n=512):
    """Coefficient bound for C(n,w) [x^B](sum C(w,b)C(n-w,b)x^b)^T."""
    if not B:return choose2(n,w)
    b,coeff=changes(w,n)
    if B==T*int(b[-1]):return choose2(n,w)+T*float(coeff[-1])
    def mean(logx):
        logs=coeff+b*logx;norm=sum2(logs)
        return float(np.dot(b,np.exp2(logs-norm)))-B/T
    logx=brentq(mean,-128,128)
    return choose2(n,w)+T*sum2(coeff+b*logx)-B*logx


def route_cap(R,area,history_bits,M=8192,conserved=False):
    """Convex-order Binomial envelope for unequal regional occupancies."""
    p=area/(R*M);budget=120+history_bits
    if p==1:return M*softplus2(32*(R-128)),0.
    b=np.arange(R+1)
    law=np.array([choose2(R,int(j)) for j in b])+b*math.log2(p)+(R-b)*math.log2(1-p)
    envelope=np.array([softplus2(R-4096+31*int(j)) for j in b])
    centered=b-R*p
    def objective(loglam):
        lam=2**loglam
        if conserved:
            # Independent Bernoulli entries conditioned on every column
            # size give the exact regional route. On that event sum b_i
            # is fixed. The tilt is NOT applied through negative association.
            def tilted_mean(a):
                logs=law+lam*(envelope-a*centered)
                return float(np.dot(centered,np.exp2(logs-sum2(logs))))
            a=brentq(tilted_mean,0,31,xtol=1e-10)
            moment=sum2(law+lam*(envelope-a*centered))
            conditioning=R*math.log2(M+1)
        else:moment=sum2(law+lam*envelope);conditioning=0
        return (M*moment+budget+conditioning)/lam
    result=minimize_scalar(objective,bounds=(-16,4),method='bounded',options=dict(xatol=1e-7))
    return float(result.fun),2**float(result.x)


def profile(R,w,B,delta=460,conserved=False):
    T=16*R-1;W=delta*(1<<21)//1000
    assert 0<=B<=T*min(w,512-w)
    area=(T+1)*w+B;intersection=T*w-B
    zeros=max(0,area-W)
    if zeros>intersection:return dict(R=R,w=w,B=B,excluded=True)
    history=history_count(T,w,B)
    ks,law=crossing_logs(512,w)
    logm=max(float(p)-choose2(512-w,int(k)) for k,p in zip(ks,law))
    cap,lam=route_cap(R,area,history,conserved=conserved)
    upper=(1+history+T*logm+B-32*(B+zeros)
           +choose2(intersection,zeros)+cap)
    return dict(R=R,w=w,B=B,history_bits=round(history,2),route_cap=round(cap,2),
                lambda_=round(lam,7),good_route_upper_bits=round(upper,2),
                conserved=conserved,float_only=True)


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--conserved',action='store_true')
    args=parser.parse_args()
    began=time.monotonic()
    for R in (192,256):
        for B in (0,1,16,256,4096,32768,131072,262144,450000):
            print(profile(R,256,B,conserved=args.conserved),flush=True)
            if time.monotonic()-began>25:raise TimeoutError('turnover screen budget')
    print(dict(seconds=time.monotonic()-began,whole_code_certificate=False))


if __name__=='__main__':main()
