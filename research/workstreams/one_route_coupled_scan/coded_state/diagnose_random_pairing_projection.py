"""Confined-support proof screen via the punctured RS enumerator.

FLOAT proposals only. Actual matching law and exact regional occupancy;
no independent-row model. Serial <=30s, no generated data files.
"""
import math
import argparse
import time
from functools import lru_cache
import numpy as np
from diagnose_random_pairing_confined import crossing_logs,choose2


def sum2(values):
    top=float(np.max(values))
    return top+math.log2(float(np.exp2(values-top).sum()))


def softplus2(x):
    return max(0.,x)+math.log1p(2**(-abs(x)))/math.log(2)


@lru_cache(maxsize=2048)
def route_proposal(R,w,n=512,M=8192):
    """Cap sum log2 of row enumerator envelopes, uniformly over masks.

    Negative association bounds the mgf by a product of row marginals.
    The route failure budget includes all C(n,w) masks and 2^17 profiles.
    """
    budget=80+17+choose2(n,w)
    if w==n:
        return dict(cap=M*softplus2(R-4096+31*R),budget=budget,lambda_=0)
    b=np.arange(R+1)
    law=np.array([choose2(R,int(j)) for j in b])+b*math.log2(w/n)+(R-b)*math.log2(1-w/n)
    envelope=np.array([softplus2(R-4096+31*int(j)) for j in b])
    best=None
    for lam in (1/1024,1/512,1/256,1/128,1/64,1/32,1/16,1/8,1/4,1/2,1.):
        cap=(M*sum2(law+lam*envelope)+budget)/lam
        if best is None or cap<best['cap']:best=dict(cap=cap,budget=budget,lambda_=lam)
    return best


def profile(R,w,threshold=None,n=512,M=8192):
    if threshold is None:threshold=460*(1<<21)//1000
    steps=R*M//n
    ks,logps=crossing_logs(n,w)
    # For each fixed external support, probability the crossing endpoints
    # equal it. Exact endpoint sets, not merely the number of crossings.
    terms=logps-np.array([choose2(n-w,int(k)) for k in ks])
    best=int(np.argmax(terms)); log_match=float(terms[best])
    zero_count=max(0,steps*w-threshold)
    if zero_count>(steps-1)*w:return dict(R=R,w=w,excluded=True)
    output=choose2((steps-1)*w,zero_count)-32*zero_count
    route=route_proposal(R,w,n,M)
    good=1+choose2(n,w)+(steps-1)*log_match+output+route['cap']
    return dict(R=R,w=w,good_route_upper_bits=round(good,3),
                route_log_enumerator_cap=round(route['cap'],3),
                match_bits_per_step=round(-log_match,6),maximizing_crossing=int(ks[best]),
                required_output_zeros=zero_count,lambda_=route['lambda_'],float_only=True)


def main():
    began=time.monotonic()
    parser=argparse.ArgumentParser();parser.add_argument('--neighbors',action='store_true')
    args=parser.parse_args()
    if args.neighbors:
        for R in (192,256):
            worst=None;failures=[]
            for w in range(1,513):
                result=profile(R,w)
                if worst is None or result['good_route_upper_bits']>worst['good_route_upper_bits']:worst=result
                if result['good_route_upper_bits']>-97:failures.append(w)
                if time.monotonic()-began>25:raise TimeoutError('neighbor screen budget')
            print(dict(R=R,profiles=512,failing_weights=failures,worst=worst),flush=True)
        print(dict(seconds=time.monotonic()-began,whole_code_certificate=False))
        return
    for R in (129,160,192,224,256):
        for w in (2,32,128,256,384,512):
            print(profile(R,w),flush=True)
            if time.monotonic()-began>25:raise TimeoutError('selected screen budget')
    print(dict(seconds=time.monotonic()-began,whole_code_certificate=False))


if __name__=='__main__':main()
