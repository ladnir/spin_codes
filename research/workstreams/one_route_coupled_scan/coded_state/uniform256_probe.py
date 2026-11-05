"""Bounded FLOAT diagnostics for the 256-state running-sum bridge.

No certificate, encoder timings, output files, or changes to frozen helpers.
"""
import math
import time
import sys
from pathlib import Path
from fractions import Fraction
from feedback_support_probe import logchoose
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'theory'))
from stratified_fractional import excess_log

N, B, M, W = 2**21, 256, 8192, 491*2**21//1000
TAUS = tuple(Fraction(32+i,32) for i in range(1,33))+(2,4,8)


def softplus2(x):
    return max(x,0)+math.log2(1+2**(-abs(x)))


def output_log(v,c,b=B):
    """Chernoff volume with one free episode and c-1 parity checks.

    P_L(z) <= P_2(z) (1+(q-1)z)^(L-2). Replace q-1 by q.
    Convex minimization only proposes a witness; this is FLOAT.
    """
    pairs = b*(c-1)
    singles = b*(v-2*(c-1))
    assert singles>=b
    def value(beta):
        return 128*beta*W+pairs*softplus2(128*(1-2*beta))+singles*softplus2(128*(1-beta))
    lo,hi=0.0,1.1
    for _ in range(55):
        a=(2*lo+hi)/3
        b=(lo+2*hi)/3
        if value(a)<value(b): hi=b
        else: lo=a
    beta=(lo+hi)/2
    return min((value(0),0),(value(beta),beta))


def point(v,c,b=B,outer=256):
    n=N//b
    region=N//outer
    g=region//b
    if (v+(2*g-2)*c)//g<=outer//2:
        return (-math.inf,'impossible',0,0)
    placements=1+logchoose(v,c)+logchoose(n-v+c,c)
    volume,beta=output_log(v,c,b)
    out=[]
    for tau in TAUS:
        bits=math.log2(tau)
        def moment(size):
            log_excess=excess_log(outer,outer//2,c,min(size,N-size)/region,math.log(tau))
            log_h=max(log_excess,0)+math.log1p(math.exp(-abs(log_excess)))
            return region*log_h/math.log(2)+max(0,size-N//2)*bits
        a=b*(v-c)
        # Complementation avoids the spurious linear occupancy variance.
        route=moment(a)
        out.append((placements+(volume/128-a)*bits+route,'inverse',float(tau),beta))
        u=b*v
        route=moment(u)
        out.append((placements-b*(c-1)*bits+route,'forward',float(tau),0))
    return min(out)


def structural_screen():
    """Exact expectation caps for all interval counts, not a distance proof.

    Convex equalization with r<=2c partial regions maximizes the one-block
    expected excess at fixed total size. Allow the final reset to be free.
    Numerators have denominator 32*r. Only c<=128 need explicit checking.
    """
    worst=(Fraction(0),None)
    started=time.monotonic()
    for c in range(1,129):
        r=2*c
        best=Fraction(0)
        for full in range(max(0,126-r),min(256-r,125)+1):
            p=4022-32*full
            d=32*r
            if not 0<=p<=d: continue
            num=sum((full+j-128)*math.comb(r,j)*p**j*(d-p)**(r-j)
                    for j in range(max(0,129-full),r+1))
            mean=Fraction(M*num,d**r)
            best=max(best,mean)
        if c==1:
            assert best==0
        else:
            ratio=best/(c-1)
            if ratio>worst[0]: worst=ratio,c
            assert ratio<256
        assert time.monotonic()-started<40, 'structural time cap'
    # For c>=129 the all-partial maximum and increasing reset budget suffice.
    assert best<256*128
    print('EXACT expectation-only reset screen: max E[D]/(c-1)=',
          float(worst[0]),'at c=',worst[1],'; below256. NOT a route-tail certificate.',flush=True)


if __name__=='__main__':
    structural_screen()
    started=time.monotonic()
    for v,c in ((8192,1),(8192,4096),(4096,1),(4096,32),(4096,64),
                (4096,128),(4096,256),(4022,128),(4000,64),(4500,64),(3000,64),(2048,128)):
        print(v,c,point(v,c),flush=True)
        assert time.monotonic()-started<40,'FLOAT screen time cap'
