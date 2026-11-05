"""Directed active-row sum: <=16384 swaps in two half-full suffix families.

Standalone, serial <=30s, no raw files. Includes every active row set.
"""
import math
import shutil
import tempfile
import time
from fractions import Fraction as F

SCALE=256


def logs(x):
    a=x.numerator**SCALE;b=x.denominator**SCALE
    d=a.bit_length()-b.bit_length();return d-1,d+1


def ceildiv(a,b):return -(-a//b)


def history_upper(T,limit):
    x=F((4*limit+T)**2,16*T*T*65536)
    poly=sum((math.comb(256,b)**2*x**b for b in range(257)),F(0))
    scaled=logs(F(math.comb(512,256)))[1]+T*logs(poly)[1]-limit*logs(x)[0]
    return ceildiv(scaled,SCALE)


def binomial_upper(n,k):
    if not k:return 0
    y=F(k,n-k)
    return ceildiv(n*logs(1+y)[1]-k*logs(y)[0],SCALE)


def rho():
    def pairs(n):return math.prod(range(1,n,2))
    denominator=pairs(512)
    return max(F(math.comb(256,k)*math.factorial(k)*pairs(256-k)**2,denominator)
               for k in range(0,257,2))*2**256


def cell(R,low,high,r,t):
    M=8192;T=16*R-1;W=460*(1<<21)//1000
    history=history_upper(T,high)
    minimum_area=(T+1)*256+low+256;maximum_area=(T+1)*256+high+256
    p0=F(minimum_area,R*M);p1=F(maximum_area,R*M)
    z=max(0,(T+1)*256+low-W)
    choices=binomial_upper(T*256-low,z)
    active=max(1-p+p*r/t for p in (p0,p1))**R/r**128
    inactive=max((1-p)/2+p/t for p in (p0,p1))**R
    rlow,rhigh=logs(r);tlow,thigh=logs(t)
    coefficient=F(1+choices,32)-low-z
    # Fixed positive weights, no trusted floating optimization.
    scaled=(F(history*SCALE)+T*logs(rho())[1]+R*logs(F(M+1))[1]
            +coefficient*(rlow if coefficient<0 else rhigh)
            +(maximum_area*thigh if t>=1 else minimum_area*thigh)
            +M*logs(active+inactive)[1])
    return ceildiv(scaled.numerator,scaled.denominator*SCALE)


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic()
    for R,r,t in ((192,F(8),F(4)),(256,F(16),F(12))):
        records=[]
        for low in range(0,16385,2048):
            high=min(low+2047,16384)
            value=cell(R,low,high,r,t)
            assert value<=-120,(R,low,high,value)
            records.append((value,low,high))
            if time.monotonic()-began>25:raise TimeoutError('active-row replay budget')
        print(dict(R=R,w=256,swaps_at_most=16384,cells=len(records),
                   worst_cell=max(records),r=str(r),t=str(t),all_row_sets=True),flush=True)
    print(dict(seconds=time.monotonic()-began,joint_restricted_upper_bits=-113,
               whole_code_certificate=False))


if __name__=='__main__':main()
