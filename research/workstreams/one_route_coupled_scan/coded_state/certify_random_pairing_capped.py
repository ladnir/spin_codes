"""Directed capped-moment certificate for <=8192 support swaps.

Standalone standard library, serial <=30s, no files; supports stdin execution.
Only two constant-weight suffix families, not the whole code.
"""
import math
import time
import shutil
import tempfile
from fractions import Fraction as F

SCALE=256


def logs(value):
    a=value.numerator**SCALE;b=value.denominator**SCALE
    d=a.bit_length()-b.bit_length()
    return d-1,d+1


def ceildiv(a,b):return -(-a//b)


def comparison_density(n,w):
    def pairs(k):return math.prod(range(1,k,2))
    total=pairs(n)
    probabilities=[F(math.comb(w,k)*math.factorial(k)*pairs(w-k)*pairs(n-w-k),total)
                   for k in range(w%2,min(w,n-w)+1,2)]
    # These are probabilities of a SPECIFIED external endpoint set; the
    # C(n-w,k) factor in the crossing-count law has been cancelled.
    return max(probabilities)*2**(n-w)


def replay(R,r):
    M=8192;n=512;w=256;T=16*R-1;limit=8192;x=F(1,8192)
    polynomial=sum((math.comb(w,b)*math.comb(n-w,b)*x**b for b in range(w+1)),F(0))
    history_scaled=logs(F(math.comb(n,w)))[1]+T*logs(polynomial)[1]-limit*logs(x)[0]
    history=ceildiv(history_scaled,SCALE)
    W=460*(1<<21)//1000;z0=max(0,(T+1)*w-W)
    choices=(math.comb(T*w,z0)-1).bit_length()
    # Initial union enlarged to the whole batch; later swaps are not
    # discarded in the proof: their area increase cancels their death cost.
    base_area=(T+1)*w+n-w
    cover=(1+F(base_area,R*M))/2
    moment=sum((math.comb(R,c)*cover**c*(1-cover)**(R-c)
                *r**max(128-c,0) for c in range(R+1)),F(0))
    coefficient=F(1+choices,32)+base_area-128*M-z0
    lower,upper=logs(r)
    assert 1<r<(1<<32)
    rlog=lower if coefficient<0 else upper
    scaled=(F(history*SCALE)+T*logs(comparison_density(n,w))[1]
            +R*logs(F(M+1))[1]+coefficient*rlog+M*logs(moment)[1])
    result=ceildiv(scaled.numerator,scaled.denominator*SCALE)
    assert result<=-120
    return dict(R=R,w=w,swaps_at_most=limit,history_upper_bits=history,
                moment_tilt=str(r),directed_upper_bits=result,
                distance_threshold=W,whole_code_certificate=False)


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic()
    for args in ((192,F(3)),(256,F(8))):
        print(replay(*args),flush=True)
        if time.monotonic()-began>25:raise TimeoutError('capped replay budget')
    print(dict(seconds=time.monotonic()-began,joint_restricted_upper_bits=-119))


if __name__=='__main__':main()
