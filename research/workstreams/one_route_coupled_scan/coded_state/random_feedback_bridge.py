"""Uniform-feedback bridge: exact integer gate and bounded FLOAT diagnostics.

No encoder benchmark, generated files, or whole-code certificate.
"""
from fractions import Fraction
import math
import time
from feedback_support_probe import logchoose, route_log

N,K,B,S = 2**21,2**20,128,64
n = N//B
W = 491*N//1000


def output_exponent(v,c,j,b,s,w,beta):
    """q exponent of an output count with support factor 2^(b(v+j))."""
    assert c>=1 and v>=2*c-1 and j>=0 and 0<=beta<=1
    p = 1-beta
    reset_charge = min(s,1+2*(b-1)*p)
    return beta*w+b*p*(v+j)-reset_charge*(c-1)-s*j


def moment(count,tau):
    if count==0:
        return 0.0
    if count==N:
        return K*math.log2(tau)
    return route_log(count-128,math.log2(tau))


def screen_point(v,c):
    """No quiet nonzero batches. Both directions bound the SAME event."""
    assert 1<=c and 2*c-1<=v<=n
    # c intervals of total v batches touch at most this many regions.
    if (v+126*c)//64<=128:
        return ('impossible outer support',)
    patterns = 1+logchoose(v,c)+logchoose(n-v+c,c)
    candidates = []
    for tau in (Fraction(17,16),Fraction(9,8),Fraction(5,4),Fraction(3,2),2,4,8,16,64,256,4096):
        bits = math.log2(tau)
        fwd = patterns-S*(c-1)*bits+moment(B*v,tau)
        candidates.append((fwd,'forward',float(tau),None))
        for beta in (Fraction(0),Fraction(1,4),Fraction(1,2),Fraction(3,4),Fraction(7,8),Fraction(1)):
            e = output_exponent(v,c,0,B,S,W,beta)
            inv = patterns+(float(e)+v-B*(v-c))*bits+moment(B*(v-c),tau)
            candidates.append((inv,'inverse',float(tau),float(beta)))
    return min(candidates)


if __name__=='__main__':
    started = time.monotonic()
    bits = 128*(K-W)-N
    assert bits==318848
    print(f'EXACT: off-state coordinates in <=128 regions: failure <=2^-{bits}.',flush=True)
    for v,c in ((n,1),(n,n//2),(8192,1),(8192,32),(8192,64),(8192,128),
                (8192,256),(4096,128),(9000,64),(8000,64),(6000,64)):
        print(f'v={v}, c={c}, j=0: {screen_point(v,c)}',flush=True)
        if time.monotonic()-started>20:
            raise RuntimeError('20-second FLOAT diagnostic cap')
    print('Selected FLOAT upper bounds only; no complete cover.',flush=True)
