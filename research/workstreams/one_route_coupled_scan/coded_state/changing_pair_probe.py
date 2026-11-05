"""Bounded FLOAT screen for full-length paths with nonzero paired updates."""
import math
import time
from functools import lru_cache
from feedback_support_probe import N,K,M,W,MASK_N,logchoose,route_log

P=8191


def pair_rank_lower(s,h0,h1,u,v):
    assert 0<=h0<s and 0<=h1<s and 0<=u<=2*s and 0<=v<=2*s
    ru,rv = min(s,u),min(s,v)
    separate = max(0,h0+ru-s)+max(0,h1+rv-s)
    joint = max(0,min(s,h0+h1)+ru-s)+max(0,h0+h1+rv-2*s)
    return max(separate,joint)


@lru_cache(maxsize=64)
def log_poly(length,bits):
    values = [logchoose(length,j)+j*bits for j in range(64)]
    peak = max(values)
    return peak+math.log2(sum(2**(v-peak) for v in values))


def log_masks(d,h):
    limit = 63*(8192+d)
    if h==0:
        return 0.0
    if h==limit:
        return (P-d)*logchoose(256,63)+(2*d+1)*logchoose(128,63)
    return min(logchoose(MASK_N,h),min(
        (P-d)*log_poly(256,b)+(2*d+1)*log_poly(128,b)-h*b
        for b in (-16,-12,-8,-6,-4,-3,-2,-1,0,1,2,4,8,16)))


def point(d,h):
    U=128*(8193+d)
    v=min(W,U-h-(64*(U-W)+126)//127)
    output=v+(math.log2(U+1)+logchoose(U,max(v,U//2)))/128
    cost=logchoose(P,d)+log_masks(d,h)
    return min((cost+math.log2(tau)*(output-K)+route_log(h,math.log2(tau)),tau)
               for tau in (2,4,8,9,10,11,12,16,32,64,128,256,1024,4096,65536))


if __name__=='__main__':
    started=time.monotonic()
    for d in (0,1,4,8,16,32,64,128,256):
        hmax=63*(8192+d)
        selected=(0,19003,100000,300000,480000,hmax-8192,hmax-2048,hmax)
        worst=max((point(d,h)[0],h,point(d,h)[1]) for h in selected)
        print(f'd={d}, selected worst(log2,H,tau)={worst}',flush=True)
        if time.monotonic()-started>25:
            raise RuntimeError('25-second FLOAT cap')
    print('Selected FLOAT points only; no certificate or full coverage.')
