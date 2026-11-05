"""Bounded FLOAT screen of a larger ideal-feedback bridge. No certificate."""
from fractions import Fraction
from functools import lru_cache
import math
import time
from random_feedback_bridge import N,K,W,output_exponent
from feedback_support_probe import logchoose


@lru_cache(maxsize=1024)
def interval_moment(count,c,tau):
    """At most 2c partially occupied regions; maximize their full-region count."""
    r = min(256,2*c)
    mean = count/8192
    bits = math.log2(tau)
    if count==0:
        return 0.0
    if count==N:
        return K*bits
    scores = []
    for full in range(max(0,math.ceil(mean-r)),min(256-r,math.floor(mean))+1):
        p = (mean-full)/r
        if p==0 or p==1:
            scores.append(max(0,full+int(p)*r-128)*bits)
            continue
        terms = [logchoose(r,j)+j*math.log2(p)+(r-j)*math.log2(1-p)
                 +max(0,full+j-128)*bits for j in range(r+1)]
        peak = max(terms)
        scores.append(peak+math.log2(sum(2**(x-peak) for x in terms)))
    return 8192*max(scores)


def point(v,c,b=256,s=256):
    n = N//b
    region_batches = 8192//b
    if (v+(2*region_batches-2)*c)//region_batches<=128:
        return ('impossible outer support',)
    placements = 1+logchoose(v,c)+logchoose(n-v+c,c)
    out = []
    for tau in (Fraction(17,16),Fraction(9,8),Fraction(5,4),Fraction(3,2),2,4,8,16,64,256):
        bits = math.log2(tau)
        out.append((placements-s*(c-1)*bits+interval_moment(b*v,c,tau),'forward',float(tau),None))
        route = interval_moment(b*(v-c),c,tau)
        for beta in (Fraction(0),Fraction(1,4),Fraction(1,2),Fraction(3,4),Fraction(7,8),Fraction(1)):
            e = output_exponent(v,c,0,b,s,W,beta)
            out.append((placements+(float(e)+b*v/128-b*(v-c))*bits+route,'inverse',float(tau),float(beta)))
    return min(out)


if __name__=='__main__':
    started = time.monotonic()
    for v,c in ((8192,1),(8192,4096),(4096,1),(4096,32),(4096,64),(4096,128),
                (4096,256),(4022,128),(4000,64),(4500,64),(3000,64),(2048,128)):
        print(f'b=s=256: v={v}, c={c}: {point(v,c)}',flush=True)
        if time.monotonic()-started>40:
            raise RuntimeError('40-second screen cap')
    print(f'Elapsed {time.monotonic()-started:.3f}s; selected FLOAT points, no certificate.')
