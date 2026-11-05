"""Bounded floating diagnostic for the p<=64 episode/route bound.

Not a certificate. Imports retained route moments without modifying them.
"""
from functools import lru_cache
import math
from pathlib import Path
import sys
import time
import numpy as np

theory = Path(__file__).resolve().parent.parent / 'theory'
if theory.is_dir():
    sys.path.insert(0, str(theory))
from stratified_fractional import binomial_logs, entropy_logchoose, LN2, NEG_INF

N, M, n = 1 << 21, 8192, 16384
MARKER_BITS = (.125, .25, .5, 1, 2, 4, 8, 16, 32, 64, 128)


def stable_excess(partials, mu, theta):
    # Same finite envelope as the retained helper, but stable for large markers.
    r = min(partials,256)
    best = NEG_INF
    for full in range(max(0,math.ceil(mu-r)),min(256-r,math.floor(mu))+1):
        probability = (mu-full)/r
        if full+r <= 128:
            continue
        if probability <= 0 or probability >= 1:
            d = max(0,full+(r if probability >= 1 else 0)-128)
            z = theta*d
            value = z+math.log(-math.expm1(-z)) if d else NEG_INF
        else:
            j = np.arange(max(0,129-full),r+1)
            z = theta*(full+j-128)
            logs = (binomial_logs(r)[j]+j*math.log(probability)
                    +(r-j)*math.log1p(-probability)+z+np.log(-np.expm1(-z)))
            peak = float(logs.max())
            value = peak+math.log(float(np.exp(logs-peak).sum()))
        best = max(best,value)
    return best


@lru_cache(maxsize=4096)
def positive_moment(v, partials, bits):
    theta = bits * LN2
    excess = stable_excess(partials, v/64, theta)
    if excess == NEG_INF:
        return NEG_INF
    log_full = (excess + math.log1p(math.exp(-excess)) if excess > 0
                else math.log1p(math.exp(excess)))
    return min(log_full, excess-math.log(-math.expm1(-theta)))


def point(p, v, c, threshold, terminal=True):
    resets = c-int(terminal)
    partials = min(2*c-int(terminal),256)
    if not (1 <= p <= 64 and 1 <= c and c+resets <= v <= n-int(not terminal)):
        return NEG_INF, None
    if (v+63)//64 + partials-1 <= 128:
        return NEG_INF, None
    count = (entropy_logchoose(M,p) + entropy_logchoose(v-resets-1,c-1)
             + entropy_logchoose(n-v+c-1,c-int(terminal)))
    deficit = max(0, 128*v-191*c-128*resets-threshold)
    witness = (deficit+1)//2
    cost = 64*resets + witness*107/128
    best = (float('inf'), None)
    for bits in MARKER_BITS:
        moment = positive_moment(v,partials,bits)
        value = count/LN2 + p*moment/LN2 - bits*(1+cost-1/64)
        if value < best[0]:
            best = value, bits
    return best


def main():
    started = time.monotonic()
    threshold = 491*N//1000
    results = []
    for c in (1,2,3,4,8,16,32,48,64,96,128,192,256):
        pivot = (threshold+319*c-128)//128
        candidates = {2*c-1,n,8192,pivot}
        for offset in (-256,-128,-64,-32,1,32,64,128,256):
            candidates.add(pivot+offset)
        for v in sorted(candidates):
            if not 2*c-1 <= v <= n:
                continue
            for p in (32,48,64):
                for terminal in (False,True):
                    value, marker = point(p,v,c,threshold,terminal)
                    results.append((value,p,v,c,terminal,marker))
            if time.monotonic()-started > 45:
                raise RuntimeError('Diagnostic time cap')
    for p in (32,48,64):
        selected = [row for row in results if row[1] == p]
        worst = max(selected)
        failures = sum(row[0] > -80 for row in selected)
        print(f"p={p}: {len(selected)} selected points; >-80 bits: {failures}; worst={worst}")
    print(f"Elapsed {time.monotonic()-started:.3f}s; selected points only, NO certificate")


if __name__ == '__main__':
    main()
