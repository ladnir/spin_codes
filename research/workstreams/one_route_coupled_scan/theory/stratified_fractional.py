"""Bounded floating diagnostics for a stratified MDS outer and uniform inner.

This is not a directed certificate. No output files or unbounded caches.
The root agent runs all numerical work serially.
"""

import argparse
import math
import time
from functools import lru_cache

import numpy as np

LN2 = math.log(2.0)
NEG_INF = float('-inf')


def logchoose(n, k):
    if k < 0 or k > n:
        return NEG_INF
    k = min(k, n-k)
    if k == 0:
        return 0.0
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)


def logsum(values):
    peak = max(values, default=NEG_INF)
    if peak == NEG_INF:
        return peak
    return peak + math.log(sum(math.exp(x-peak) for x in values))


@lru_cache(maxsize=257)
def binomial_logs(r):
    return np.array([logchoose(r, j) for j in range(r+1)])


def excess_log(B, k, c, mu, theta):
    """Upper log E[exp(theta*(T-(B-k))_+)-1], retaining <=2c partial regions."""
    r = min(2*c, B)
    best = NEG_INF
    for f in range(max(0, math.ceil(mu-r)), min(B-r, math.floor(mu))+1):
        p = (mu-f)/r
        if f+r <= B-k:
            continue
        if p <= 0 or p >= 1:
            d = max(0, f+(r if p >= 1 else 0)-(B-k))
            val = math.log(math.expm1(theta*d)) if d else NEG_INF
        else:
            j = np.arange(max(0, B-k-f+1), r+1)
            if not len(j):
                continue
            z = theta*(f+j-(B-k))
            logs = (binomial_logs(r)[j]+j*math.log(p)
                    +(r-j)*math.log1p(-p)+z+np.log(-np.expm1(-z)))
            peak = float(logs.max())
            val = peak+math.log(float(np.exp(logs-peak).sum()))
        best = max(best, val)
    return best


def route_log(M, excess, theta):
    if excess == NEG_INF:
        return NEG_INF
    if excess < -700:
        return math.log(M)+excess-math.log(-math.expm1(-theta))
    logH = math.log1p(math.exp(excess)) if excess < 700 else excess
    z = M*logH
    corrected = (math.log(math.expm1(z)) if z < 700 else z)
    corrected -= math.log(-math.expm1(-theta))
    return min(z, corrected)


def representation_log(N, L, u, c):
    free = u-L*(c-1)
    if free < c or u > N or u < 1:
        return NEG_INF
    return logchoose(free-1, c-1)+logchoose(N-u+c, c)


def output_dimension(free, c, e, field_bits, simple=False):
    """Upper log_q count of outputs, with c prescribed nonzero activations."""
    if c > e or free < c:
        return NEG_INF
    if e >= free:
        return float(free)
    if simple:
        return min(float(free), e+free/field_bits)
    # Replace q-1 by q. The geometric tail contributes <2N/q in natural logs.
    return e+(logchoose(free-c, e-c)+2*free/(2.0**field_bits))/(field_bits*LN2)


def point(N, B, L, delta, u, c, field_bits=128, simple=False):
    e = math.floor(delta*N)
    count = representation_log(N, L, u, c)
    if count == NEG_INF:
        return NEG_INF, 0.0
    # Every block contributes at most two boundary coordinates per interval.
    if u*B/N+2*c <= B//2:
        return NEG_INF, 0.0
    free = u-L*(c-1)
    w = output_dimension(free, c, e, field_bits, simple)
    if w == NEG_INF:
        return NEG_INF, 0.0
    best = (float('inf'), 0.0)
    for tau_bits in (.125, .25, .5, 1, 2, 3, 4, 6, 8):
        theta = tau_bits*LN2
        # Fixed cardinality removes the linear fluctuation for large U:
        # D(U)=u-N/2+D(U complement), at rate one half.
        shifted = max(0, u-N//2)
        mu = min(u, N-u)*B/N
        excess = excess_log(B, B//2, c, mu, theta)
        if shifted:
            # D(U)>0 is already guaranteed; do not impose D(complement)>0.
            logH = 0.0 if excess == NEG_INF else np.logaddexp(0.0, excess)
            routed = (N//B)*float(logH)
        else:
            routed = route_log(N//B, excess, theta)
        value = count+theta*(c-u+w+shifted)+routed
        if value < best[0]:
            best = value, tau_bits
    return best[0]/LN2, best[1]


def entropy_logchoose(n, k):
    if k < 0 or k > n:
        return NEG_INF
    if k == 0 or k == n:
        return 0.0
    return n*math.log(n)-k*math.log(k)-(n-k)*math.log(n-k)


def box(N, B, L, e, ul, uh, cl, ch):
    """Floating proposal for a uniform upper bound over all feasible points."""
    ul = max(ul, L*(cl-1)+cl)
    ch = min(ch, (uh+L)//(L+1))
    if ul > uh or cl > ch or uh*B/N+2*ch <= B//2:
        return NEG_INF, (5,4)
    free = uh-L*(cl-1)
    n1 = free-1
    n2 = N-ul+ch
    count = entropy_logchoose(n1, min(ch-1, n1//2))
    count += entropy_logchoose(n2, min(ch, n2//2))
    w = min(float(free), e+free/128)
    coefficient = ch-min(ul, N//2)+w
    closest = min(uh, N//2) if ul <= N//2 else N-ul
    mu = closest*B/N
    best = (float('inf'), (5,4))
    for num, den in ((17,16), (9,8), (5,4), (3,2), (2,1), (4,1), (8,1)):
        theta = math.log(num/den)
        excess = excess_log(B, B//2, ch, mu, theta)
        if uh > N//2:
            routed = (N//B)*float(np.logaddexp(0.0, excess))
        else:
            routed = route_log(N//B, excess, theta)
        value = (count+theta*coefficient+routed)/LN2
        if value < best[0]:
            best = value, (num, den)
    return best


def propose_cover(B, L, delta, threshold=-100, seconds=120):
    N = 2**21
    e = math.floor(delta*N)
    queue = [(1, N, 1, (N+L)//(L+1))]
    passed = []
    started = time.monotonic()
    checked = 0
    worst = (NEG_INF, None)
    while queue:
        if time.monotonic()-started > seconds or checked > 50000:
            raise RuntimeError(f'Cover proposal cap reached: {checked} checks, {len(queue)} pending')
        ul, uh, cl, ch = queue.pop()
        ul = max(ul, L*(cl-1)+cl)
        ch = min(ch, (uh+L)//(L+1))
        if ul > uh or cl > ch:
            continue
        checked += 1
        value, tau = box(N, B, L, e, ul, uh, cl, ch)
        if value <= threshold:
            passed.append((ul, uh, cl, ch, tau, value))
            if value > worst[0]:
                worst = value, (ul, uh, cl, ch, tau)
            continue
        if ul == uh and cl == ch:
            raise RuntimeError(f'Point does not close: {(ul,cl,value,tau)}')
        if uh-ul >= L*(ch-cl) and ul < uh:
            middle = (ul+uh)//2
            queue.extend(((ul,middle,cl,ch),(middle+1,uh,cl,ch)))
        else:
            middle = (cl+ch)//2
            queue.extend(((ul,uh,cl,middle),(ul,uh,middle+1,ch)))
    print(dict(kind='floating_complete_cover_not_certificate', B=B,L=L,delta=delta,
               boxes=len(passed),checks=checked,worst=worst,
               seconds=round(time.monotonic()-started,3)), flush=True)
    return passed


def screen(B, L, delta, simple=False):
    N = 2**21
    cs = set(range(1, min(B//2+2, 72)))
    cs.update((96, 128, 192, 256, 384, 512, 768, 1024, 2048,
               4096, 8192, 16384, (N+L)//(L+1)))
    aa = (.01, .05, .1, .2, .3, .34, .35, .36, .38, .4,
          .42, .44, .45, .46, .47, .48, .49, .5, .52, .55,
          .6, .7, .8, .9, 1.0)
    worst = (NEG_INF, None)
    checked = 0
    started = time.monotonic()
    for c in sorted(cs):
        if time.monotonic()-started > 120:
            raise RuntimeError('Diagnostic exceeded its 120-second cap')
        minimum = L*(c-1)+c
        if minimum > N:
            continue
        us = {max(minimum, round(a*N)) for a in aa}
        us.update((minimum, min(N, minimum+round(.01*N)),
                   min(N, math.floor(delta*N)+L*(c-1))))
        for u in sorted(us):
            value, marker = point(N, B, L, delta, u, c, simple=simple)
            checked += 1
            if value > worst[0]:
                worst = value, dict(u=u, c=c, tau_bits=marker)
    print(dict(kind='floating_grid_not_certificate', B=B, L=L, delta=delta,
               points=checked, worst_log2=worst[0], at=worst[1],
               seconds=round(time.monotonic()-started, 3)), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--B', type=int, default=256, choices=(64,128,256))
    ap.add_argument('--L', type=int, default=256)
    ap.add_argument('--delta', type=float, default=.45)
    ap.add_argument('--simple-volume', action='store_true')
    ap.add_argument('--cover', action='store_true')
    args = ap.parse_args()
    if args.cover:
        propose_cover(args.B, args.L, args.delta)
    else:
        screen(args.B, args.L, args.delta, args.simple_volume)
