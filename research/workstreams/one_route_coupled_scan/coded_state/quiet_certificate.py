"""All-quiet branch for B_t=B*diag(nonzero GF256 coefficients).

This is a distinct candidate, NOT a certificate for fixed B or mixed trajectories.
See QUIET_BRANCH.md. Bounded, in-memory floating proposal and directed replay.
"""
import argparse
from decimal import Decimal
from fractions import Fraction
import math
from pathlib import Path
import sys
import time

theory = Path(__file__).resolve().parent.parent / 'theory'
if theory.is_dir():
    sys.path.insert(0,str(theory))
from stratified_fractional import entropy_logchoose, route_log, LN2
from stratified_certificate import exact_excess, Upward
from route_episode_probe import stable_excess

N, M, n = 1 << 21,8192,16384
KAPPA = Fraction(1279,20)
MARKERS = ((17,16),(9,8),(5,4),(3,2),(2,1),(4,1),(8,1),(16,1),(256,1))


def bounds(lo,hi):
    count_at = min(hi,(n-1)//2) if lo<=(n-1)//2 else lo
    closest = hi+1 if hi+1<=n//2 else n-lo-1 if lo+1>=n//2 else n//2
    def coefficient(j):
        return max(0,128*(j+1)-N//2)-KAPPA*j
    return count_at,closest,max(coefficient(lo),coefficient(hi)),hi+1<=n//2


def float_box(lo,hi):
    count_at,closest,coefficient,positive = bounds(lo,hi)
    count = entropy_logchoose(n-1,count_at)
    best = float('inf'),None
    for marker in MARKERS:
        theta = math.log(marker[0]/marker[1])
        excess = stable_excess(256,closest/64,theta)
        if positive:
            routed = route_log(M,excess,theta)
        else:
            log_h = excess+math.log1p(math.exp(-excess)) if excess>0 else math.log1p(math.exp(excess))
            routed = M*log_h
        value = (count+theta*float(coefficient)+routed)/LN2
        if value<best[0]:
            best = value,marker
    return best


def cover(weight,seconds=30,term_bits=80):
    # 255^-64 <=256^-KAPPA, proved by an integer comparison.
    assert 255**160 >= 2**1279
    maximum = weight//65
    started = time.monotonic()
    stack,accepted = [(128,maximum)],[]
    checked = 0
    while stack:
        if time.monotonic()-started>seconds or checked>40000:
            raise RuntimeError('Quiet cover cap')
        lo,hi = stack.pop()
        value,marker = float_box(lo,hi)
        checked += 1
        if value<=-term_bits:
            accepted.append((lo,hi,marker,value))
        elif lo==hi:
            raise RuntimeError(f'Quiet cell does not close: {lo,value,marker}')
        else:
            middle = (lo+hi)//2
            stack.extend(((middle+1,hi),(lo,middle)))
    ordered = sorted(accepted)
    assert ordered[0][0] ==128 and ordered[-1][1] ==maximum
    assert all(a[1]+1==b[0] for a,b in zip(ordered,ordered[1:]))
    print(f'Quiet floating full cover: {len(accepted)} boxes, {checked} checks, '
          f'j=128..{maximum} plus a free terminal batch, worst={max(accepted,key=lambda x:x[3])}, '
          f'seconds={time.monotonic()-started:.3f}',flush=True)
    return accepted,maximum-127


def directed(boxes,precision,seconds=60):
    started = time.monotonic()
    up = Upward(precision)
    worst = Decimal('-Infinity')
    for lo,hi,marker,_ in boxes:
        if time.monotonic()-started>seconds:
            raise RuntimeError('Quiet directed cap')
        count_at,closest,coefficient,positive = bounds(lo,hi)
        excess = exact_excess(256,128,128*closest,N,*marker)
        value = up.entropy_choose(n-1,count_at)
        value = up.up.add(value,up.times_log(coefficient,Fraction(*marker)))
        value = up.up.add(value,up.routed(M,excess,Fraction(*marker),positive))
        log2 = up.up.divide(value,up.log_int(2)[1 if value<0 else 0])
        assert log2<=-80,(lo,hi,marker,log2)
        worst = max(worst,log2)
    print(f'Quiet directed {precision}: PASS; worst upper log2={worst}; '
          f'seconds={time.monotonic()-started:.3f}',flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--float-only',action='store_true')
    args = parser.parse_args()
    boxes,cells = cover(491*N//1000)
    if not args.float_only:
        for precision in (80,120):
            directed(boxes,precision)
        margin = 80-(cells-1).bit_length()
        print(f'All-quiet scaled-B branch: >49.1% except with probability <=2^-{margin}. '
              'Mixed active/quiet trajectories remain open.')
