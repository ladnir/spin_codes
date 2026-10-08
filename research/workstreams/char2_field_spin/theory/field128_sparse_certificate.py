"""Directed active-row cover for K20 GRS[128,87], width16/state8.

The packet map, diagonal distribution, and rank cover are unchanged from
the certified GRS256 reference. Only the outer and its regional geometry
change. This replay covers every A in [1,1105], not the whole code alone.
Fixed rational marks are valid witnesses; no optimizer participates in
acceptance. Adaptive bisection checks a complete positive interval union.
"""
from fractions import Fraction as F
from hashlib import sha256
import argparse
import json

from flint import arb, ctx

from cauchy_firstmoment_screen import geometry
from joint_outer_screen import WITNESSES_128_87, packet
from sparse_active_interval_screen import anchor_zones, interpolated_anchors
from state8_sparse_certificate import (ball, integer_bits, polynomial_power,
                                      rational_marks, source_maximum,
                                      terminal_power)


def region(coefficients, active, rows, tau):
    """Positive binomial-ratio sum, including outward-rounded recurrence."""
    if not 0 <= active <= rows:
        raise ValueError('invalid regional occupancy')
    weights = [arb(1)]
    t = ball(tau)
    for h in range(active):
        weights.append(weights[-1]*t*(active-h)/(rows-h))
    return tuple(sum((p[h]*weights[h] for h in range(active+1)), arb(0))
                 for p in coefficients)


def complete_cover(entries, first, last):
    if not entries or entries[0]['lower'] != first or entries[-1]['upper'] != last:
        raise ValueError('wrong active-row endpoints')
    if any(e['lower'] > e['upper'] for e in entries) or any(
            a['upper']+1 != b['lower'] for a,b in zip(entries, entries[1:])):
        raise ValueError('active-row cover has a gap or overlap')


def certify_zone(lower, upper, anchor, point, g, cell_bits=64, precision=8192):
    """R_A <= R_upper entrywise and an exact source-factor maximum."""
    if not 1 <= lower <= upper <= g['rows'] or cell_bits < 40:
        raise ValueError('invalid replay zone')
    y,tau = rational_marks(anchor,point,g['n'])
    previous = ctx.prec
    try:
        ctx.prec = precision
        local = [packet(w,y,tau,joint=True) for w in range(17)]
        coefficients = polynomial_power(local,g['rows']//16,upper)
        output_factor = ball(y)**(-g['threshold'])
        target = arb(2)**(-cell_bits)
        cached = {}
        def bound(lo,hi):
            if hi not in cached:
                cached[hi] = terminal_power(region(coefficients,hi,g['rows'],tau),128)
            return (source_maximum(lo,hi,g['rows'],41,tau)*output_factor
                    *cached[hi]*(hi-lo+1))
        stack,entries,total = [(lower,upper)],[],arb(0)
        while stack:
            lo,hi = stack.pop()
            value = bound(lo,hi)
            endpoint = value.upper()
            if endpoint.is_finite() and 0 < endpoint <= target:
                entries.append(dict(lower=lo,upper=hi,anchor=anchor,
                                    certified_integer_bits=integer_bits(endpoint)))
                total += endpoint
            elif lo == hi:
                raise ArithmeticError(f'active={lo} not certified with anchor={anchor}: {value.str(30)}')
            else:
                mid = (lo+hi)//2
                stack.extend(((mid+1,hi),(lo,mid)))
        complete_cover(entries,lower,upper)
        return total.upper(),dict(anchor=anchor,y=str(y),tau=str(tau),
                                 precision_bits=precision,entries=entries,
                                 certified_integer_bits=integer_bits(total.upper()))
    finally:
        ctx.prec = previous


def certify(progress=None, last=1105, precision=8192, only_anchor=None):
    if last != 1105:
        raise ValueError('this checkpoint fixes active rows1..1105')
    g = geometry(20,128,87)
    # Extra fixed marks keep each replay zone close to its anchor. Reusing
    # the last tilt at larger anchors is a witness choice, not extrapolating
    # a previously measured margin: every interval is recomputed below.
    raw = dict(WITNESSES_128_87)
    raw.update({a: raw[512] for a in (768,1024,1105)})
    anchors = interpolated_anchors(interpolated_anchors(raw))
    zones = anchor_zones(1,last,anchors)
    if only_anchor is not None:
        zones = [z for z in zones if z[2] == only_anchor]
        if not zones:
            raise ValueError('anchor has no replay zone')
    previous = ctx.prec
    try:
        ctx.prec = 256
        entries,receipts,total = [],[],arb(0)
        for lo,hi,anchor in zones:
            upper,receipt = certify_zone(lo,hi,anchor,anchors[anchor],g,
                                         precision=precision)
            total += upper
            entries.extend(receipt['entries'])
            receipts.append(receipt)
            if progress:
                progress(receipt)
        if only_anchor is None:
            complete_cover(entries,1,last)
        endpoint = total.upper()
        bits = integer_bits(endpoint)
        if bits < 40:
            raise ArithmeticError('sparse union fails forty bits')
        q = endpoint.fmpq()
        digest = sha256(''.join(f"{e['lower']}:{e['upper']}:{e['anchor']}\n"
                               for e in entries).encode('ascii')).hexdigest()
        return dict(outer=[128,87],q=256,width=16,state=8,k_log=20,
                    requested_k=2**20,rows=g['rows'],k=g['k'],n=g['n'],
                    output_threshold=g['threshold'],active_range=[1,last]
                    if only_anchor is None else None,
                    complete_sparse_range=only_anchor is None,
                    interval_count=len(entries),partition_sha256=digest,
                    certified_integer_bits=bits,failure_upper_ball=str(endpoint),
                    failure_upper_dyadic=[str(q.p),str(q.q)],
                    rigorous=True,whole_code=False,
                    scope='specified active-row classes only',receipts=receipts)
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--precision',type=int,default=8192)
    p.add_argument('--anchor',type=int)
    p.add_argument('--quiet',action='store_true')
    args = p.parse_args()
    value = certify(None if args.quiet else lambda r: print(json.dumps(r),flush=True),
                    precision=args.precision,only_anchor=args.anchor)
    del value['receipts']
    print(json.dumps(value),flush=True)
