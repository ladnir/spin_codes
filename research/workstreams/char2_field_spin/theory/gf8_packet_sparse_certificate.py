"""Directed sparse replay for the repaired cubic width8/state4 candidate.

The default profile is GRS[128,83], requested K=2^20, one independent
regional route, and independent nonzero GF256 input/state diagonals.
Exact packet minima refer to the eight phases audited by gf8_packet_maps.
This certifies its explicit active-row interval only, never a whole code.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math

from flint import arb, arb_poly, ctx

from cauchy_firstmoment_screen import geometry
from field128_sparse_certificate import complete_cover, region
from gf8_packet_width_screen import local_matrix
from sparse_active_interval_screen import anchor_zones, interpolated_anchors
from state8_sparse_certificate import (ball, integer_bits, poly_product,
    rational_marks, source_maximum, terminal_power)


CONSTRUCTION = 'gf256-repaired-cubic-w8-s4-eta2-c3-regional-nonzero-diagonals'
MINIMA = (8, 6, 5, 4, 3, 2, 2, 1)
# Witness choices only: accepted bounds use their rational y,tau conversions.
ANCHORS = {
    1: (6.582496741565871, 5.541263545158426),
    4: (6.821146348026014, 5.354341280773075),
    8: (6.65, 3.05),
    16: (6.521550051392082, 2.449849392695374),
    32: (6.47, 2.23),
    64: (6.438181414836547, 2.164220696196884),
    96: (6.434007319113216, 2.146977127569882),
    128: (6.434830759728649, 2.144541008148496),
    160: (6.437102619529023, 2.1470413060449576),
    192: (6.439962833219732, 2.1517723678626974),
    256: (6.446410320362117, 2.164395241407062),
}


def metadata(k_log=20, dimension=83):
    if k_log not in (16, 18, 20) or not 1 <= dimension < 128:
        raise ValueError('invalid profile')
    g = geometry(k_log, 128, dimension)
    return dict(outer=[128, dimension], width=8, state=4, q=256,
                construction=CONSTRUCTION, k_log=k_log, requested_k=2**k_log,
                padded_k=g['k'], n=g['n'], rows=g['rows'],
                output_threshold=g['threshold'], eta=2, injection_constant=3,
                return_minima=list(MINIMA),
                rigorous=True, whole_code=False)


def polynomial_power(local, packets, degree):
    if len(local) != 9 or packets < 0 or degree < 0:
        raise ValueError('width-eight local matrices and nonnegative exponent/degree required')
    result = (arb_poly([1]), arb_poly([]), arb_poly([]), arb_poly([1]))
    factor = tuple(arb_poly([ball(local[w][i][j])*math.comb(8, w)
                            for w in range(min(8, degree)+1)])
                   for i, j in ((0, 0), (0, 1), (1, 0), (1, 1)))
    while packets:
        if packets & 1:
            result = poly_product(result, factor, degree)
        packets >>= 1
        if packets:
            factor = poly_product(factor, factor, degree)
    return result


def certify_zone(lower, upper, anchor, point, k_log=20, dimension=83,
                 cell_bits=64, precision=4096, split=True):
    g = metadata(k_log, dimension)
    if not 1 <= lower <= upper <= g['rows'] or cell_bits < 1:
        raise ValueError('invalid active interval or cell target')
    y, tau = rational_marks(anchor, point, g['n'])
    previous = ctx.prec
    try:
        ctx.prec = precision
        local = [local_matrix(w, y, tau, 8, 4, MINIMA) for w in range(9)]
        coefficients = polynomial_power(local, g['rows']//8, upper)
        output_factor = ball(y)**(-g['output_threshold'])
        target = arb(2)**(-cell_bits)
        cache = {}
        def bound(lo, hi):
            if hi not in cache:
                cache[hi] = terminal_power(region(coefficients, hi, g['rows'], tau), 128)
            return (source_maximum(lo, hi, g['rows'], 128-dimension, tau)
                    *output_factor*cache[hi]*(hi-lo+1))
        stack, entries, total = [(lower, upper)], [], arb(0)
        while stack:
            lo, hi = stack.pop()
            endpoint = bound(lo, hi).upper()
            if not endpoint.is_finite() or not endpoint > 0:
                raise ArithmeticError('increase Arb precision')
            if endpoint <= target or not split:
                entries.append(dict(lower=lo, upper=hi, anchor=anchor,
                                    certified_integer_bits=integer_bits(endpoint)))
                total += endpoint
            elif lo == hi:
                raise ArithmeticError(f'active={lo}: only {integer_bits(endpoint)} bits at anchor={anchor}')
            else:
                mid = (lo+hi)//2
                stack.extend(((mid+1, hi), (lo, mid)))
        complete_cover(entries, lower, upper)
        endpoint = total.upper()
        exact = endpoint.fmpq()
        return dict(anchor=anchor, y=str(y), tau=str(tau), precision_bits=precision,
                    entries=entries, certified_integer_bits=integer_bits(endpoint),
                    failure_upper_dyadic_hex=[hex(int(exact.p)), hex(int(exact.q))])
    finally:
        ctx.prec = previous


def certify(last=1057, k_log=20, dimension=83, cell_bits=64, precision=4096, progress=None):
    g = metadata(k_log, dimension)
    if not 1 <= last <= g['rows']:
        raise ValueError('active endpoint outside row range')
    raw = dict(ANCHORS)
    for anchor in (512, 768, 1024, last):
        if anchor > 256:
            raw[anchor] = ANCHORS[256]
    anchors = interpolated_anchors(interpolated_anchors(raw))
    entries, receipts, total = [], [], F(0)
    for lo, hi, anchor in anchor_zones(1, last, anchors):
        receipt = certify_zone(lo, hi, anchor, anchors[anchor], k_log, dimension,
                               cell_bits, precision)
        receipts.append(receipt)
        entries.extend(receipt['entries'])
        p, q = (int(value, 16) for value in receipt['failure_upper_dyadic_hex'])
        total += F(p, q)
        if progress:
            progress(dict(zone=[lo, hi], anchor=anchor,
                          certified_integer_bits=receipt['certified_integer_bits'],
                          intervals=len(receipt['entries'])))
    complete_cover(entries, 1, last)
    previous = ctx.prec
    try:
        ctx.prec = 256
        bits = integer_bits(ball(total).upper())
    finally:
        ctx.prec = previous
    digest = sha256(''.join(f"{e['lower']}:{e['upper']}:{e['anchor']}\n" for e in entries).encode('ascii')).hexdigest()
    return dict(**g, active_range=[1, last], complete=True,
                margin_bits_lower=bits, interval_count=len(entries),
                partition_sha256=digest, receipts=receipts,
                failure_upper_dyadic_hex=[hex(total.numerator), hex(total.denominator)],
                scope='complete stated active-row interval only')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--last', type=int, default=1057)
    parser.add_argument('--dimension', type=int, default=83)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--cell-bits', type=int, default=64)
    parser.add_argument('--precision', type=int, default=4096)
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args()
    result = certify(args.last, args.k_log, args.dimension, args.cell_bits, args.precision,
                     None if args.quiet else lambda value: print(json.dumps(value), flush=True))
    del result['receipts'], result['failure_upper_dyadic_hex']
    print(json.dumps(result), flush=True)
