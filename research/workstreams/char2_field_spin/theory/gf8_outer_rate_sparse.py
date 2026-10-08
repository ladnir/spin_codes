"""Directed sparse-class replay for variable-dimension length128 GRS outers.

The fixed width16/state8 packet law and independent regional routing are
unchanged. Every accepted bound uses rational marks and outward Arb balls.
Point mode certifies only the selected active-row class; cover mode checks
every active count in its explicit range. Neither mode is a whole-code proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math

from flint import arb, arb_mat, arb_poly, ctx

from cauchy_firstmoment_screen import geometry
from field128_sparse_certificate import complete_cover, region
from joint_outer_screen import WITNESSES_128_87, packet
from sparse_active_interval_screen import anchor_zones, interpolated_anchors
from state8_sparse_certificate import (ball, integer_bits, polynomial_power,
    rational_marks, source_maximum, terminal_power)


def rate_geometry(k_log, dimension):
    if not 1 <= dimension < 128 or k_log < 1:
        raise ValueError('positive size and dimension below128 required')
    return geometry(k_log, 128, dimension)


def metadata(k_log, dimension):
    g = rate_geometry(k_log, dimension)
    return dict(outer=[128, dimension], q=256, width=16, state=8,
        construction='symmetric-cauchy-w16-s8-eta2-regional-nonzero-diagonals',
        k_log=k_log, requested_k=2**k_log, padded_k=g['k'], rows=g['rows'],
        n=g['n'], output_threshold=g['threshold'], rigorous=True,
        whole_code=False, scope='explicit active-row classes only')


def witnesses(last):
    raw = dict(WITNESSES_128_87)
    # These are valid mark choices, not extrapolated probability bounds.
    for anchor in (768, 1024, 1536, 2048, last):
        if anchor > 512:
            raw[anchor] = raw[512]
    return interpolated_anchors(interpolated_anchors(raw))


def rank_polynomial_power(local, packets, degree):
    """Positive matrix-polynomial products, retaining every needed coefficient."""
    size = len(local[0])
    result = [[arb_poly([int(i == j)]) for j in range(size)] for i in range(size)]
    factor = [[arb_poly([ball(local[w][i][j])*math.comb(16, w)
                         for w in range(min(16, degree)+1)])
               for j in range(size)] for i in range(size)]
    def product(left, right):
        return [[sum((left[i][k]*right[k][j] for k in range(size)), arb_poly([])).truncate(degree+1)
                 for j in range(size)] for i in range(size)]
    while packets:
        if packets & 1:
            result = product(result, factor)
        packets >>= 1
        if packets:
            factor = product(factor, factor)
    return tuple(value for row in result for value in row)


def rank_terminal_power(flat, exponent):
    size = math.isqrt(len(flat))
    if size*size != len(flat):
        raise ValueError('square flattened matrix required')
    matrix = arb_mat([[flat[i*size+j] for j in range(size)] for i in range(size)])
    result = arb_mat([[int(i == j) for j in range(size)] for i in range(size)])
    while exponent:
        if exponent & 1:
            result = result*matrix
        exponent >>= 1
        if exponent:
            matrix = matrix*matrix
    return sum((result[0, j] for j in range(size)), arb(0))


def certify_zone(lower, upper, anchor, point, k_log, dimension,
                 cell_bits=48, precision=4096, split=True, rank_history=False):
    g = rate_geometry(k_log, dimension)
    if not 1 <= lower <= upper <= g['rows'] or cell_bits < 1:
        raise ValueError('invalid replay interval')
    y, tau = rational_marks(anchor, point, g['n'])
    previous = ctx.prec
    try:
        ctx.prec = precision
        if rank_history:
            from gf8_rank_history_screen import local_matrix
            local = [local_matrix(w, y, tau, state=8) for w in range(17)]
            coefficients = rank_polynomial_power(local, g['rows']//16, upper)
            terminal = rank_terminal_power
        else:
            local = [packet(w, y, tau, joint=True) for w in range(17)]
            coefficients = polynomial_power(local, g['rows']//16, upper)
            terminal = terminal_power
        output_factor = ball(y)**(-g['threshold'])
        target = arb(2)**(-cell_bits)
        cached = {}
        def bound(lo, hi):
            if hi not in cached:
                cached[hi] = terminal(region(coefficients, hi, g['rows'], tau), 128)
            return (source_maximum(lo, hi, g['rows'], 128-dimension, tau)
                    * output_factor * cached[hi] * (hi-lo+1))
        stack, entries, total = [(lower, upper)], [], arb(0)
        while stack:
            lo, hi = stack.pop()
            endpoint = bound(lo, hi).upper()
            if not endpoint.is_finite() or not endpoint > 0:
                raise ArithmeticError('insufficient Arb precision')
            if endpoint <= target or not split:
                entries.append(dict(lower=lo, upper=hi, anchor=anchor,
                                    certified_integer_bits=integer_bits(endpoint)))
                total += endpoint
            elif lo == hi:
                raise ArithmeticError(f'active={lo}: {integer_bits(endpoint)} bits below cell target{cell_bits}')
            else:
                mid = (lo+hi)//2
                stack.extend(((mid+1, hi), (lo, mid)))
        complete_cover(entries, lower, upper)
        endpoint = total.upper()
        exact = endpoint.fmpq()
        return dict(anchor=anchor, y=str(y), tau=str(tau), precision_bits=precision,
                    rank_history=rank_history,
                    entries=entries, certified_integer_bits=integer_bits(endpoint),
                    failure_upper_dyadic=[str(exact.p), str(exact.q)])
    finally:
        ctx.prec = previous


def certify_point(active, k_log, dimension, precision=4096, point=None, rank_history=False):
    marks = witnesses(active)
    anchor = min(marks, key=lambda a: abs(a-active))
    if point is None:
        point = marks[anchor]
    # Match the supplied tilt at the requested active count.
    anchor = active
    receipt = certify_zone(active, active, anchor, point, k_log, dimension,
                           precision=precision, split=False, rank_history=rank_history)
    return dict(**metadata(k_log, dimension), active_range=[active, active],
                complete=True, **receipt)


def certify_cover(last, k_log, dimension, precision=4096, cell_bits=48, progress=None, rank_history=False):
    g = rate_geometry(k_log, dimension)
    if not 1 <= last <= g['rows']:
        raise ValueError('cover endpoint outside row range')
    anchors = witnesses(last)
    entries, receipts, total = [], [], F(0)
    for lo, hi, anchor in anchor_zones(1, last, anchors):
        receipt = certify_zone(lo, hi, anchor, anchors[anchor], k_log, dimension,
                               precision=precision, cell_bits=cell_bits, rank_history=rank_history)
        receipts.append(receipt)
        entries.extend(receipt['entries'])
        p, q = map(int, receipt['failure_upper_dyadic'])
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
    return dict(**metadata(k_log, dimension), active_range=[1, last], complete=True,
                certified_integer_bits=bits, margin_bits_lower=bits,
                interval_count=len(entries), partition_sha256=digest,
                failure_upper_dyadic=[str(total.numerator), str(total.denominator)],
                receipts=receipts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension', type=int, nargs='+', default=(88, 89, 90))
    parser.add_argument('--k-log', type=int, nargs='+', default=(20,))
    parser.add_argument('--active', type=int, nargs='+', default=(1, 16, 64, 128, 256, 512))
    parser.add_argument('--cover', type=int)
    parser.add_argument('--precision', type=int, default=4096)
    parser.add_argument('--cell-bits', type=int, default=48)
    parser.add_argument('--quiet', action='store_true')
    parser.add_argument('--rank-history', action='store_true')
    args = parser.parse_args()
    for k_log in args.k_log:
        for dimension in args.dimension:
            if args.cover is not None:
                result = certify_cover(args.cover, k_log, dimension, args.precision,
                    args.cell_bits, None if args.quiet else lambda r: print(json.dumps(r), flush=True),
                    args.rank_history)
                del result['receipts']
                print(json.dumps(result), flush=True)
            else:
                for active in args.active:
                    result = certify_point(active, k_log, dimension, args.precision, rank_history=args.rank_history)
                    del result['entries'], result['failure_upper_dyadic']
                    print(json.dumps(result), flush=True)
