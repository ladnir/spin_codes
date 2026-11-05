"""Bounded replay of the retained Toeplitz-window certificate at higher distance.

No construction or probability bound is changed. A FLOAT cover proposes boxes;
each directed stage regenerates that cover and verifies every box independently.
Run one stage at a time, with an external wall-clock/memory limit as well.
"""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import time

import stratified_batch_certificate as retained
from stratified_batch_certificate import cover, BatchUpward


def check_partition(boxes, n, m, deadline):
    """Check exact adjacency on every c slice, independently of cover's traversal."""
    cmax = (n + m) // (m + 1)
    canonical = []
    for vl, vh, cl, ch, marker, _ in boxes:
        tau = Fraction(*marker)
        if not 1 <= vl <= vh <= n or not 1 <= cl <= ch <= cmax or not 1 < tau <= 2**128:
            raise ArithmeticError('Invalid certificate box')
        canonical.append((vl, vh, cl, ch, tau.numerator, tau.denominator))
    total = 0
    for c in range(1, cmax + 1):
        if time.monotonic() > deadline:
            raise RuntimeError('Partition check exceeded time limit')
        lower = (m + 1) * c - m
        spans = sorted((max(vl, lower), vh) for vl, vh, cl, ch, _, _ in canonical
                       if cl <= c <= ch and vh >= lower)
        cursor = lower
        for lo, hi in spans:
            if lo != cursor:
                raise ArithmeticError('Certificate partition has a hole or overlap')
            total += hi - lo + 1
            cursor = hi + 1
        if cursor != n + 1:
            raise ArithmeticError('Certificate partition omits a suffix')
    expected = cmax * (n + m + 1) - (m + 1) * cmax * (cmax + 1) // 2
    if total != expected:
        raise ArithmeticError('Certificate pair count disagrees with exact formula')
    digest = hashlib.sha256(json.dumps(sorted(canonical), separators=(',', ':')).encode()).hexdigest()
    return total, digest


def check_negative_controls(boxes, n, m, deadline):
    changed = list(boxes)
    vl, vh, cl, ch, marker, value = changed[0]
    changed[0] = (vl + 1, vh, cl, ch, marker, value)
    for bad in (boxes[1:], boxes + boxes[:1], changed):
        try:
            check_partition(bad, n, m, deadline)
        except ArithmeticError:
            continue
        raise ArithmeticError('Partition checker accepted a corrupted cover')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--L', type=int, default=256)
    parser.add_argument('--batch', type=int, default=64)
    parser.add_argument('--distance', type=Fraction, required=True)
    parser.add_argument('--term-bits', type=int)
    parser.add_argument('--precision', type=int, choices=(80, 120))
    parser.add_argument('--fine-markers', action='store_true')
    parser.add_argument('--check-controls', action='store_true')
    parser.add_argument('--point', type=int, nargs=2, metavar=('V', 'C'))
    parser.add_argument('--seconds', type=float, default=45)
    args = parser.parse_args()
    N, B = 2**21, 256
    if (args.L < 1 or args.batch < 1 or args.L % args.batch or N % args.batch
            or not 0 < args.distance < Fraction(1, 2)
            or (args.term_bits is not None and args.term_bits < 1)
            or not 0 < args.seconds <= 50):
        parser.error('Invalid dimensions, distance, precision, or time limit')
    n, m = N // args.batch, args.L // args.batch
    states = n * ((n + m) // (m + 1))
    if args.term_bits is None:
        args.term_bits = 40 + (states - 1).bit_length()
    margin = args.term_bits - (states - 1).bit_length()
    if margin < 40:
        parser.error('Per-term threshold does not provide a 40-bit union bound')
    if args.fine_markers:
        retained.MARKERS = tuple((32 + i, 32) for i in range(1, 33)) + ((4, 1), (8, 1))
    e = args.distance.numerator * N // args.distance.denominator
    if args.point is not None:
        v, c = args.point
        if not 1 <= c <= (n + m) // (m + 1) or not m * (c - 1) + c <= v <= n:
            parser.error('Infeasible point')
        value, marker = retained.float_box(N, B, args.L, args.batch, e,
                                          v, v, c, c, entropy_volume=True)
        print(dict(kind='floating_point_not_certificate', L=args.L, batch=args.batch,
                   distance=str(args.distance), v=v, c=c, log2=value, marker=marker), flush=True)
        return
    started = time.monotonic()
    boxes = cover(B, args.L, args.batch, args.distance, args.term_bits,
                  entropy_volume=True, seconds=args.seconds)
    deadline = started + args.seconds
    pairs, digest = check_partition(boxes, n, m, deadline)
    if args.check_controls:
        check_negative_controls(boxes, n, m, deadline)
    print(dict(kind='exact_integer_partition', pairs=pairs, sha256=digest,
               corrupted_cover_controls=args.check_controls), flush=True)
    if args.precision is None:
        return
    bound = BatchUpward(args.precision)
    worst = Decimal('-Infinity')
    with localcontext(bound.up):
        for item in boxes:
            if time.monotonic() - started > args.seconds:
                raise RuntimeError('Combined cover/replay time limit reached; no certificate')
            value = bound.rectangle(N, B, args.L, args.batch, e, item,
                                    entropy_volume=True)
            if not value <= -args.term_bits:
                raise ArithmeticError(f'Directed rectangle did not pass: {item}, {value}')
            worst = max(worst, value)
    if time.monotonic() > deadline:
        raise RuntimeError('Final directed rectangle exceeded time limit; no certificate')
    print(dict(kind='directed_complete_high_distance_cover', B=B, L=args.L,
               batch=args.batch, N=N, threshold=e, distance=str(args.distance),
               precision=args.precision, boxes=len(boxes), covered_pairs=pairs,
               partition_sha256=digest,
               worst_log2_upper=str(worst), margin_at_least_bits=margin,
               seconds=round(time.monotonic() - started, 3)), flush=True)


if __name__ == '__main__':
    main()
