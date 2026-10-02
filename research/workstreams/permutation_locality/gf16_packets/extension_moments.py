"""Exact low-degree orthogonal-array bounds for shared BCH union supports.

A binary code extended to GF(16) identifies its words with ordered binary
quadruples. Symbol weight is their union support. A dual distance exceeding
2r makes every degree-2r weight moment equal to the ambient-binomial moment.
The squared reproducing kernel below is nonnegative at every weight; its
minimum on the requested support interval is checked over the integers.
This is an optional counting experiment, not a new distance certificate.
"""
import argparse
import json
from fractions import Fraction as Q
from math import comb, log2
from pathlib import Path


def krawtchouk(n, q, x, degree):
    if (any(type(v) is not int for v in (n, q, x, degree))
            or n < 1 or q < 2 or not 0 <= x <= n or not 0 <= degree <= n):
        raise ValueError('valid integer q-ary polynomial parameters required')
    row = [1]
    if degree:
        row.append((q-1)*n-q*x)
    for j in range(1, degree):
        numerator = ((q-1)*n-q*x-(q-2)*j)*row[j]-(q-1)*(n-j+1)*row[j-1]
        value, remainder = divmod(numerator, j+1)
        assert remainder == 0
        row.append(value)
    return row


def support_cap(n, q, size, dual_distance, minimum, cutoff, degree=None):
    """Return a rationally verified integer upper bound, assuming the premises.

    No saved count is a premise. The caller must establish size, minimum
    distance, and dual distance independently. No assertion of tightness.
    """
    if (any(type(v) is not int for v in (size, dual_distance, minimum, cutoff))
            or size < 1 or not 1 <= minimum <= cutoff <= n
            or not 1 <= dual_distance <= n+1):
        raise ValueError('valid code size and distance parameters required')
    degree = (dual_distance-1)//2 if degree is None else degree
    if type(degree) is not int or not 0 <= degree <= n or 2*degree >= dual_distance:
        raise ValueError('twice the polynomial degree must be below dual distance')
    anchor = krawtchouk(n, q, cutoff, degree)
    norms = [comb(n, j)*(q-1)**j for j in range(degree+1)]
    coefficients = [Q(v, norm) for v, norm in zip(anchor, norms)]
    moment = sum((Q(v*v, norm) for v, norm in zip(anchor, norms)), Q(0))
    def polynomial(x):
        return sum((c*v for c, v in zip(coefficients, krawtchouk(n, q, x, degree))), Q(0))
    lower = min(polynomial(x)**2 for x in range(minimum, cutoff+1))
    if not lower:
        return size-1
    budget = size*moment-polynomial(0)**2
    if budget < 0:
        raise ArithmeticError('code premises contradict the exact moment identity')
    bound = budget/lower
    return min(size-1, bound.numerator//bound.denominator)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supports', type=int, nargs='+', default=[96, 104, 108, 112, 116, 120, 128])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    # shared_support loads the existing dependency paths and authenticates
    # both the production binary generator and the BCH spectral premises.
    from shared_support import shared_counts
    from dual_shortening import verify_bch_premise
    cdf, _ = shared_counts(True, True)
    verify_bch_premise()
    rows = []
    for u in args.supports:
        candidate = support_cap(256, 16, 1 << 512, 30, 38, u)
        rows.append(dict(support=u, existing_cap=str(cdf[u]), moment_cap=str(candidate),
                         improvement_bits=max(0., log2(cdf[u])-log2(candidate))))
        print('GF16 MOMENT support', u, 'existing log2', log2(cdf[u]),
              'moment log2', log2(candidate), 'gain', rows[-1]['improvement_bits'], flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(dict(schema='shared-gf16-extension-moments-1',
            degree=14, dual_distance=30, rows=rows), indent=2)+'\n')


if __name__ == '__main__':
    main()
