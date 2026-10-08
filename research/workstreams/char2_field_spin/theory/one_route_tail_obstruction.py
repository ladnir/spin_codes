"""Exact late-activation obstruction for one uniform global permutation.

This bounds actual code-failure probability from below, not the expected
number of bad words. Assumptions: a repeated MDS[b,k] outer, a uniform
permutation of all N coordinates, and a zero-initialized causal linear inner.
The final window must begin on a whole inner-step boundary. Nonzero input
diagonals and arbitrary causal state updates do not remove this obstruction.
"""

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import comb, log2


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def row_tail(total, window, block, minimum):
    """Probability that a fixed block places at least minimum slots in window."""
    if not 0 <= window <= total or not 0 < block <= total:
        raise ValueError('invalid population')
    return sum((Fraction(choose(window, j)*choose(total-window, block-j),
                         comb(total, block))
                for j in range(minimum, block+1)), Fraction(0))


def pair_tail(total, window, block, minimum):
    """Joint event for two disjoint labeled outer blocks."""
    if total < 2*block:
        raise ValueError('two disjoint blocks required')
    result = Fraction(0)
    for j in range(minimum, block+1):
        first = Fraction(choose(window, j)*choose(total-window, block-j),
                         comb(total, block))
        if first:
            result += first*row_tail(total-block, window-j, block, minimum)
    return result


def obstruction(k=2**16, block=16, dimension=8, step=16,
                distance=Fraction(1, 10)):
    if not 0 < dimension < block or k % dimension or step < 1:
        raise ValueError('outer dimensions or step width')
    rows = k//dimension
    n = rows*block
    if n % step or not 0 < distance < 1:
        raise ValueError('aligned output length and relative distance required')
    target = -(-n*distance.numerator//distance.denominator)
    window = ((target-1)//step)*step
    minimum = block-dimension+1
    p = row_tail(n, window, block, minimum)
    joint = pair_tail(n, window, block, minimum) if rows >= 2 else Fraction(0)
    # A subset of rows is enough for a lower bound. Optimize the quadratic
    # Bonferroni expression exactly when using all rows would become vacuous.
    center = min(rows, max(1, int(p/joint))) if joint else rows
    candidates = {1, rows, center, min(rows, center+1)}
    used, lower = max(((a, a*p-comb(a, 2)*joint) for a in candidates),
                      key=lambda item: item[1])
    return dict(k=k, n=n, block=block, dimension=dimension, step=step,
                distance=distance, target=target, window=window,
                minimum_tail_coordinates=minimum, rows=rows, used_rows=used,
                one_row_probability=p, two_row_probability=joint,
                code_failure_lower=lower)


def report(result):
    out = {key: str(value) if isinstance(value, Fraction) else value
           for key, value in result.items()}
    lower = result['code_failure_lower']
    out['code_failure_lower_decimal_diagnostic'] = float(lower)
    out['maximum_possible_margin_bits_diagnostic'] = (
        log2(lower.denominator)-log2(lower.numerator) if lower else None)
    out['scope'] = 'actual code-failure lower bound; uniform global permutation only'
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, default=16)
    parser.add_argument('--blocks', type=int, nargs='+', default=[8, 16, 32, 64])
    parser.add_argument('--step', type=int, default=16)
    args = parser.parse_args()
    for block in args.blocks:
        print(json.dumps(report(obstruction(2**args.k_log, block, block//2,
                                            args.step)), indent=2))
