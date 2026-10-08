"""One-route design screens and support operators; not full certificates."""

import argparse
from fractions import Fraction as F
import json
import math


def logadd(values):
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def logchoose(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def causal_suffix_probability(n, forbidden_weight, length, dimension):
    """Exact fixed-row event forcing failure of any causal inner.

    A GRS row has a nonzero shortening inside the suffix whenever at
    least length-dimension+1 of its coordinates land there.
    This is one-row probability, not a union over repeated rows.
    """
    minimum = length - dimension + 1
    numerator = sum(math.comb(forbidden_weight, count)
                    * math.comb(n - forbidden_weight, length - count)
                    for count in range(minimum, min(length, forbidden_weight) + 1))
    return F(numerator, math.comb(n, length))


def one_scan_run_upper(n, s, e, g):
    """Exact rational upper count H_g(s,e), using monotone placements.

    No normalization or source multiplicity is included. Closed runs
    use binom(n,r) in place of binom(n-h,r); open runs similarly.
    """
    total = F(0)
    for runs in range(1, s // 2 + 1):
        if s - runs <= e:
            total += (math.comb(s - runs - 1, runs - 1) * math.comb(n, runs)
                      * math.comb(e, s - runs) * g ** runs)
    for runs in range(1, (s + 1) // 2 + 1):
        if s - runs + 1 <= e:
            total += (math.comb(s - runs, runs - 1) * math.comb(n, runs - 1)
                      * math.comb(e, s - runs + 1) * g ** (runs - 1))
    return total


def log_one_scan_run_upper(n, s, e, g):
    values = []
    for runs in range(1, s // 2 + 1):
        values.append(logchoose(s - runs - 1, runs - 1) + logchoose(n, runs)
                      + logchoose(e, s - runs) + runs * math.log(g))
    for runs in range(1, (s + 1) // 2 + 1):
        values.append(logchoose(s - runs, runs - 1) + logchoose(n, runs - 1)
                      + logchoose(e, s - runs + 1) + (runs - 1) * math.log(g))
    return logadd(values)


def one_scan_q1_screen(k=2 ** 16, length=128, dimension=64, alphabet=255,
                       distance=0.1):
    assert k % dimension == 0
    rows = k // dimension
    n, e = rows * length, math.ceil(distance * rows * length) - 1
    best = None
    for index in range(801):
        tau = math.exp(index * math.log(4096) / 800)
        g = 1 / alphabet + (1 - 1 / alphabet) / tau
        terms = []
        for weight in range(length - dimension + 1, length + 1):
            terms.append(logchoose(length, weight)
                         + (weight - (length - dimension) - 1) * math.log(tau)
                         + log_one_scan_run_upper(n, weight, e, g)
                         - logchoose(n, weight))
        bits = -(math.log(rows) + logadd(terms)) / math.log(2)
        if best is None or bits > best['margin_bits']:
            best = dict(k=k, outer=[length, dimension], alphabet=alphabet,
                        tau=tau, margin_bits=bits, rigorous=False,
                        scope='one-active-row rank bound only')
    return best


def scan_step(previous, input_nonzero):
    """Admissible (new state, number of returns) of one scalar scan."""
    if previous == 0:
        return [(input_nonzero, 0)]
    if input_nonzero == 0:
        return [(1, 0)]
    return [(0, 1), (1, 0)]


def cascade_operator(x, z, g):
    """Two same-order weighted scans, without an intermediate permutation.

    States are 00,01,10,11 for the first/second scan's previous outputs.
    x marks original inputs, z final outputs, g each return rank charge.
    The entries are support-path weights, not probabilities.
    """
    matrix = [[0 for _ in range(4)] for _ in range(4)]
    for old in range(4):
        first, second = old >> 1, old & 1
        for source in (0, 1):
            for next_first, first_returns in scan_step(first, source):
                for next_second, second_returns in scan_step(second, next_first):
                    matrix[old][2 * next_first + next_second] += (
                        x ** source * z ** next_second * g ** (first_returns + second_returns))
    return matrix


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen', action='store_true')
    parser.add_argument('--causal-suffix', action='store_true')
    args = parser.parse_args()
    if args.screen:
        for length in (16, 32, 64, 128):
            print(json.dumps(one_scan_q1_screen(length=length, dimension=length // 2)), flush=True)
    if args.causal_suffix:
        k, n, e = 2 ** 16, 2 ** 17, 13107
        for length in (16, 32, 64, 128):
            dimension = length // 2
            probability = causal_suffix_probability(n, e, length, dimension)
            # The second quantity is an expected count of bad rows,
            # not an asserted probability when row events overlap.
            expected = (k // dimension) * probability
            print(json.dumps(dict(outer=[length, dimension],
                fixed_row_probability=float(probability),
                expected_bad_rows=float(expected),
                expected_bad_rows_log2=math.log2(expected))), flush=True)
