"""Exact obstruction for a uniform interleaver followed by a causal inner.

For L disjoint MDS[b,k] outer blocks, N=L*b. If at least b-k+1
coordinates of one block fall in the last e positions, an outer word and
therefore its causal-inner output are supported in that suffix.

The block-event indicators have nonpositive pairwise covariance. Their
second moment gives P[bad code] >= L*p / (1+(L-1)*p), where p is the
exact hypergeometric block-event probability. This is a LOWER bound on
failure, not a distance certificate. Decimal/log displays are diagnostic;
comparisons with a proposed failure target use exact integers.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
import json
import math


@lru_cache(maxsize=128)
def tail_probability(length, suffix, block, minimum_support):
    """P[at least minimum_support of block labels lie in a fixed suffix]."""
    if not (0 <= suffix <= length and 1 <= block <= length
            and 1 <= minimum_support <= block):
        raise ValueError('invalid suffix or block geometry')
    lo = max(minimum_support, block-(length-suffix))
    hi = min(block, suffix)
    if lo > hi:
        return Fraction(0)
    term = math.comb(suffix, lo)*math.comb(length-suffix, block-lo)
    numerator = term
    for h in range(lo, hi):
        product = term*(suffix-h)*(block-h)
        denominator = (h+1)*(length-suffix-block+h+1)
        term, remainder = divmod(product, denominator)
        if remainder:
            raise ArithmeticError('hypergeometric integer recurrence failed')
        numerator += term
    return Fraction(numerator, math.comb(length, block))


def obstruction(rows, block, dimension, suffix):
    if not (rows >= 1 and 1 <= dimension <= block):
        raise ValueError('invalid MDS outer geometry')
    length = rows*block
    p = tail_probability(length, suffix, block, block-dimension+1)
    failure = rows*p/(1+(rows-1)*p)
    return p, failure


def margin_ceiling(failure):
    if failure == 0:
        return math.inf
    return math.log2(failure.denominator)-math.log2(failure.numerator)


def binomial_tail(block, minimum_support, probability):
    return sum((Fraction(math.comb(block, h))*probability**h
                *(1-probability)**(block-h)
                for h in range(minimum_support, block+1)), Fraction(0))


def success_bound_exponent(rows, block, dimension, suffix):
    """Exact integer a with P[no suffix witness] <= min(1,(N+1)*2^-a).

    Under independent Bernoulli e/N coordinate sampling, block events are
    independent. Conditioning on exactly e sampled coordinates costs at
    most N+1 because e is a binomial mode. Use 1-p<=exp(-p) and ln(2)<7/10;
    no floating logarithm is needed to obtain the displayed dyadic bound.
    """
    length = rows*block
    if not (rows >= 1 and 1 <= dimension <= block and 0 <= suffix <= length):
        raise ValueError('invalid suffix or MDS outer geometry')
    p = binomial_tail(block, block-dimension+1, Fraction(suffix, length))
    exponent = rows*p*Fraction(10, 7)
    return exponent.numerator//exponent.denominator


def screen(k, block, dimension, delta, bits=40):
    rows = (k+dimension-1)//dimension
    length = rows*block
    # A target d_min >= ceil(delta*N) forbids all weights <= e below it.
    e = (delta.numerator*length+delta.denominator-1)//delta.denominator-1
    p, failure = obstruction(rows, block, dimension, e)
    exponent = success_bound_exponent(rows, block, dimension, e)
    return {'requested_k': k, 'k': rows*dimension, 'n': length,
            'block': block, 'dimension': dimension, 'rows': rows,
            'forbidden_suffix': e, 'block_event_probability': float(p),
            'expected_bad_blocks': float(rows*p),
            'failure_lower_bound': float(failure),
            'maximum_possible_margin_bits': margin_ceiling(failure),
            'no_suffix_witness_probability_upper': f'min(1, {length+1} * 2^-{exponent})',
            'rules_out_target': failure > Fraction(1, 1 << bits)}


def first_length_not_excluded(k, rate, delta, bits=40, maximum=1024):
    """Necessary-only length screen: passing says nothing about sufficiency."""
    for block in range(rate.denominator, maximum+1, rate.denominator):
        row = screen(k, block, block*rate.numerator//rate.denominator, delta, bits)
        if not row['rules_out_target']:
            return row
    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log-k', type=int, nargs='+', default=[16, 18, 20])
    parser.add_argument('--delta', type=Fraction, default=Fraction(3, 10))
    parser.add_argument('--rate', type=Fraction, default=Fraction(1, 2))
    parser.add_argument('--blocks', type=int, nargs='+', default=[16, 64, 128, 256])
    parser.add_argument('--bits', type=int, default=40)
    parser.add_argument('--find-minimum', action='store_true')
    args = parser.parse_args()
    if not (0 < args.delta < 1 and 0 < args.rate < 1 and args.bits >= 0):
        parser.error('require 0<delta,rate<1 and bits>=0')
    for exponent in args.log_k:
        if args.find_minimum:
            print(json.dumps(first_length_not_excluded(
                1 << exponent, args.rate, args.delta, args.bits)))
        else:
            for block in args.blocks:
                if block % args.rate.denominator:
                    parser.error('each block length must support the stated rate exactly')
                print(json.dumps(screen(1 << exponent, block,
                                        block*args.rate.numerator//args.rate.denominator,
                                        args.delta, args.bits)))
