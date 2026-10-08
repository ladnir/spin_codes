"""Exact dense-input class certificate for one regional routing stage.

This is not a whole-code certificate. It covers inverse-inner inputs above
the displayed density, for any invertible inner independent of the route
and the external nonzero coordinate diagonals. Floating optimization only
chooses rational tilts; every accepted inequality is checked with integers.
No native encoding or timing is performed.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
import json
import math

from scipy.optimize import minimize_scalar
from scipy.special import logsumexp


PRESETS = ((32, 24, 60), (33, 24, 56), (32, 23, 54),
           (32, 22, 48), (34, 24, 50), (64, 48, 60), (36, 24, 50))
GRID = 200


def interval_base(length, dimension, lower, upper, tilt, q=256):
    """Unnormalized row bracket, rounded upward across one density interval."""
    if not 0 <= lower <= upper <= 1 or tilt < 1:
        raise ValueError('require an ordered density interval and tilt >= 1')
    return ((1-lower)**length
            + Fraction(1, q**(length-dimension))
            * (1+(Fraction(q, q-1)*tilt-1)*upper)**length)


def choose_tilt(length, dimension, index):
    """Choose a rational witness; approximation errors cannot certify a bound."""
    lower, upper = index/GRID, (index+1)/GRID
    def objective(log_tilt):
        tilt = math.exp(log_tilt)
        return (float(logsumexp([
            length*math.log1p(-lower),
            -(length-dimension)*math.log(256)
            + length*math.log1p((256/255*tilt-1)*upper)]))
            - length*lower*log_tilt)
    result = minimize_scalar(objective, bounds=(0, 6), method='bounded')
    return max(Fraction(1), Fraction(math.exp(result.x)).limit_denominator(16))


def certify_density_interval(length, dimension, index, tilt):
    """Prove membership exponent >= 1.28 on [index/200,(index+1)/200]."""
    lower, upper = Fraction(index, GRID), Fraction(index+1, GRID)
    base = interval_base(length, dimension, lower, upper, tilt)
    # (tilt^(-lower*b) * base)^(1/b) <= 2^(-32/25).
    # Raising to 200 removes all fractional powers.
    return base**GRID * 2**(256*length) <= tilt**(index*length)


def certify_dense_class(length, dimension, first_interval, period=64):
    tilts = []
    for index in range(first_interval, GRID):
        tilt = choose_tilt(length, dimension, index)
        if not certify_density_interval(length, dimension, index, tilt):
            raise ArithmeticError(f'uncertified density interval {index}/{GRID}')
        tilts.append((index, tilt.numerator, tilt.denominator))
    # At delta <= .1, the q-ary ball has size at most
    # [(10/9) * 2295^(1/10)]^N <= 2^(1.27*N).
    assert Fraction(10, 9)**100 * 2295**10 < 2**127
    sizes = []
    for k_log in (16, 18, 20):
        multiple = dimension*period
        k = ((2**k_log+multiple-1)//multiple)*multiple
        rows, n = k//dimension, k*length//dimension
        # (rows+1)^b * 2^(-N/100) <= 2^(-margin).
        margin = n//100-length*rows.bit_length()
        assert margin >= 40
        sizes.append(dict(k_log=k_log, k=k, rows=rows, n=n,
                          dense_class_margin_bits_lower=margin))
    return dict(outer=[length, dimension], q=256, period=period,
                inverse_input_density_min=str(Fraction(first_interval, GRID)),
                membership_exponent='1.28', ball_exponent='1.27',
                rigorous=True, scope='dense-inverse-input-class-only',
                sizes=sizes, rational_tilts=tilts)


def self_test():
    """Exhaust small exact profile probabilities independently of the tilt proof."""
    # The GF(4) [3,2] parity code x+y+z=0 has these prescribed-support counts.
    q, length, dimension, rows = 4, 3, 2, 3
    support_counts = (1, 0, 3, 6)
    subsets = {weight: tuple(combinations(range(rows), weight))
               for weight in range(rows+1)}
    checked = 0
    for profile in product(range(rows+1), repeat=length):
        total, choices = Fraction(0), 0
        for columns in product(*(subsets[weight] for weight in profile)):
            probability = Fraction(1)
            for row in range(rows):
                weight = sum(row in column for column in columns)
                probability *= Fraction(support_counts[weight], (q-1)**weight)
            total += probability
            choices += 1
        actual = total/choices
        densities = tuple(Fraction(weight, rows) for weight in profile)
        weight = sum(profile)
        alpha = Fraction(weight, length*rows)
        for tilt in (Fraction(1, 2), Fraction(1), Fraction(2), Fraction(5)):
            exact_profile_bracket = (math.prod(1-a for a in densities)
                + Fraction(1, q**(length-dimension))
                * math.prod(1-a+Fraction(q, q-1)*tilt*a for a in densities))
            symmetric_bracket = ((1-alpha)**length
                + Fraction(1, q**(length-dimension))
                * (1-alpha+Fraction(q, q-1)*tilt*alpha)**length)
            profile_bound = ((rows+1)**length * tilt**(-weight)
                             * exact_profile_bracket**rows)
            symmetric_bound = ((rows+1)**length * tilt**(-weight)
                               * symmetric_bracket**rows)
            assert actual <= profile_bound <= symmetric_bound
            checked += 1
    assert Fraction(10, 9)**100 * 2295**10 < 2**127
    return checked


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--show-tilts', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(dict(exact_profile_checks=self_test())), flush=True)
    for length, dimension, start in PRESETS:
        result = certify_dense_class(length, dimension, start)
        if not args.show_tilts:
            result['certified_density_intervals'] = len(result.pop('rational_tilts'))
        print(json.dumps(result), flush=True)
