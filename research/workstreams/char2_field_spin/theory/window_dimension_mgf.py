"""Window-family upper bounds from negative association (not distance).

All rows may have different outer-coordinate supports. For each fixed set
of windows, D=sum_i (hits_i-parity)_+ has an exponential moment bounded by
the product of its one-row moments. Smoothing the window fractions gives
a binomial moment. Unioning all aligned starts and length compositions
covers the stated one-window-per-used-region capacity event only.

The optimizer is FLOATING and cannot certify a margin. Rational replay of
fixed marks is available in moment() and exact_fixed_profile_moment().
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
import json
import math

from scipy.optimize import minimize_scalar


def moment(regions, parity, probability, mark):
    """One-row binomial exponential moment, exact for Fraction arguments."""
    if not (0 <= parity < regions and 0 <= probability <= 1 and mark >= 1):
        raise ValueError('invalid binomial moment')
    return 1+sum(math.comb(regions, count)*probability**count
                 *(1-probability)**(regions-count)
                 *(mark**(count-parity)-1) for count in range(parity+1, regions+1))


def profile_moment(parity, probabilities, mark):
    distribution = [Fraction(1)]
    for probability in probabilities:
        following = [Fraction(0)]*(len(distribution)+1)
        for count, mass in enumerate(distribution):
            following[count] += mass*(1-probability)
            following[count+1] += mass*probability
        distribution = following
    return sum(mass*mark**max(count-parity, 0) for count, mass in enumerate(distribution))


def exact_fixed_profile_moment(rows, lengths, parity, mark):
    """Tiny exact enumeration of independent uniform subsets per region.

    This is deliberately exponential and used only for independent tests.
    Window positions do not affect the distribution for a fixed profile.
    """
    total, cases = Fraction(0), 0
    choices = [tuple(combinations(range(rows), length)) for length in lengths]
    for selected in product(*choices):
        counts = [sum(row in subset for subset in selected) for row in range(rows)]
        total += mark**sum(max(count-parity, 0) for count in counts)
        cases += 1
    return total/cases


def log_moment(regions, parity, probability, log_mark):
    # Do not subtract probabilities whose sum is numerically almost one.
    extra = sum(math.comb(regions, count)*probability**count
                *(1-probability)**(regions-count)
                *math.expm1(log_mark*(count-parity))
                for count in range(parity+1, regions+1))
    return math.log1p(extra)


def screen(rows, regions, dimension, state=4, alignment=16):
    if not (rows % alignment == 0 and 0 < dimension < regions and state > 0):
        raise ValueError('invalid geometry')
    parity = regions-dimension
    threshold = (rows*regions+9)//10-1
    budget = threshold//alignment
    region_units = rows//alignment
    entries = []
    for used in range(parity+1, regions+1):
        # All used lengths are positive integers. Number of ordered positive
        # compositions with total at most budget is C(budget,used).
        log_choices = (math.log(math.comb(regions, used))
                       +math.log(math.comb(budget, used))
                       +used*math.log(region_units))
        probability = min(1., alignment*budget/(rows*used))
        crossing = state*used+1
        objective = lambda lam: (rows*log_moment(used, parity, probability, lam)
                                 -crossing*lam)
        result = minimize_scalar(objective, bounds=(0., min(20., 650./dimension)),
                                 method='bounded', options={'xatol': 1e-11})
        log_upper = min(0., float(result.fun))+log_choices
        entries.append(dict(used_regions=used, margin_bits=-log_upper/math.log(2),
                            fixed_profile_margin_bits=-float(result.fun)/math.log(2),
                            log2_window_choices=log_choices/math.log(2),
                            log_mark=float(result.x)))
    peak = max(-entry['margin_bits'] for entry in entries)
    log2_upper = peak+math.log2(sum(2**(-entry['margin_bits']-peak) for entry in entries))
    return dict(outer=[regions,dimension], rows=rows, state=state,
                bad_output_limit=threshold, alignment=alignment,
                margin_bits=-log2_upper, entries=entries, rigorous=False,
                scope='one aligned window per used region; excess dimension only',
                whole_code=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', type=int, required=True)
    parser.add_argument('--outer', type=int, nargs=2, required=True)
    parser.add_argument('--state', type=int, default=4)
    args = parser.parse_args()
    print(json.dumps(screen(args.rows,*args.outer,args.state)), flush=True)
