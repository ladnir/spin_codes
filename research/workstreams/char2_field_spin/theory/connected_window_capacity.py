"""Bound one-window-per-region dimension obstructions, including joins.

This is NOT a minimum-distance certificate. The event is D>s*(c-tail),
where D is the shortened outer dimension, c counts connected intervals,
and tail is one when the support reaches the global stream end.

Floating optimization selects rational Chernoff marks only. --replay uses
outward decimal intervals for every probability and the complete union.
"""

import argparse
from fractions import Fraction
import json
import math

from scipy.optimize import minimize_scalar

from subfield_certificate import Interval, integer_margin
from window_dimension_mgf import log_moment


def choose(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def shape_count(used, components, singles, tail_kind=0):
    """Ordered positive component lengths; 0=no tail, 1=single, 2=multi."""
    multis = components-singles
    if not (0 <= singles <= components <= used and components > 0):
        return 0
    if multis == 0:
        lengths = int(used == components)
    else:
        lengths = choose(used-components-1, multis-1)
    if tail_kind == 0:
        positions = choose(components, singles)
    elif tail_kind == 1:
        positions = choose(components-1, singles-1)
    elif tail_kind == 2:
        positions = choose(components-1, singles)
    else:
        raise ValueError('unknown tail kind')
    return positions*lengths


def families(regions, parity, units, budget):
    """Finite covering counts; length caps and non-touching exclusions dropped.

    Each multi-region component has two variable endpoints and full
    interior regions. At the global tail its final endpoint is full too.
    Single-region components have at most units choices of start, except
    a final suffix whose start is fixed by its length.
    """
    for used in range(parity+1, regions+1):
        for components in range(1, used+1):
            for singles in range(max(0, 2*components-used), components+1):
                for tail_kind in range(3):
                    shapes = shape_count(used, components, singles, tail_kind)
                    if not shapes:
                        continue
                    tail = int(tail_kind != 0)
                    variables = 2*components-singles-int(tail_kind == 2)
                    full = used-variables
                    available = budget-full*units
                    if available < variables or variables <= 0:
                        continue
                    placements = choose(regions-used+components-tail, components-tail)
                    count = (shapes*placements*choose(available, variables)
                             *units**(singles-int(tail_kind == 1)))
                    if count:
                        yield dict(used=used, components=components, singles=singles,
                                   tail_kind=tail_kind, full=full, variables=variables,
                                   available=available, count=count)


def rational_moment(trials, parity, probability, mark):
    """Positive binomial sum in exact Fractions or outward intervals."""
    return sum(choose(trials, hits)*probability**hits
               *(1-probability)**(trials-hits)*mark**max(hits-parity, 0)
               for hits in range(trials+1))


def evaluate_family(family, rows, parity, state, units, replay=False):
    trials = family['variables']
    residual_parity = parity-family['full']
    probability = min(Fraction(1), Fraction(family['available'], units*trials))
    crossing = state*(family['components']-int(family['tail_kind'] != 0))+1
    if residual_parity < 0:
        # Not needed for the present 10% / rate<=3/4 candidates. Retain a
        # safe unoptimized bound rather than silently drop these families.
        mark = Fraction(1)
        log_fixed = 0.
    else:
        objective = lambda lam: (rows*log_moment(trials, residual_parity,
                                                float(probability), lam)-crossing*lam)
        result = minimize_scalar(objective, bounds=(0., min(20., 650./trials)),
                                 method='bounded', options={'xatol': 1e-10})
        mark = Fraction(math.exp(float(result.x))).limit_denominator(10**6)
        log_fixed = min(0., objective(math.log(float(mark))))
    entry = dict(family, probability=[probability.numerator, probability.denominator],
                 mark=[mark.numerator, mark.denominator], crossing=crossing,
                 margin_bits=-(math.log(family['count'])+log_fixed)/math.log(2))
    if not replay:
        return entry, None
    pi = Interval.rational(probability.numerator, probability.denominator)
    zi = Interval.rational(mark.numerator, mark.denominator)
    fixed = rational_moment(trials, residual_parity, pi, zi)**rows/zi**crossing
    if fixed.hi > 1:
        fixed = Interval.integer(1)
    upper = fixed*family['count']
    entry['failure_upper'] = str(upper.hi)
    return entry, upper


def screen(rows, regions, dimension, state, alignment=16, replay=False, detailed=False):
    if not (rows > 0 and rows % alignment == 0 and 0 < dimension < regions and state > 0):
        raise ValueError('invalid geometry')
    parity, units = regions-dimension, rows//alignment
    threshold = (rows*regions+9)//10-1
    budget = threshold//alignment
    entries = []
    total = Interval.integer(0)
    for family in families(regions, parity, units, budget):
        entry, upper = evaluate_family(family, rows, parity, state, units, replay)
        entries.append(entry)
        if upper is not None:
            total += upper
    if entries:
        peak = max(-entry['margin_bits'] for entry in entries)
        margin = -peak-math.log2(sum(2**(-entry['margin_bits']-peak) for entry in entries))
    else:
        margin = math.inf
    output = dict(outer=[regions, dimension], rows=rows, state=state,
                  bad_output_limit=threshold, alignment=alignment, margin_bits=margin,
                  family_count=len(entries), rigorous=replay, whole_code=False,
                  scope='shortened dimension exceeds connected-window exit capacity',
                  worst_entries=sorted(entries, key=lambda entry: entry['margin_bits'])[:5])
    if replay:
        output['failure_upper'] = str(total.hi)
        output['certified_integer_bits'] = integer_margin(total.hi) if total.hi else None
        output['zero_event'] = not entries
    if detailed:
        output['entries'] = entries
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', type=int, required=True)
    parser.add_argument('--outer', type=int, nargs=2, required=True)
    parser.add_argument('--state', type=int, default=8)
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--detailed', action='store_true')
    args = parser.parse_args()
    print(json.dumps(screen(args.rows, *args.outer, args.state,
                            replay=args.replay, detailed=args.detailed)), flush=True)
