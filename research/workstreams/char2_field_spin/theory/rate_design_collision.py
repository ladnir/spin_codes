"""Exact finite two-row collision lower bound for a periodic no-refresh inner.

The outer is repeated systematic GRS[b,k] over GF(q), with one independent
nonzero diagonal per outer symbol and independent symbol permutations in
the b regions. A regional injection label repeats with period B, where
B divides regional length L and is a multiple of inner width. Matching labels cause equal
injections. A pair with equal diagonal ratios on a minimum support returns
to zero within each active region; its output is confined between inputs.

The lower bound uses all row pairs and all minimum supports, then applies
the first two Bonferroni terms. This is an exact rational lower bound on
setup failure for the stated no-refresh periodic construction. It is not
an upper distance certificate and does not apply to independent GL refresh.
"""

import argparse
from fractions import Fraction
import itertools
import json
import math


def choose(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def distance_sum_count(copies, weight, maximum):
    """Ordered-pair position count with sum of positive gaps <= maximum.

    One gap d has 2*(copies-d) ordered choices. The positive polynomial is
    W(x)=sum_(d=1)^(copies-1) (copies-d)*x^d, and
    (1-x)^2 W(x)=(copies-1)x-copies*x^2+x^(copies+1).
    Expanding this short numerator gives an exact cumulative coefficient.
    """
    if copies < 2 or weight < 1:
        raise ValueError('require at least two copies and positive weight')
    if maximum < weight:
        return 0
    if maximum >= weight*(copies-1):
        return (copies*(copies-1))**weight
    total = 0
    for high in range(weight+1):
        for second in range(weight-high+1):
            first = weight-high-second
            degree = first+2*second+(copies+1)*high
            coefficient = choose(weight, high)*choose(weight-high, second)
            coefficient *= (copies-1)**first*(-copies)**second
            total += coefficient*choose(maximum-degree+2*weight, 2*weight)
    return 2**weight*total


def collision_lower(k_log=20, period=16, length=16, dimension=12, q=256,
                    distance=Fraction(1, 10), region_rows=None, width=16):
    if not 0 < dimension < length <= q or period < 1 or period % width:
        raise ValueError('invalid geometry')
    rows = region_rows or (((1 << k_log)+dimension*period-1)//(dimension*period))*period
    if rows % period or rows//period < 2:
        raise ValueError('period must divide regional length with at least two copies')
    copies = rows//period
    n, weight = rows*length, length-dimension+1
    threshold = (distance.numerator*n+distance.denominator-1)//distance.denominator-1
    maximum = (threshold-weight)//period
    gap = Fraction(distance_sum_count(copies, weight, maximum), (copies*(copies-1))**weight)
    pair_label = Fraction(copies-1, rows-1)
    triple_label = Fraction((copies-1)*(copies-2), (rows-1)*(rows-2))
    two_pair_label = pair_label*Fraction((copies-2)*(copies-3)+(period-1)*copies*(copies-1),
                                         (rows-2)*(rows-3))
    diagonal = Fraction(1, (q-1)**(weight-1))
    supports = choose(length, weight)
    single = choose(rows, 2)*supports*pair_label**weight*diagonal*gap
    same_pair, overlapping_pairs, disjoint_pairs = Fraction(0), Fraction(0), Fraction(0)
    for overlap in range(weight+1):
        ordered_support_pairs = supports*choose(weight, overlap)*choose(length-weight, weight-overlap)
        if not ordered_support_pairs:
            continue
        if overlap < weight:
            rank = 2*(weight-1) if overlap == 0 else 2*weight-overlap-1
            same_pair += (Fraction(choose(rows, 2)*ordered_support_pairs, 2)
                          *pair_label**(2*weight-overlap)/((q-1)**rank))
        overlapping_pairs += (3*choose(rows, 3)*ordered_support_pairs
                              *triple_label**overlap*pair_label**(2*(weight-overlap))*diagonal**2)
        disjoint_pairs += (3*choose(rows, 4)*ordered_support_pairs
                           *two_pair_label**overlap*pair_label**(2*(weight-overlap))*diagonal**2)
    second = same_pair+overlapping_pairs+disjoint_pairs
    lower = max(Fraction(0), single-second)
    return dict(k=rows*dimension, n=n, rows=rows, copies=copies, period=period,
                outer=[length, dimension], q=q, distance=str(distance),
                minimum_support=weight, bad_weight_threshold=threshold, maximum_gap_sum=maximum,
                geometry_probability=gap, first_bonferroni=single, same_pair_second=same_pair,
                overlapping_pairs_second=overlapping_pairs, disjoint_pairs_second=disjoint_pairs,
                failure_lower=lower, exceeds_2_to_minus_40=(lower > Fraction(1, 1 << 40)))


def paired_window_upper(k_log=20, period=16, length=32, dimension=24, q=256):
    """Upper bound for the paired-return window family, not all bad words.

    All active coordinates of both rows have matching labels, and both
    words have the same outer support. Every active region returns to zero.
    The existence probability for a fixed support is clipped at one before
    averaging route geometry; this avoids a large message-count loss.
    """
    from q1_batched_screen import fixed_support_multiplicity
    base = collision_lower(k_log, period, length, dimension, q)
    rows, copies, threshold = base['rows'], base['copies'], base['bad_weight_threshold']
    pair_label = Fraction(copies-1, rows-1)
    total, entries = Fraction(0), []
    for weight in range(length-dimension+1, length+1):
        count = fixed_support_multiplicity(length, dimension, weight, q)
        diagonal_upper = min(Fraction(1), Fraction(count**2, (q-1)**(weight+1)))
        maximum = (threshold-weight)//period
        gap = Fraction(distance_sum_count(copies, weight, maximum), (copies*(copies-1))**weight)
        contribution = choose(rows, 2)*choose(length, weight)*pair_label**weight*diagonal_upper*gap
        total += contribution
        entries.append(dict(weight=weight, contribution=contribution))
    return dict(k=base['k'], n=base['n'], rows=rows, outer=[length, dimension], period=period,
                family='matching-label-paired-return-windows', failure_upper=total, entries=entries)


def joint_match_tail(length, minimum, marginal, joint):
    """Exact P(two regional match counts >= minimum), no geometry constraint."""
    denominator = math.lcm(marginal.denominator, joint.denominator)
    both = int(joint*denominator)
    only = int((marginal-joint)*denominator)
    neither = denominator-2*only-both
    current = [[1]]
    for region in range(length):
        following = [[0]*(region+2) for _ in range(region+2)]
        for a, row in enumerate(current):
            for b, value in enumerate(row):
                following[a][b] += value*neither
                following[a+1][b] += value*only
                following[a][b+1] += value*only
                following[a+1][b+1] += value*both
        current = following
    return Fraction(sum(current[a][b] for a in range(minimum, length+1)
                        for b in range(minimum, length+1)), denominator**length)


def kernel_pair_lower(k_log=20, period=16, length=32, dimension=24):
    """Rank-based pair event, independent of diagonals and regional GLs.

    On 2p+1 matched coordinates, two shortened MDS words have dimension
    2(p+1)>2p+1. Therefore a nonzero cancelling pair exists for every D.
    Require exactly this many matching coordinates to avoid counting the
    same row pair once for each support subset. Geometry confines output.
    """
    base = collision_lower(k_log, period, length, dimension)
    rows, copies, threshold = base['rows'], base['copies'], base['bad_weight_threshold']
    minimum = 2*(length-dimension)+1
    if minimum > length:
        raise ValueError('this positive-dimensional paired kernel requires rate above one half')
    pair = Fraction(copies-1, rows-1)
    triple = Fraction((copies-1)*(copies-2), (rows-1)*(rows-2))
    double = pair*Fraction((copies-2)*(copies-3)+(period-1)*copies*(copies-1), (rows-2)*(rows-3))
    gap = Fraction(distance_sum_count(copies, minimum, (threshold-minimum)//period),
                   (copies*(copies-1))**minimum)
    per_pair = choose(length, minimum)*pair**minimum*(1-pair)**(length-minimum)*gap
    first = choose(rows, 2)*per_pair
    shared = 3*choose(rows, 3)*joint_match_tail(length, minimum, pair, triple)
    disjoint = 3*choose(rows, 4)*joint_match_tail(length, minimum, pair, double)
    lower = max(Fraction(0), first-shared-disjoint)
    return dict(k=base['k'], n=base['n'], rows=rows, period=period, outer=[length, dimension],
                matched_coordinates=minimum, kernel_dimension_lower=1, geometry_probability=gap,
                first_bonferroni=first, overlapping_pairs_second=shared, disjoint_pairs_second=disjoint,
                failure_lower=lower, exceeds_2_to_minus_40=(lower > Fraction(1, 1 << 40)))


def report(record):
    result = dict(record)
    for key, value in record.items():
        if isinstance(value, Fraction):
            result[key] = str(value)
            result[key+'_float'] = float(value)
            if value:
                result[key+'_negative_log2'] = math.log2(value.denominator)-math.log2(value.numerator)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', type=int, nargs='+', default=[16, 18, 20])
    parser.add_argument('--period', type=int, nargs='+', default=[16, 256, 1024, 2048])
    args = parser.parse_args()
    for k_log, period in itertools.product(args.k_log, args.period):
        print(json.dumps(report(collision_lower(k_log, period))), flush=True)
