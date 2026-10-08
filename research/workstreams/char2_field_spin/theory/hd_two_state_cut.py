"""Exact finite mixed-interval cut for the two-state regional candidate.

No field arithmetic or Monte Carlo is used. The dynamic program chooses
disjoint intervals in fixed blocks, maximizing 3*gain-cost, where an
interval with m safe and u unsafe positions has gain=m-2 and cost=u+2.
Moment sums count every binary word, then use exact hypergeometric weights.
"""
import argparse
from fractions import Fraction as F
from math import comb
import json


def interval_policy(bits):
    """Return the chosen half-open intervals, using a fixed DP tie rule."""
    outside, inside = (0, []), (-8, [])
    for i, bit in enumerate(bits):
        if bit not in (0, 1):
            raise ValueError("not a bit")
        previous_outside, previous_inside = outside, inside
        if previous_inside[0] > previous_outside[0]:
            outside = previous_inside
        weight = 3 if bit else -1
        if previous_inside[0] > previous_outside[0]-8:
            intervals = previous_inside[1]
            inside = (previous_inside[0]+weight,
                      intervals[:-1]+[(intervals[-1][0], i+1)])
        else:
            inside = (previous_outside[0]-8+weight,
                      previous_outside[1]+[(i, i+1)])
    return (inside if inside[0] > outside[0] else outside)[1]


def interval_statistics(bits, intervals):
    safe = sum(sum(bits[first:last]) for first, last in intervals)
    size = sum(last-first for first, last in intervals)
    count = len(intervals)
    return safe-2*count, size-safe+2*count


def _advance(moment, count, a, keep, shift):
    # Sums of O,D,O^2,D^2,OD. Update O'=O+aD, D'=keep*D+shift.
    o, d, oo, dd, od = moment
    return (o+a*d, keep*d+shift*count,
            oo+2*a*od+a*dd,
            keep*dd+2*keep*shift*d+shift*shift*count,
            keep*od+shift*o+a*shift*d)


def word_moments(length):
    """By Hamming weight: counts and sums G,C,G^2,C^2 over all words.

    D is the integer score difference (inside minus outside). Its finite
    range is [-9,3]. The first and second statistic moments retain exact
    correlations with D; no expectation is passed through a max operation.
    """
    if length < 1:
        raise ValueError("empty block")
    # count, [O,D,O^2,D^2,OD] for gain, then for cost.
    states = {(0, -8): (1, 0, -2, 0, 4, 0, 0, 2, 0, 4, 0)}
    for _ in range(length):
        next_states = {}
        for (ones, difference), values in states.items():
            count = values[0]
            a = int(difference > 0)
            keep = int(-8 < difference <= 0)
            restart = int(difference <= -8)
            clipped = min(0, max(-8, difference))
            for bit in (0, 1):
                key = ones+bit, clipped+(3 if bit else -1)
                gain = _advance(values[1:6], count, a, keep, bit-2*restart)
                cost = _advance(values[6:11], count, a, keep, 1-bit+2*restart)
                updated = (count,)+gain+cost
                previous = next_states.get(key)
                next_states[key] = (updated if previous is None else
                                    tuple(x+y for x, y in zip(previous, updated)))
        states = next_states
    result = [[0]*5 for _ in range(length+1)]
    for (ones, difference), values in states.items():
        a = int(difference > 0)
        count = values[0]
        o, d, oo, dd, od = values[1:6]
        g, gg = o+a*d, oo+2*a*od+a*dd
        o, d, oo, dd, od = values[6:11]
        c, cc = o+a*d, oo+2*a*od+a*dd
        for i, value in enumerate((count, g, c, gg, cc)):
            result[ones][i] += value
    return result


def _choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def regional_moments(rows, safe, length, sums, column):
    """Exact mean/variance of the sum over all blocks of one region."""
    if rows % length or rows < 2*length or not 0 <= safe <= rows:
        raise ValueError("need at least two whole blocks and a valid safe count")
    blocks = rows//length
    denominator = comb(rows, safe)
    first = F(sum(row[column]*_choose(rows-length, safe-j)
                  for j, row in enumerate(sums)), denominator)
    second = F(sum(row[column+2]*_choose(rows-length, safe-j)
                   for j, row in enumerate(sums)), denominator)
    convolution = [0]*(2*length+1)
    for j, row in enumerate(sums):
        for k, other in enumerate(sums):
            convolution[j+k] += row[column]*other[column]
    pair = F(sum(value*_choose(rows-2*length, safe-j)
                 for j, value in enumerate(convolution)), denominator)
    mean = blocks*first
    variance = blocks*second+blocks*(blocks-1)*pair-mean*mean
    if variance < 0:
        raise ArithmeticError("negative exact variance")
    return mean, variance


def finite_bound(rows=32768, b=64, k=32, a=8192, pre_order=2,
                 block=256, distance=F(2, 5)):
    if b != 2*k or not 0 <= pre_order <= a <= rows:
        raise ValueError("invalid outer geometry or source cut")
    if not 0 < distance < 1:
        raise ValueError("invalid distance")
    n, h, safe = b*rows, k*a, a-pre_order
    target = distance*n
    e = (target.numerator+target.denominator-1)//target.denominator-1
    cap = e-h-1
    sums = word_moments(block)
    mean_g, var_g = regional_moments(rows, safe, block, sums, 1)
    mean_c, var_c = regional_moments(rows, safe, block, sums, 2)
    mean_g, var_g, mean_c, var_c = b*mean_g, b*var_g, b*mean_c, b*var_c
    gap_g, gap_c = mean_g-h, cap-mean_c
    if gap_g <= 0 or gap_c <= 0:
        raise ValueError("this interval policy has no positive mean gaps")
    bad_g, bad_c = var_g/gap_g**2, var_c/gap_c**2
    lower = max(F(0), 1-bad_g-bad_c)
    return dict(rows=rows, b=b, k=k, n=n, a=a, safe=safe, block=block,
                check_capacity=h, forbidden_cutoff=e, cost_cap=cap,
                mean_gain=mean_g, mean_cost=mean_c,
                variance_gain=var_g, variance_cost=var_c,
                gap_gain=gap_g, gap_cost=gap_c,
                bad_gain_upper=bad_g, bad_cost_upper=bad_c,
                failure_lower=lower, exceeds_19_over_20=lower > F(19, 20),
                rigorous_arithmetic=True,
                bound_type="whole-code failure lower bound",
                routing="independent uniform permutations inside regions")


def direct_window_bound(states, rows=32768, b=64, k=32, distance=F(2, 5)):
    """Actual obstruction for no precode and an r-state causal postscan.

    Shorten outer rows to a prefix of each routed region, then impose zero
    postscan state at each prefix end. Only the expectation of the total
    shortened dimension is used; distinct outer rows need not be independent.
    """
    if not 0 < k < b or states < 1 or rows < 1 or not 0 < distance < 1:
        raise ValueError("invalid direct-window geometry")
    target = distance*rows
    width = (target.numerator+target.denominator-1)//target.denominator-1
    probability = F(width, rows)
    mean_dimension = rows*sum(F(comb(b, t))*probability**t*(1-probability)**(b-t)
                              *(t-(b-k)) for t in range(b-k+1, b+1))
    constraints, maximum = b*states, k*rows
    lower = (max(F(0), (mean_dimension-constraints)/(maximum-constraints))
             if constraints < maximum else F(0))
    return dict(states=states, width=width, output_weight_bound=b*width,
                mean_shortened_dimension=mean_dimension,
                state_constraints=constraints, failure_lower=lower,
                exceeds_2_to_minus_40=lower > F(1, 1 << 40),
                bound_type="whole-code failure lower bound; no precode")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block", type=int, default=256)
    parser.add_argument("--exact", action="store_true",
                        help="also print potentially long exact fractions")
    args = parser.parse_args()
    result = finite_bound(block=args.block)
    output = {name: ({"display": float(value),
                      **({"exact": str(value)} if args.exact else {})}
                     if isinstance(value, F) else value)
              for name, value in result.items()}
    print(json.dumps(output, indent=2))
