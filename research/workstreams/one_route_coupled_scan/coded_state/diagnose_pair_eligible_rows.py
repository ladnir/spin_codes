"""RS-zero-aware clean-episode diagnostic; fixed research construction.

Exact rational route moments and integer upper bounds. No generated files.
Restricted family only: disjoint crossing episodes, followed by a live tail.
"""
import argparse
import itertools
import math
import time
from fractions import Fraction as F

from diagnose_pair_region_rank import count_cap, upper_log_fraction, SCALE

M = 8192


def elementary(weights):
    out = [F(1)]
    for weight in weights:
        out.append(F(0))
        for j in range(len(out)-1, 0, -1):
            out[j] += weight*out[j-1]
    return out


def eligible_majorant(lengths, dimension):
    """A row needs >=129 routed coordinates inside the allowed intervals.

    Give the tail full occupancy. Among 2m private early regions it must
    occupy at least 2m-(dimension-1). Equalizing the two fragment fractions
    of an episode increases every elementary symmetric polynomial.
    """
    weights = [F(length+1, 32) for length in lengths for _ in range(2)]
    need = len(weights)-dimension+1
    if need <= 0:
        return F(1)
    coefficients = elementary(weights)
    return min([F(1)] + [coefficients[j]/math.comb(need, j)
                          for j in range(1, need+1)])


def reset_prices(lengths, dimension):
    # w_i >= (n_left+n_right)/(length+1). An active RS row
    # can omit <=dimension-1 of the touched regions. Removing the largest
    # prices bounds its remaining reset charge from below.
    prices = sorted([(32-(length-1).bit_length())//(length+1)
                     for length in lengths for _ in range(2)], reverse=True)
    linear = sum(prices[dimension-1:])
    absolute = sum((32-(length-1).bit_length())*(1 << length) for length in lengths)
    return linear, absolute


def case(lengths, dimension, bits):
    # Retain every available shortened-dimension threshold jointly. Actual
    # words may have extra zeros: their weighted enumerator contributes a
    # small correction below; available dimension does not fix word weight.
    caps = [count_cap(eligible_majorant(lengths, dimension-j+1),
                      bits+(dimension-1).bit_length())
            for j in range(1,dimension+1)]
    assert all(x >= y for x,y in zip(caps,caps[1:]))
    cap = caps[0]
    minimum = 1 << max(lengths)
    if cap < minimum:
        return dict(upper=-(1 << 60), cap=cap, minimum=minimum, excluded=True)
    prices = [reset_prices(lengths, dimension-j+1)[0]
              for j in range(1,dimension+1)]
    _, absolute = reset_prices(lengths, dimension)
    assert all(0 <= y-x < 32 for x,y in zip(prices,prices[1:]))
    beta_price = upper_log_fraction(eligible_majorant(lengths, dimension))
    terms = []
    for active in range(minimum, cap+1):
        levels = [min(active,c) for c in caps] + [0]
        message_dimension = sum(levels)
        linear = sum((levels[j]-levels[j+1])*prices[j] for j in range(dimension))
        messages = (math.comb(cap, active)-1).bit_length() + 32*message_dimension
        # Each extra zero costs a message degree (32 bits), has <=256
        # possible locations, and saves <=16 reset bits. The weighted
        # nonzero row sum is <=(256/255)*q^d*2^-P_d. Since
        # (256/255)^128<2, all s rows cost at most ceil(s/128) bits.
        correction = (active+127)//128
        conditional = messages-max(linear-correction, absolute)
        # A second bound keeps route rarity in the first moment. It uses
        # all M rows, not the good-route cap, to avoid double counting it.
        route_price = -((-active*beta_price)//SCALE)
        unconditional = ((math.comb(M,active)-1).bit_length()
                         + 32*dimension*active + route_price
                         - max(prices[0]*active, absolute))
        terms.append(min(conditional, unconditional))
    # One factor two for the accepted pair-map density, then the row-count union.
    upper = 1 + max(terms) + (len(terms)-1).bit_length()
    return dict(upper=upper, cap=cap, minimum=minimum, excluded=False,
                caps=caps, prices=prices, absolute=absolute)


def run(max_dimension=3, seconds=35):
    began = time.monotonic()
    closed = 5
    placement = 26*(closed+1)
    # <=792 length multisets per dimension. Endpoints already determine the
    # ordered lengths, so only a maximum is needed inside the placement union.
    bits = placement+100
    for dimension in range(1, max_dimension+1):
        cases = excluded = failed = 0
        worst = -(1 << 60)
        witness = None
        for lengths in itertools.combinations_with_replacement(range(1,9), closed):
            result = case(lengths, dimension, bits)
            cases += 1
            excluded += result['excluded']
            failed += result['upper']+placement > -80
            if result['upper'] > worst:
                worst, witness = result['upper'], (lengths, result)
            if time.monotonic()-began > seconds:
                raise TimeoutError('eligible-row diagnostic cap; incomplete')
        good = worst+placement
        print(dict(regions=128+dimension, dimension=dimension, cases=cases,
                   excluded=excluded, cells_above_minus80=failed, good_route_upper=good,
                   route_exception_upper=-100, family_upper=max(good,-100)+1,
                   worst=witness, whole_code_certificate=False,
                   seconds=time.monotonic()-began), flush=True)


def test():
    began = time.monotonic()
    assert 256**128 < 2*255**128
    majorants = prices = balancing = nested = 0
    for weights in itertools.product((F(1,4), F(1,2), F(3,4)), repeat=4):
        coefficients = elementary(weights)
        for need in range(1,5):
            actual = sum(math.prod(weights[i] if mask>>i&1 else 1-weights[i]
                                   for i in range(4))
                         for mask in range(16) if mask.bit_count() >= need)
            for degree in range(1,need+1):
                assert actual <= coefficients[degree]/math.comb(need,degree)
                majorants += 1
    for lengths in itertools.product(range(1,5), repeat=2):
        base = [(32-(l-1).bit_length())//(l+1) for l in lengths for _ in range(2)]
        for dimension in range(1,5):
            linear, _ = reset_prices(lengths, dimension)
            for mask in range(16):
                if 4-mask.bit_count() <= dimension-1:
                    assert sum(base[j] for j in range(4) if mask>>j&1) >= linear
                    prices += 1
        for split1 in range(1,lengths[0]+1):
            for split2 in range(1,lengths[1]+1):
                actual = elementary([F(split1,16), F(lengths[0]+1-split1,16),
                                     F(split2,16), F(lengths[1]+1-split2,16)])
                balanced = elementary([F(l+1,32) for l in lengths for _ in range(2)])
                assert all(x <= y for x,y in zip(actual,balanced))
                balancing += 1
    # For fixed caps, filling every allowed dimension level maximizes both
    # entropy-minus-linear-reset and entropy-minus-absolute-reset. Check
    # their minimum against every small degree assignment, not only bounds
    # on their two separate pieces.
    for active in range(1,5):
        for cap2 in range(active+1):
            for cap3 in range(cap2+1):
                levels = [active,cap2,cap3,0]
                p = [8,16,24]
                saturated = 32*sum(levels)-max(100, sum(
                    (levels[j]-levels[j+1])*p[j] for j in range(3)))
                for degrees in itertools.product(range(1,4),repeat=active):
                    if sum(d>=2 for d in degrees)>cap2 or sum(d>=3 for d in degrees)>cap3:
                        continue
                    actual = 32*sum(degrees)-max(100,sum(p[d-1] for d in degrees))
                    assert actual <= saturated
                    nested += 1
    print(dict(factorial_majorant_checks=majorants, reset_zero_checks=prices,
               balancing_checks=balancing, nested_dimension_checks=nested,
               seconds=time.monotonic()-began), flush=True)
    assert time.monotonic()-began < 10


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--max-dimension', type=int, default=3)
    args = parser.parse_args()
    assert 1 <= args.max_dimension <= 6
    if args.test:
        test()
    else:
        run(args.max_dimension)
