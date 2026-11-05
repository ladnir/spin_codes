"""Exact common-source route partition; no distance bound or production gate.

A routed source pair spans emission blocks (0,1) or (2,3), not two
coordinates of one block. Each coupled half contains block_length pair
slots. The scalar good-event probability supports production parameters.
The six-type matrix builder is intentionally restricted to tiny fixtures;
it specifies the conditional operator interface without an exponential
production implementation. All coefficients are nonnegative rationals.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial


ZERO = (F(0),) * 4
IDENTITY = (F(1), F(0), F(0), F(1))
TYPE_NAMES = ('outside', 'single_inactive', 'single_active',
              'full_inactive', 'full_one_active', 'full_both_active')


def _dimensions(block_length, packets, cap):
    if (any(type(x) is not int for x in (block_length, packets, cap))
            or block_length < 1 or packets < 1 or not 0 <= cap <= block_length):
        raise ValueError('Positive blocks/packets and cap in 0..block_length required')


def coupled_half_polynomial(*, block_length=16, cap=4):
    """Integer coefficients of H(X,Y), before the singleton 2^-s factor."""
    _dimensions(block_length, 1, cap)
    result = {}
    for full in range(cap+1):
        for left in range(cap-full+1):
            for right in range(cap-full+1):
                if left+right+full > block_length:
                    continue
                count = (comb(block_length, full)*comb(block_length-full, left)
                         *comb(block_length-full-left, right))
                key = left+right, full
                result[key] = result.get(key, 0)+count
    return result


def _scalar_product(left, right, target):
    result = {}
    for (s, f), a in left.items():
        for (t, g), b in right.items():
            key = s+t, f+g
            if key[0] <= target[0] and key[1] <= target[1]:
                result[key] = result.get(key, 0)+a*b
    return result


def good_probability(singletons, full_pairs, *, block_length=16, packets=128, cap=4):
    """Exact probability every potential-G emission-block occupancy <=cap.

    Odd/even cyclic source matchings change only (singletons,full_pairs).
    Full-pair swaps are invisible to this event; singleton swaps contribute
    the exact 2^singletons denominator. No packet independence is assumed:
    coefficient extraction conditions on the fixed global source histogram.
    """
    _dimensions(block_length, packets, cap)
    pair_count = 2*block_length*packets
    if (any(type(x) is not int for x in (singletons, full_pairs))
            or min(singletons, full_pairs) < 0 or singletons+full_pairs > pair_count):
        raise ValueError('Invalid fixed source-pair histogram')
    target = singletons, full_pairs
    base = {key: value for key, value in coupled_half_polynomial(
        block_length=block_length, cap=cap).items()
            if key[0] <= singletons and key[1] <= full_pairs}
    total, exponent = {(0, 0): 1}, 2*packets
    while exponent:
        if exponent & 1:
            total = _scalar_product(total, base, target)
        exponent >>= 1
        if exponent:
            base = _scalar_product(base, base, target)
    denominator = 2**singletons*comb(pair_count, singletons)*comb(pair_count-singletons, full_pairs)
    return F(total.get(target, 0), denominator)


def six_type_slot_terms():
    """(type index, potential bits, actual bits, orientation probability).

    Both masks are two-bit masks in the destination coupled half. The
    potential mask belongs to the SAME fixed G for all outer mask terms.
    """
    half = F(1, 2)
    return ((0, 0, 0, F(1)),
            (1, 1, 0, half), (1, 2, 0, half),
            (2, 1, 1, half), (2, 2, 2, half),
            (3, 3, 0, F(1)),
            (4, 3, 1, half), (4, 3, 2, half),
            (5, 3, 3, F(1)))


def matrix_product(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h


def matrix_add(left, right):
    return tuple(a+b for a, b in zip(left, right))


def _matrix_polynomial_product(left, right, target):
    result = {}
    for degrees, a in left.items():
        for other, b in right.items():
            key = tuple(x+y for x, y in zip(degrees, other))
            if all(x <= y for x, y in zip(key, target)):
                result[key] = matrix_add(result.get(key, ZERO), matrix_product(a, b))
    return result


def tiny_packet_polynomials(caps, target, *, block_length=2, cap=1):
    """Good/bad six-degree packet polynomials, with raw ordered matrix caps.

    This reference enumerates at most 9^4 oriented slot outcomes. A future
    production builder must compress it, not lift this tiny guard. There
    is no alpha here and the caps need not commute or be probabilities.
    """
    _dimensions(block_length, 1, cap)
    if block_length > 2 or len(target) != 6 or any(type(x) is not int or x < 0 for x in target):
        raise ValueError('Tiny block length <=2 and six nonnegative counts required')
    good, bad = {}, {}
    for slots in product(six_type_slot_terms(), repeat=2*block_length):
        counts, potential, actual, mass = [0]*6, [0]*4, [0]*4, F(1)
        for slot, (kind, potential_bits, actual_bits, probability) in enumerate(slots):
            counts[kind] += 1
            mass *= probability
            # Native pair index 2*J+half feeds lanes 2*half and 2*half+1.
            base = 2*(slot % 2)
            for lane in range(2):
                potential[base+lane] += (potential_bits >> lane) & 1
                actual[base+lane] += (actual_bits >> lane) & 1
        key = tuple(counts)
        if any(x > y for x, y in zip(key, target)):
            continue
        matrix = tuple(F(value) for value in caps[tuple(actual)])
        if len(matrix) != 4 or any(value < 0 for value in matrix):
            raise ValueError('Four nonnegative unpowered cap entries required')
        destination = good if max(potential) <= cap else bad
        destination[key] = matrix_add(destination.get(key, ZERO), tuple(mass*x for x in matrix))
    return good, bad


def tiny_regional_partition(caps, target, *, block_length=2, packets=2, cap=1):
    """Unconditional event-weighted and conditional cap means for fixed U⊆G.

    target counts the six joint (G,U) pair types. The positive bad-flag
    recurrence preserves packet chronology and performs no subtraction.
    Zero-probability conditional branches are returned as None.
    """
    _dimensions(block_length, packets, cap)
    if len(target) != 6 or any(type(x) is not int or x < 0 for x in target):
        raise ValueError('Six nonnegative joint source/active-pair counts required')
    if sum(target) != 2*block_length*packets:
        raise ValueError('Joint histogram must fill every pair slot')
    good_packet, bad_packet = tiny_packet_polynomials(caps, target, block_length=block_length, cap=cap)
    all_packet = dict(good_packet)
    for key, value in bad_packet.items():
        all_packet[key] = matrix_add(all_packet.get(key, ZERO), value)
    good, bad = {(0,)*6: IDENTITY}, {}
    for _ in range(packets):
        next_bad = _matrix_polynomial_product(bad, all_packet, target)
        newly_bad = _matrix_polynomial_product(good, bad_packet, target)
        for key, value in newly_bad.items():
            next_bad[key] = matrix_add(next_bad.get(key, ZERO), value)
        good = _matrix_polynomial_product(good, good_packet, target)
        bad = next_bad
    denominator = factorial(sum(target))
    for count in target:
        denominator //= factorial(count)
    mean_good = tuple(x/denominator for x in good.get(tuple(target), ZERO))
    mean_bad = tuple(x/denominator for x in bad.get(tuple(target), ZERO))
    probability = good_probability(target[1]+target[2], sum(target[3:]),
                                   block_length=block_length, packets=packets, cap=cap)
    return dict(good_probability=probability, bad_probability=1-probability,
                good_weighted=mean_good, bad_weighted=mean_bad,
                good_conditional=tuple(x/probability for x in mean_good) if probability else None,
                bad_conditional=tuple(x/(1-probability) for x in mean_bad) if probability != 1 else None,
                alpha_applied=False, distance_certificate=False)
