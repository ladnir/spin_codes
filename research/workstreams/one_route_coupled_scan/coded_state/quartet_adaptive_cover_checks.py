"""Exact tiny checks for QUARTET_ADAPTIVE_HALL_NORMALIZATION.md.

Reverse state: A=V_(t+1), W=V_(t+1) intersect V_(t+2), W subset A.
For a quartet put u=|A|, w=|W|, j=|J|, r=|J intersect W|. The
selected-quartet score increment lies in [b-4-u,b], b=2j-r. Set
a=max(0,4+u-b). Adaptive asymmetric weights give hole probability
 g=y^-b (1-y^-a)/(1-y^-(a+b))
and selection probability p=1-g. j=0 forces a hole; a=0 forces selection.
These probabilities MUST be recomputed from the actual future branch.

The tests compare the reverse decision tree with physical Hall scores,
including exact balance, unbalanced offsets, forced input-1 selection,
formal full future holes, and all15 local coefficient types. No large
operator, scalar optimizer, or numerical experiment is implemented.

An independent integer check certifies a limitation of the proposed
coefficientwise maximum across15 types at y=2. It concerns that scalar
comparison only, not the adaptive cover moment or the construction.
No jobs or large arrays are created on import.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import comb
import os
import shutil
import threading
import time


def _ac_choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(maxsize=None)
def adaptive_cover_weights(y, u, w, j, r):
    """Return (hole probability, selection probability, absent weight,a,b)."""
    y = Fraction(y)
    if y <= 1 or not 0 <= w <= u <= 4 or not 0 <= j <= 4:
        raise ValueError('Invalid marker or quartet counts')
    if not 0 <= r <= min(w, j) or j-r > 4-w:
        raise ValueError('Invalid J/W overlap')
    b, a = 2*j-r, max(0, 4+u-(2*j-r))
    if j == 0:
        return Fraction(1), Fraction(0), Fraction(1), a, b
    if a == 0:
        return Fraction(0), Fraction(1), Fraction(0), a, b
    A, B = y**(-a), y**(-b)
    absent, present = (1-A)/(1-A*B), (1-B)/(1-A*B)
    return B*absent, present, absent, a, b


def adaptive_cover_local_coefficients(y, alpha, beta, u, w):
    """Fixed-old-type local action divided by alpha^u beta^w, in J degree."""
    y, alpha, beta = map(Fraction, (y, alpha, beta))
    if alpha <= 0 or beta <= 0 or not 0 <= w <= u <= 4:
        raise ValueError('Use positive potentials and a valid nested type')
    coefficients = []
    for j in range(5):
        value = Fraction(0)
        for r in range(max(0, j-(4-w)), min(w, j)+1):
            hole, present, _, _, _ = adaptive_cover_weights(y, u, w, j, r)
            value += (comb(w, r)*comb(4-w, j-r)*y**(w-r)
                      *(present+hole*alpha**4*beta**u))
        coefficients.append(value/(alpha**u*beta**w))
    return tuple(coefficients)


def _ac_geometry(partitions, width):
    masks, lookup = [], []
    for partition in partitions:
        assert sorted(coordinate for quartet in partition for coordinate in quartet) == list(range(width))
        assert all(len(quartet) == 4 for quartet in partition)
        row, reverse = [], [0]*width
        for quartet in partition:
            mask = sum(1 << coordinate for coordinate in quartet)
            row.append(mask)
            for coordinate in quartet:
                reverse[coordinate] = mask
        masks.append(tuple(row))
        lookup.append(tuple(reverse))
    return tuple(masks), tuple(lookup)


def _ac_predecessors(C, lookup):
    """Independent direct own/next/two-ahead physical input predecessor union."""
    result = list(C)
    for t, current in enumerate(C):
        while current:
            bit = current & -current
            quartet = lookup[t][bit.bit_length()-1]
            if t >= 1:
                result[t-1] |= quartet
            if t >= 2:
                result[t-2] |= quartet
            current -= bit
    return tuple(result)


def _ac_input_subsets(J):
    positions = [(t, 1 << coordinate) for t, mask in enumerate(J)
                 for coordinate in range(mask.bit_length()) if mask >> coordinate & 1]
    current = [0]*len(J)
    yield tuple(current)
    for counter in range(1, 1 << len(positions)):
        t, bit = positions[(counter & -counter).bit_length()-1]
        current[t] ^= bit
        yield tuple(current)


def _ac_instance(J, partitions, width, counters):
    full, updates = (1 << width)-1, len(J)
    assert updates == len(partitions) and all(0 <= mask <= full for mask in J)
    quartets, lookup = _ac_geometry(partitions, width)
    candidates = [(t, quartet) for t in range(1, updates)
                  for quartet in quartets[t] if quartet & J[t]]
    assert len(candidates) <= 12
    all_records, scores = [], []
    for bits in range(1 << len(candidates)):
        U, C = [0]*updates, [J[0]]+[0]*(updates-1)
        for i, (t, quartet) in enumerate(candidates):
            if bits >> i & 1:
                U[t] |= quartet
                C[t] |= J[t] & quartet
        V = [full ^ mask for mask in U]+[full, full]
        predecessor = _ac_predecessors(C, lookup)
        score = 2*sum(mask.bit_count() for mask in C)-sum(mask.bit_count() for mask in predecessor)
        future = [(V[t+1], V[t+1] & V[t+2]) for t in range(updates)]
        cover_cost = sum(((full ^ J[t]) & future[t][1]).bit_count() for t in range(updates))
        records, absent_b = [], 0
        for t in range(updates-1, 0, -1):
            A, W = future[t]
            for quartet in quartets[t]:
                u, w = (A & quartet).bit_count(), (W & quartet).bit_count()
                j, r = (J[t] & quartet).bit_count(), (J[t] & W & quartet).bit_count()
                selected = bool(U[t] & quartet)
                assert not selected or j > 0
                records.append((u, w, j, r, selected))
                if not selected:
                    absent_b += 2*j-r
        imbalance = 2*sum(mask.bit_count() for mask in J)-width*updates
        assert score+absent_b == cover_cost+imbalance
        assert C[0] == J[0]  # No first-layer cover coin or selectable U1.
        all_records.append((records, cover_cost, future))
        scores.append(score)
        counters['physical_branch_identities'] += 1
    # The exact input-subset maximum independently checks canonical saturation.
    maximum = -width*updates
    for C in _ac_input_subsets(J):
        score = 2*sum(mask.bit_count() for mask in C)-sum(mask.bit_count() for mask in _ac_predecessors(C, lookup))
        maximum = max(maximum, score)
        counters['physical_input_subsets'] += 1
    assert maximum == max(scores)
    counters['canonical_maxima'] += 1
    for bits, score in enumerate(scores):
        future = all_records[bits][2]
        for i, (t, quartet) in enumerate(candidates):
            if bits >> i & 1:
                continue
            A, W = future[t]
            u = (A & quartet).bit_count()
            j, r = (J[t] & quartet).bit_count(), (J[t] & W & quartet).bit_count()
            b = 2*j-r
            difference = scores[bits | (1 << i)]-score
            assert b-4-u <= difference <= b
            assert difference >= -max(0, 4+u-b)
            counters['adaptive_physical_flip_bounds'] += 1
    for y in (Fraction(2), Fraction(4)):
        total_probability = weighted = cover_moment = Fraction(0)
        for score, (records, cover_cost, _) in zip(scores, all_records):
            probability = weight = Fraction(1)
            for u, w, j, r, selected in records:
                hole, present, absent, a, b = adaptive_cover_weights(y, u, w, j, r)
                assert hole+present == 1 and absent+y**(-a)*present == 1
                assert y**(-b)*absent+present == 1
                probability *= present if selected else hole
                weight *= present if selected else absent
            total_probability += probability
            weighted += weight*y**score
            cover_moment += probability*y**cover_cost
        assert total_probability == 1
        assert y**maximum <= weighted
        assert weighted == y**imbalance*cover_moment
        counters['adaptive_probability_normalizations'] += 1
        counters['adaptive_decision_tree_majorants'] += 1
        if imbalance:
            assert weighted != cover_moment
            counters['unbalanced_offset_controls'] += 1
        else:
            assert weighted == cover_moment
            counters['balanced_cover_identities'] += 1


def _ac_poly_product(left, right):
    result = [Fraction(0)]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return result


def _ac_poly_power(poly, exponent):
    result = [Fraction(1)]
    for _ in range(exponent):
        result = _ac_poly_product(result, poly)
    return result


def _ac_local_checks(counters):
    for y in (Fraction(2), Fraction(4)):
        for alpha, beta in ((Fraction(1), Fraction(1)), (Fraction(3, 2), Fraction(2, 3)),
                            (Fraction(2, 3), Fraction(3, 2))):
            for u in range(5):
                for w in range(u+1):
                    formula = adaptive_cover_local_coefficients(y, alpha, beta, u, w)
                    # Fixed labeled nested masks; no old-state multiplicity enters a row.
                    A, W = (1 << u)-1, (1 << w)-1
                    physical = [Fraction(0)]*5
                    for J in range(16):
                        j, r = J.bit_count(), (J & W).bit_count()
                        hole, present, _, _, _ = adaptive_cover_weights(y, u, w, j, r)
                        for V, probability in ((0, present), (15, hole)):
                            next_A, next_W = V, V & A
                            physical[j] += (probability*y**((15 ^ J) & W).bit_count()
                                            *alpha**next_A.bit_count()*beta**next_W.bit_count()
                                            /(alpha**u*beta**w))
                    assert tuple(physical) == formula
                    counters['local_type_coefficients'] += 5
                    # Forced input-1 endpoint keeps all labels and no adaptive coin.
                    endpoint = [Fraction(0)]*5
                    for J in range(16):
                        endpoint[J.bit_count()] += y**((15 ^ J) & W).bit_count()/(alpha**u*beta**w)
                    exact = [value/(alpha**u*beta**w) for value in
                             _ac_poly_product(_ac_poly_power((y, 1), w), _ac_poly_power((1, 1), 4-w))]
                    assert endpoint == exact
                    E = (max(Fraction(1), 1/alpha, y/(alpha*beta)),
                         max(Fraction(1), 1/alpha, 1/(alpha*beta)))
                    bound = _ac_poly_power(E, 4)
                    assert all(value <= upper for value, upper in zip(endpoint, bound))
                    counters['forced_endpoint_coefficients'] += 5
            lam = [max(adaptive_cover_local_coefficients(y, alpha, beta, u, w)[j]
                       for u in range(5) for w in range(u+1)) for j in range(5)]
            # A width8 partition has two independently labeled old quartet types.
            for first, second in (((0, 0), (4, 4)), ((3, 1), (2, 2)), ((4, 0), (1, 1))):
                actual = _ac_poly_product(adaptive_cover_local_coefficients(y, alpha, beta, *first),
                                          adaptive_cover_local_coefficients(y, alpha, beta, *second))
                comparison = _ac_poly_power(lam, 2)
                assert all(a <= b for a, b in zip(actual, comparison))
                counters['coefficientwise_partition_comparisons'] += 9


def _ac_log2_ratio_floor(numerator, denominator):
    """Exact floor(log2(numerator/denominator)), without floating logarithms."""
    if numerator <= 0 or denominator <= 0:
        raise ValueError('Use positive integers')
    exponent = numerator.bit_length()-denominator.bit_length()
    below = numerator < denominator << exponent if exponent >= 0 else numerator << (-exponent) < denominator
    return exponent-int(below)


def adaptive_cover_scalar_obstruction():
    """Exact y2 obstruction for the coefficient-max scalar family, all alpha,beta.

    Full/full old type and j1 has a potential-cancelling hole monomial
    4*2^3*(127/255). Empty/empty old type j2,j3 forces selection and gives
    6,4. Thus lambda(x)>=L(x)=(4064/255)x+6x^2+4x^3 coefficientwise.
    The endpoint satisfies E(x)>=(2+x)/(alpha*beta), cancelling the
    initial potential. Single multinomial terms suffice below.
    """
    hole, _, _, a, b = adaptive_cover_weights(Fraction(2), 4, 4, 1, 1)
    assert (a, b, hole) == (7, 1, Fraction(127, 255))
    for j in (2, 3):
        hole, present, _, a, _ = adaptive_cover_weights(Fraction(2), 0, 0, j, 0)
        assert a == 0 and hole == 0 and present == 1
    for alpha, beta in ((Fraction(1, 8), Fraction(4)), (Fraction(2), Fraction(1, 3)),
                        (Fraction(16), Fraction(32))):
        assert adaptive_cover_local_coefficients(2, alpha, beta, 4, 4)[1] >= Fraction(4064, 255)
        empty = adaptive_cover_local_coefficients(2, alpha, beta, 0, 0)
        assert empty[2] >= 6 and empty[3] >= 4
    # Degree683+2*682+3*683=4096 in2048 quartic factors.
    bulk_numerator = comb(2048, 683)*comb(1365, 682)*4064**683*6**682*4**683
    bulk_denominator = 255**683
    # Degree640+2*640+3*640+256=4096 in1920 factors plus512 labels.
    endpoint_numerator = (comb(1920, 640)*comb(1280, 640)*4064**640*6**640*4**640
                          *comb(512, 256)*2**256)
    endpoint_denominator = 255**640
    bulk_bits = _ac_log2_ratio_floor(bulk_numerator, bulk_denominator)
    endpoint_bits = _ac_log2_ratio_floor(endpoint_numerator, endpoint_denominator)
    quota_bits_upper = comb(8192, 4096).bit_length()
    exponent = 255*bulk_bits+endpoint_bits-255*quota_bits_upper-104860
    assert exponent > 0
    return dict(marker_y=2, universal_bulk_log2_floor=bulk_bits,
                initial_times_endpoint_log2_floor=endpoint_bits,
                regional_binomial_log2_upper=quota_bits_upper,
                scalar_comparison_tail_log2_lower=exponent,
                exact_integer_arithmetic=True, all_positive_alpha_beta=True,
                only_coefficient_max_scalar_family=True,
                adaptive_cover_moment_obstruction=False,
                actual_failure_probability_lower_bound=False)


def adaptive_cover_checks():
    counters = {name: 0 for name in (
        'physical_branch_identities', 'physical_input_subsets', 'canonical_maxima',
        'adaptive_physical_flip_bounds', 'adaptive_probability_normalizations',
        'adaptive_decision_tree_majorants', 'unbalanced_offset_controls',
        'balanced_cover_identities', 'local_type_coefficients',
        'forced_endpoint_coefficients', 'coefficientwise_partition_comparisons')}
    single = ((0, 1, 2, 3),)
    for updates in (1, 2, 3):
        for positions in combinations(range(4*updates), 2*updates):
            J = [0]*updates
            for position in positions:
                J[position//4] |= 1 << (position % 4)
            _ac_instance(tuple(J), (single,)*updates, 4, counters)
    for J in ((0, 0), (15, 15), (1, 0, 3), (15, 7, 15)):
        _ac_instance(J, (single,)*len(J), 4, counters)
    partitions = (((0, 1, 2, 3), (4, 5, 6, 7)),
                  ((0, 2, 4, 6), (1, 3, 5, 7)),
                  ((0, 3, 4, 7), (1, 2, 5, 6)),
                  ((0, 1, 4, 5), (2, 3, 6, 7)))
    # Each pair of updates has exactly eight non-H slots; partitions vary.
    for J in ((15, 240, 85, 170), (51, 204, 102, 153), (0, 255, 255, 0)):
        _ac_instance(J, partitions, 8, counters)
    _ac_local_checks(counters)
    for numerator, denominator, expected in ((1, 1, 0), (3, 2, 0), (1, 3, -2),
                                             (8, 1, 3), (7, 8, -1), (256, 255, 0)):
        assert _ac_log2_ratio_floor(numerator, denominator) == expected
    print(dict(**counters, exact_log2_floor_checks=6,
               scalar_family_obstruction=adaptive_cover_scalar_obstruction(),
               full_code_certificate=False), flush=True)


def adaptive_cover_checks_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks', action='store_true', required=True)
    parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        adaptive_cover_checks()
        print(dict(seconds=round(time.monotonic()-started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    adaptive_cover_checks_main()
