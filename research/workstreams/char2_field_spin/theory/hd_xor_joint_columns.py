"""Joint XOR2 column-category bounds with an aggregate RS support floor.

Each active group has two nonzero RS[32,16] rows. The exact column-category
contractions are retained. An information-set count bounds prescribed RS
supports; only the sum of the constituent support floors is relaxed.
Directed replay checks every labelled-region multiset case before evaluating
a positive univariate polynomial tail. No whole-code certificate is claimed.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement, islice
from math import ceil, comb, log2
from time import perf_counter

import numpy as np

from hd_xor_column_potential import (FULL, PAIRS, LIMIT, Q, REGION_LENGTH,
    arb_region_matrices, factored_column_vectors, factored_word_vector,
    chronological_words, mask_code, multiset_code_batches, streaming_multiset_screen)
from hd_xor_multigroup import region_operators

# The packet profile formula covers every degree through 64 and regional
# normalization covers 0..8192. Thirty-two is a bounded evaluator limit,
# not a change to the occupancy or support-count proof.
MAX_ENVELOPE_GROUPS = 32


def categories(groups):
    return tuple((single, double) for single in range(groups+1)
                 for double in range(groups-single+1))


def joint_multiset_code_batches(groups, single, double, batch_size=1024):
    """Keep the existing r<=5 path; extend the same bounded generator to r=6."""
    if groups <= 5:
        yield from multiset_code_batches(groups, single, double, batch_size)
        return
    if groups != 6 or min(single, double) < 0 or single+double > groups or batch_size < 1:
        raise ValueError('Invalid six-group category or batch size')
    base = groups+1
    full = mask_code(FULL, base)
    pair_shifts = np.array([mask_code(pair, base)-full for pair in PAIRS], dtype=np.int32)
    flags = np.array([[branch >> index & 1 for branch in range(1 << double)]
                      for index in range(double)], dtype=np.int32).reshape(double, 1 << double)
    for rows in combinations_with_replacement(range(8), single):
        initial = (single+double)*full-sum(base**row for row in rows)
        edges = combinations_with_replacement(range(28), double)
        while True:
            batch = list(islice(edges, batch_size))
            if not batch:
                break
            indices = np.array(batch, dtype=np.int32).reshape(len(batch), double)
            yield initial+pair_shifts[indices] @ flags


def six_group_category_screen(z, log_potential):
    """One streaming FLOAT r=6 screen, with at most 1024*64 branches live."""
    if not 0 < z < 1 or log_potential > 0:
        raise ValueError('Invalid weight mark or potential')
    start = perf_counter()
    prefixes = chronological_words(region_operators(z, 6), 4)
    suffixes = np.logaddexp2(prefixes[:, :, 0], prefixes[:, :, 1]+log_potential)
    split = len(prefixes)
    costs, checked, branches = {}, 0, 0
    for single, double in categories(6):
        weights = np.array([branch.bit_count()*log2(1/(Q-1))
                            +(double-branch.bit_count())*log2((Q-2)/(Q-1))
                            for branch in range(1 << double)])
        maximum = -np.inf
        for codes in joint_multiset_code_batches(6, single, double):
            left, right = prefixes[codes % split], suffixes[codes // split]
            vectors = np.logaddexp2(left[:, :, :, 0]+right[:, :, 0, None],
                                   left[:, :, :, 1]+right[:, :, 1, None])
            contractions = np.logaddexp2.reduce(vectors+weights[None, :, None], axis=1)
            contractions[:, 1] -= log_potential
            maximum = max(maximum, float(np.max(contractions)))
            checked += len(codes)
            branches += codes.size
        costs[single, double] = maximum
    if checked != comb(42, 6):
        raise ArithmeticError('The six-group support screen is incomplete')
    return dict(category_log2_costs={str(key): value for key, value in costs.items()},
                checked_multiset_cases=checked, cancellation_branches=branches,
                cached_prefix_matrices=split, elapsed_seconds=perf_counter()-start)


def column_polynomial_coefficients(caps, groups):
    """Exact source-support multiplicities for each total column degree.

    This counts source support arrays for fixed labelled constituent pairs;
    the independent choice of one of 28 constituent pairs per group is not
    part of this polynomial. Group multisets are only inequality checks.
    """
    if groups < 1 or set(caps) != set(categories(groups)):
        raise ValueError('A complete family of category caps is required')
    coefficients = [0 for degree in range(2*groups+1)]
    for (single, double), cap in caps.items():
        count = comb(groups, single)*comb(groups-single, double)*2**single
        coefficients[single+2*double] += count*cap
    return coefficients


def joint_log_count(costs, groups, n=32, k=16, q=256):
    """FLOAT positive polynomial tail, retaining total support at least 2r*d."""
    if not 0 < k <= n or q < 2 or set(costs) != set(categories(groups)):
        raise ValueError('Invalid code dimensions or category costs')
    polynomial = np.full(2*groups+1, -np.inf)
    for (single, double), cost in costs.items():
        degree = single+2*double
        multiplicity = comb(groups, single)*comb(groups-single, double)*2**single
        polynomial[degree] = np.logaddexp2(polynomial[degree], cost+log2(multiplicity))
    coefficients = np.array([0.])
    for column in range(n):
        following = np.full(len(coefficients)+2*groups, -np.inf)
        for degree, coefficient in enumerate(polynomial):
            following[degree:degree+len(coefficients)] = np.logaddexp2(
                following[degree:degree+len(coefficients)], coefficients+coefficient)
        coefficients = following
    floor = 2*groups*(n-k+1)
    coefficients += np.arange(len(coefficients))*log2(q-1)
    return float(np.logaddexp2.reduce(coefficients[floor:])-2*groups*(n-k)*log2(q-1))


def joint_count_arb(caps, groups, n=32, k=16, q=256):
    """Directed tail sum; caller sets precision and supplies positive Arb caps."""
    from flint import arb, arb_poly
    if not 0 < k <= n or q < 2:
        raise ValueError('Invalid code dimensions')
    coefficients = column_polynomial_coefficients(caps, groups)
    polynomial = arb_poly([coefficient*(q-1)**degree
                           for degree, coefficient in enumerate(coefficients)])
    expanded = polynomial**n
    floor = 2*groups*(n-k+1)
    return sum((expanded[degree] for degree in range(floor, 2*groups*n+1)), arb(0))/(q-1)**(2*groups*(n-k))


def resource_envelope_float(regions, log_potential, lower, upper, length=8):
    """Coordinatewise upper envelope for every labelled box-and-sum path."""
    if not 0 <= lower <= upper < len(regions) or length < 0:
        raise ValueError('Invalid occupancy interval or word length')
    width = upper-lower
    values = np.array([[0., log_potential]])
    for region in range(length):
        following = np.full((len(values)+width, 2), -np.inf)
        for delta in range(width+1):
            matrix = regions[lower+delta]
            action = np.stack([np.logaddexp2(matrix[row, 0]+values[:, 0],
                                             matrix[row, 1]+values[:, 1]) for row in range(2)], axis=1)
            following[delta:delta+len(values)] = np.maximum(following[delta:delta+len(values)], action)
        values = following
    return values


def arb_upper_max(left, right):
    """An exact endpoint bounding both inputs, including overlapping balls."""
    left, right = left.upper(), right.upper()
    return left if left >= right else right


def resource_envelope_arb(regions, potential, lower, upper, length=8):
    """Positive-matrix recurrence with explicit directed upper-endpoint maxima."""
    from flint import arb
    if not 0 <= lower <= upper < len(regions) or length < 0:
        raise ValueError('Invalid occupancy interval or word length')
    zero, one = arb(0), arb(1)
    width = upper-lower
    values = [(one, potential)]
    for region in range(length):
        following = [(zero, zero) for degree in range(len(values)+width)]
        for delta in range(width+1):
            a, b, c, d = regions[lower+delta]
            for degree, (first, second) in enumerate(values):
                left, right = following[degree+delta]
                following[degree+delta] = (arb_upper_max(left, a*first+b*second),
                                           arb_upper_max(right, c*first+d*second))
        values = following
    return values


def occupancy_envelope_caps_float(groups, z, log_potential):
    """Average the exact cancellation-count law after a box-and-sum envelope."""
    if not 1 <= groups <= MAX_ENVELOPE_GROUPS or not 0 < z < 1 or log_potential > 0:
        raise ValueError('One to thirty-two groups and valid marks required')
    regions = region_operators(z, groups)
    @lru_cache(None)
    def envelope(lower, upper):
        return resource_envelope_float(regions, log_potential, lower, upper)
    costs = {}
    for single, double in categories(groups):
        values = np.full(2, -np.inf)
        for cancelled in range(double+1):
            weight = (log2(comb(double, cancelled))+cancelled*log2(1/(Q-1))
                      +(double-cancelled)*log2((Q-2)/(Q-1)))
            action = envelope(double-cancelled, single+double)[7*single+2*cancelled]
            values = np.logaddexp2(values, action+weight)
        costs[single, double] = float(np.max(values-np.array([0., log_potential])))
    return costs


def occupancy_envelope_screen(groups=6, z=.99814, log_potential=-3997/64):
    start = perf_counter()
    costs = occupancy_envelope_caps_float(groups, z, log_potential)
    count = joint_log_count(costs, groups)
    bound = (log2(comb(REGION_LENGTH, groups))+groups*log2(28)-log2(Q-1)
             -LIMIT*log2(z)-log_potential+count)
    return dict(active_groups=groups, mark=z, log2_potential=log_potential,
                joint_log2_upper_screen=bound,
                category_log2_caps={str(key): value for key, value in costs.items()},
                source_multisets_enumerated=0, resource_tables=comb(groups+2, 2),
                resource_envelope=True, directed=False, whole_code_certificate=False,
                elapsed_seconds=perf_counter()-start)


def occupancy_envelope_receipt(groups=6, z='49907/50000', log_potential='-3997/64', precision=256,
                               common_pair=False):
    """Directed polynomial-time category envelope and aggregate-floor count."""
    from flint import arb, ctx
    z, log_potential = Fraction(z), Fraction(log_potential)
    if not 1 <= groups <= MAX_ENVELOPE_GROUPS or not 0 < z < 1 or log_potential > 0 or precision < 128:
        raise ValueError('One to thirty-two groups, valid marks, and at least 128 bits required')
    start = perf_counter()
    with ctx.workprec(precision):
        zero, one = arb(0), arb(1)
        potential = arb(2)**(arb(log_potential.numerator)/log_potential.denominator)
        regions = arb_region_matrices(z, groups)
        @lru_cache(None)
        def envelope(lower, upper):
            return resource_envelope_arb(regions, potential, lower, upper)
        caps = {}
        for single, double in categories(groups):
            first, second = zero, zero
            for cancelled in range(double+1):
                weight = arb(comb(double, cancelled)*(Q-2)**(double-cancelled))/(Q-1)**double
                left, right = envelope(double-cancelled, single+double)[7*single+2*cancelled]
                first = first+weight*left
                second = second+weight*right
            if (single, double) == (0, 0):
                if not first == one or not second < potential:
                    raise ArithmeticError('The empty column failed its directed contraction')
                caps[single, double] = one
            else:
                caps[single, double] = arb_upper_max(first, second/potential)
        count = joint_count_arb(caps, groups)
        zz = arb(z.numerator)/z.denominator
        pair_choices = 28 if common_pair else 28**groups
        bound = comb(REGION_LENGTH, groups)*pair_choices*count*zz**(-LIMIT)/(Q-1)/potential
        logarithm = bound.log()/arb(2).log()
        if not logarithm < -40:
            raise ArithmeticError(f'Directed envelope did not reach 40 bits: {logarithm}')
        return dict(active_groups=groups, active_rows_per_group=2, mark=str(z),
                    log2_potential=str(log_potential), failure_log2_upper=str(logarithm),
                    total_support_floor=34*groups, polynomial_degree=64*groups,
                    source_multisets_enumerated=0, covered_source_multisets=comb(36+groups, groups),
                    resource_tables=envelope.cache_info().currsize,
                    resource_state_vectors=sum((groups-width+1)*(8*width+1) for width in range(groups+1)),
                    cancellation_count_mixture=True, labelled_region_order=True,
                    coordinatewise_resource_envelope=True, individual_support_floors_relaxed=True,
                    exact_category_multiplicities=True, source_row_pair_choices=pair_choices,
                    constituent_pair_geometry='one shared pair' if common_pair else 'arbitrary pairs',
                    full_homogeneous_group_class=not common_pair,
                    projective_normalizations=1, precision_bits=precision,
                    directed=True, whole_code_certificate=False, elapsed_seconds=perf_counter()-start)


SAVED_ENVELOPE_WITNESSES = {
    6: ('49907/50000', '-3997/64'),
    7: ('997990/1000000', '-4002/64'),
    8: ('997840/1000000', '-4007/64'),
    9: ('997740/1000000', '-4004/64'),
    10: ('997590/1000000', '-4001/64'),
    11: ('997440/1000000', '-3998/64'),
    12: ('997290/1000000', '-3995/64'),
    13: ('997140/1000000', '-3992/64'),
    14: ('996990/1000000', '-3997/64'),
    15: ('996840/1000000', '-4002/64'),
    16: ('99669/100000', '-4063/64'),
}


def replay_envelope(groups=6, precision=256):
    return occupancy_envelope_receipt(groups, *SAVED_ENVELOPE_WITNESSES[groups], precision=precision)


SAVED_COMMON_PAIR_WITNESSES = {17: ('996548217/1000000000', '-127/2')}


def replay_common_pair(groups=17, precision=256):
    """Only the subclass in which every active group uses one shared row pair."""
    return occupancy_envelope_receipt(groups, *SAVED_COMMON_PAIR_WITNESSES[groups],
                                      precision=precision, common_pair=True)


def replay_envelope_family(precision=256):
    """Disjoint r=6,...,16 class union, using only literal fixed witnesses."""
    from flint import arb, ctx
    with ctx.workprec(precision):
        results = [replay_envelope(groups, precision) for groups in range(6, 17)]
        total = sum((arb(2)**arb(result['failure_log2_upper']) for result in results), arb(0))
        logarithm = total.log()/arb(2).log()
        return dict(active_rows_per_group=2, first_active_groups=6, last_active_groups=16,
                    failure_log2_upper=str(logarithm), directed=True,
                    precision_bits=precision, whole_code_certificate=False, classes=results)


def fixed_screen(groups=5, z=.99829, log_potential=-62.375):
    """One bounded FLOAT screen; no construction or parameter search."""
    result = (six_group_category_screen(z, log_potential) if groups == 6 else
              streaming_multiset_screen(groups, z, log_potential))
    costs = {category: result['category_log2_costs'][str(category)] for category in categories(groups)}
    count = joint_log_count(costs, groups)
    bound = (log2(comb(REGION_LENGTH, groups))+groups*log2(28)-log2(Q-1)
             -LIMIT*log2(z)-log_potential+count)
    return dict(active_groups=groups, active_rows_per_group=2, mark=z,
                log2_potential=log_potential, joint_log2_upper_screen=bound,
                separable_log2_upper_screen=result.get('log2_upper_screen'),
                joint_log2_message_count=count, total_support_floor=34*groups,
                polynomial_degree=64*groups,
                category_log2_caps={str(key): value for key, value in costs.items()},
                checked_support_cases=result['checked_multiset_cases'],
                cancellation_branches=result.get('cancellation_branches'),
                elapsed_seconds=result['elapsed_seconds'],
                exact_category_multiplicities=True, individual_support_floors_relaxed=True,
                support_count_bound='shortened-code information set',
                cancellation_before_potential=True, column_chronology=True,
                directed=False, whole_code_certificate=False)


def receipt(groups, z, log_potential, log_caps, precision=256):
    """Check all support mixtures, then the joint positive polynomial tail."""
    from flint import arb, ctx
    if groups not in range(1, 7) or precision < 128:
        raise ValueError('One to six groups and at least 128-bit precision required')
    z, log_potential = Fraction(z), Fraction(log_potential)
    log_caps = {category: Fraction(value) for category, value in log_caps.items()}
    if not 0 < z < 1 or log_potential > 0 or set(log_caps) != set(categories(groups)) or log_caps[0, 0] != 0:
        raise ValueError('Invalid marks or incomplete category caps')
    def rational(value):
        return arb(value.numerator)/value.denominator
    start = perf_counter()
    with ctx.workprec(precision):
        zero, one = arb(0), arb(1)
        potential = arb(2)**rational(log_potential)
        caps = {category: arb(2)**rational(value) for category, value in log_caps.items()}
        prefixes, suffixes = factored_column_vectors(arb_region_matrices(z, groups), potential)
        # Coin mixtures repeat many occupancy vectors. This bounded cache keeps
        # their directed actions, without constructing the full (r+1)^8 table.
        @lru_cache(16384)
        def action(code):
            return factored_word_vector(code, prefixes, suffixes)
        checked = 0
        for single, double in categories(groups):
            cap, terminal = caps[single, double], caps[single, double]*potential
            weights = [arb((Q-2)**(double-branch.bit_count()))/(Q-1)**double
                       for branch in range(1 << double)]
            index = 0
            for batch in joint_multiset_code_batches(groups, single, double):
                for codes in batch:
                    first, second = zero, zero
                    for code, weight in zip(codes, weights):
                        left, right = action(int(code))
                        first = first+weight*left
                        second = second+weight*right
                    if (single, double) == (0, 0):
                        valid = first == one and second < potential
                    else:
                        valid = first < cap and second < terminal
                    if not valid:
                        raise ArithmeticError('Directed joint category failed for '
                                              f'category {(single, double)}, case {index}')
                    checked += 1
                    index += 1
        if checked != comb(36+groups, groups):
            raise ArithmeticError('The directed support checks are incomplete')
        count = joint_count_arb(caps, groups)
        bound = comb(REGION_LENGTH, groups)*28**groups*count*rational(z)**(-LIMIT)/(Q-1)/potential
        logarithm = bound.log()/arb(2).log()
        if not logarithm < -40:
            raise ArithmeticError(f'Directed joint class bound did not reach 40 bits: {logarithm}')
        return dict(active_groups=groups, active_rows_per_group=2, mark=str(z),
                    log2_potential=str(log_potential),
                    category_log2_caps={str(key): str(value) for key, value in log_caps.items()},
                    failure_log2_upper=str(logarithm), checked_support_cases=checked,
                    source_row_pair_choices=28**groups, total_support_floor=34*groups,
                    polynomial_degree=64*groups,
                    cached_prefix_matrices=len(prefixes), occupancy_cache_limit=16384,
                    exact_category_multiplicities=True, individual_support_floors_relaxed=True,
                    support_count_bound='shortened-code information set',
                    cancellation_before_potential=True, column_chronology=True,
                    group_permutation_reduction=True, region_permutation_reduction=False,
                    projective_normalizations=1, precision_bits=precision,
                    directed=True, whole_code_certificate=False, elapsed_seconds=perf_counter()-start)


def propose_and_replay(groups=5, precision=256, z=Fraction(99829, 100000), potential=Fraction(-499, 8)):
    z, potential = Fraction(z), Fraction(potential)
    screen = fixed_screen(groups, float(z), float(potential))
    logs = {category: Fraction(ceil(screen['category_log2_caps'][str(category)]*4096)+1, 4096)
            for category in categories(groups)}
    logs[0, 0] = Fraction(0)
    return receipt(groups, z, potential, logs, precision)


SAVED_WITNESSES = {
    5: ('99829/100000', '-499/8', {
        (0, 0): '0', (0, 1): '-92257/1024', (0, 2): '-411519/4096',
        (0, 3): '-451621/4096', (0, 4): '-122509/1024', (0, 5): '-527129/4096',
        (1, 0): '-536603/4096', (1, 1): '-534427/4096', (1, 2): '-266133/2048',
        (1, 3): '-66265/512', (1, 4): '-131997/1024',
        (2, 0): '-534693/4096', (2, 1): '-532531/4096', (2, 2): '-530383/4096',
        (2, 3): '-264125/2048', (3, 0): '-133199/1024', (3, 1): '-265323/2048',
        (3, 2): '-528511/4096', (4, 0): '-530909/4096',
        (4, 1): '-132193/1024', (5, 0): '-264517/2048',
    }),
    6: ('49907/50000', '-3997/64', {
        (0, 0): '0', (0, 1): '-375744/4096', (0, 2): '-418753/4096',
        (0, 3): '-459372/4096', (0, 4): '-498297/4096', (0, 5): '-535921/4096',
        (0, 6): '-572506/4096', (1, 0): '-583959/4096', (1, 1): '-581865/4096',
        (1, 2): '-579787/4096', (1, 3): '-577722/4096', (1, 4): '-575672/4096',
        (1, 5): '-573635/4096', (2, 0): '-582122/4096', (2, 1): '-580041/4096',
        (2, 2): '-577975/4096', (2, 3): '-575923/4096', (2, 4): '-573885/4096',
        (3, 0): '-580296/4096', (3, 1): '-578228/4096', (3, 2): '-576174/4096',
        (3, 3): '-574135/4096', (4, 0): '-578481/4096', (4, 1): '-576426/4096',
        (4, 2): '-574384/4096', (5, 0): '-576677/4096', (5, 1): '-574634/4096',
        (6, 0): '-574884/4096',
    }),
}


def replay(groups=5, precision=256):
    z, potential, caps = SAVED_WITNESSES[groups]
    return receipt(groups, z, potential, caps, precision)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups', type=int, choices=range(1, MAX_ENVELOPE_GROUPS+1), default=5)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--directed', action='store_true')
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--envelope', action='store_true')
    parser.add_argument('--common-pair', action='store_true', help='Restrict to one pair shared by all active groups')
    parser.add_argument('--replay-envelope-family', action='store_true')
    parser.add_argument('--mark', default='99829/100000')
    parser.add_argument('--potential', default='-499/8')
    args = parser.parse_args()
    if args.groups > 6 and not (args.envelope or args.replay_envelope_family):
        parser.error('More than six groups requires the occupancy envelope')
    if args.common_pair and (not args.envelope or args.replay_envelope_family):
        parser.error('The common-pair restriction requires an individual envelope receipt')
    if args.common_pair and not (args.directed or args.replay):
        parser.error('The common-pair option is only implemented for directed receipts')
    result = (replay_envelope_family(args.precision) if args.replay_envelope_family else
              replay_common_pair(args.groups, args.precision) if args.envelope and args.common_pair and args.replay else
              replay_envelope(args.groups, args.precision) if args.envelope and args.replay else
              occupancy_envelope_receipt(args.groups, args.mark, args.potential, args.precision, args.common_pair) if args.envelope and args.directed else
              occupancy_envelope_screen(args.groups, float(Fraction(args.mark)), float(Fraction(args.potential))) if args.envelope else
              replay(args.groups, args.precision) if args.replay else
              propose_and_replay(args.groups, args.precision, args.mark, args.potential) if args.directed else
              fixed_screen(args.groups, float(Fraction(args.mark)), float(Fraction(args.potential))))
    print(json.dumps(result), flush=True)
