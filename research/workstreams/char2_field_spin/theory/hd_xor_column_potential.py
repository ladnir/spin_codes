"""Eight-region potentials for up to four XOR2 groups with two active rows.

Every column's cancellation mixture is averaged before the potential bound.
All source-row-pair assignments are covered. Exact RS overlap counts retain
both support floors. Directed checks preserve labelled region chronology;
no whole-code certificate is claimed.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, islice
from math import ceil, comb, log2
from time import perf_counter

import numpy as np
from scipy.optimize import minimize_scalar

from hd_xor_column_probe import mds_prescribed_support
from hd_xor_multigroup import (LIMIT, Q, REGION_LENGTH, region_operators,
                              nonzero_arb_coefficients)

PAIRS = tuple((1 << a) | (1 << b) for a, b in combinations(range(8), 2))
FULL = 255


def mask_code(mask, base=3):
    return sum(base**index for index in range(8) if mask >> index & 1)


def output_support_law(source):
    """Exact mask law for a source support of size zero, one, or two."""
    if not 0 <= source <= FULL or source.bit_count() > 2:
        raise ValueError('A source support must contain at most two rows')
    if not source:
        return ((0, Fraction(1)),)
    if source.bit_count() == 1:
        return ((FULL ^ source, Fraction(1)),)
    return ((source, Fraction(1, Q-1)), (FULL, Fraction(Q-2, Q-1)))


@lru_cache(None)
def case_plans(groups=2):
    """All distinct input-support cases, grouped by (singletons,doubles).

    Singleton rows range over all eight positions. Double supports range
    over all 28 pairs; the double/double category keeps all 28^2 ordered
    assignments. No source-row-choice union factor belongs in these plans.
    """
    if groups not in (1, 2):
        raise ValueError('Only one or two groups are implemented')
    singles = tuple(1 << index for index in range(8))
    sources = {(0, 0): ((0, 0),), (1, 0): tuple((s, 0) for s in singles),
               (0, 1): tuple((s, 0) for s in PAIRS)}
    if groups == 2:
        sources.update({(2, 0): tuple((s, t) for s in singles for t in singles),
                        (1, 1): tuple((s, t) for s in singles for t in PAIRS),
                        (0, 2): tuple((s, t) for s in PAIRS for t in PAIRS)})
    result = {}
    for category, assignments in sources.items():
        cases = []
        for first, second in assignments:
            law = {}
            for left, p in output_support_law(first):
                for right, q in output_support_law(second):
                    code = mask_code(left)+mask_code(right)
                    law[code] = law.get(code, Fraction(0))+p*q
            cases.append((first, second, tuple(sorted(law.items()))))
        result[category] = tuple(cases)
    return result


def chronological_words(regions, length=8):
    """All fixed-length words, in base len(regions), retaining region order."""
    if len(regions.shape) != 3 or regions.shape[1:] != (2, 2):
        raise ValueError('A consecutive family of 2x2 regional matrices is required')
    words = np.array([[[0., -np.inf], [-np.inf, 0.]]])
    for region in range(length):
        following = []
        for occupancy in range(len(regions)):
            matrix = regions[occupancy]
            following.append(np.logaddexp2(words[:, :, 0, None]+matrix[None, 0, None, :],
                                            words[:, :, 1, None]+matrix[None, 1, None, :]))
        words = np.concatenate(following, axis=0)
    return words


def column_matrices(z, groups=2):
    """Average each exact XOR support mixture before taking any potential."""
    words = chronological_words(region_operators(z, 2))
    result = {}
    for category, plans in case_plans(groups).items():
        matrices = []
        for first, second, law in plans:
            terms = [words[code]+log2(probability) for code, probability in law]
            matrices.append(np.logaddexp2.reduce(terms))
        result[category] = np.array(matrices)
    return result


def category_costs(matrices, log_potential):
    """Worst log contraction, retaining each case's complete mask mixture."""
    potential = np.array([0., log_potential])
    costs, worst = {}, {}
    for category, values in matrices.items():
        contractions = np.logaddexp2(values[:, :, 0], values[:, :, 1]+log_potential)-potential
        index = np.unravel_index(np.argmax(contractions), contractions.shape)
        costs[category] = float(contractions[index])
        worst[category] = (int(index[0]), int(index[1]))
    return costs, worst


@lru_cache(None)
def overlap_counts(q=Q, n=32, k=16):
    """Exact integers indexed by (singleton columns, double columns)."""
    counts = {}
    for first in range(n-k+1, n+1):
        left = mds_prescribed_support(q, n, k, first)
        for second in range(n-k+1, n+1):
            right = mds_prescribed_support(q, n, k, second)
            for overlap in range(max(0, first+second-n), min(first, second)+1):
                multiplicity = comb(n, first)*comb(first, overlap)*comb(n-first, second-overlap)
                key = first+second-2*overlap, overlap
                counts[key] = counts.get(key, 0)+multiplicity*left*right
    return tuple(sorted(counts.items()))


@lru_cache(None)
def overlap_arrays():
    counts = overlap_counts()
    return np.array([degree for degree, _ in counts]), np.array([log2(count) for _, count in counts])


def group_log_moment(singleton_mark, double_mark):
    degrees, counts = overlap_arrays()
    return float(np.logaddexp2.reduce(counts+degrees[:, 0]*singleton_mark+degrees[:, 1]*double_mark))


def best_support_marks(costs, groups):
    first, second = costs[1, 0], costs[0, 1]
    if costs[0, 0] > 1e-12:
        raise ArithmeticError('The empty column must be a contraction')
    if groups == 1:
        return first, second
    first, second = max(first, costs[2, 0]/2), max(second, costs[0, 2]/2)
    if first+second >= costs[1, 1]:
        return first, second
    total = costs[1, 1]
    chosen = minimize_scalar(lambda single: group_log_moment(single, total-single),
                             bounds=(first, total-second), method='bounded', options={'xatol': 1e-9})
    return float(chosen.x), total-float(chosen.x)


def potential_bound(groups, z, log_potential, matrices=None):
    if groups not in (1, 2) or not 0 < z < 1:
        raise ValueError('One or two groups and a mark in (0,1) are required')
    if matrices is None:
        matrices = column_matrices(z, groups)
    costs, worst = category_costs(matrices, log_potential)
    first, second = best_support_marks(costs, groups)
    bound = (log2(comb(REGION_LENGTH, groups))+groups*log2(28)-log2(Q-1)
             -LIMIT*log2(z)+max(0., -log_potential)+groups*group_log_moment(first, second))
    return bound, first, second, costs, worst


def proposal(groups=2, iterations=40):
    start = perf_counter()
    best = None
    def objective(ell):
        nonlocal best
        z = 2**(-ell)
        matrices = column_matrices(z, groups)
        chosen = minimize_scalar(lambda p: potential_bound(groups, z, p, matrices)[0],
                                 bounds=(-150., 0.), method='bounded', options={'xatol': 1e-7})
        value, first, second, costs, worst = potential_bound(groups, z, float(chosen.x), matrices)
        if best is None or value < best['log2_upper_screen']:
            best = dict(active_groups=groups, active_rows_per_group=2, mark=z,
                        log2_potential=float(chosen.x), log2_singleton_mark=first,
                        log2_double_mark=second, log2_upper_screen=value,
                        category_log2_costs={str(key): number for key, number in costs.items()},
                        worst_cases={str(key): list(case_plans(groups)[key][index][:2])+[state]
                                     for key, (index, state) in worst.items()},
                        all_source_row_pairs=True, exact_rs_overlap=True,
                        cancellation_before_potential=True, column_chronology=True,
                        alpha=1., directed=False, whole_code_certificate=False)
        return value
    # A bounded scalar search; no construction changes are considered.
    minimize_scalar(objective, bounds=(1e-6, .03), method='bounded',
                    options={'maxiter': iterations, 'xatol': 1e-9})
    best['elapsed_seconds'] = perf_counter()-start
    return best


def arb_region_matrices(z, maximum):
    """Directed R_0,...,R_maximum with exact regional normalization."""
    from flint import arb, arb_poly
    one = arb(1)
    base = [arb_poly(entry) for entry in nonzero_arb_coefficients(z, maximum)]
    result = [arb_poly([one]), arb_poly([]), arb_poly([]), arb_poly([one])]
    def product(left, right):
        return [(left[2*i]*right[j]+left[2*i+1]*right[2+j]).truncate(maximum+1)
                for i in range(2) for j in range(2)]
    exponent = 128
    while exponent:
        if exponent & 1:
            result = product(result, base)
        exponent >>= 1
        if exponent:
            base = product(base, base)
    return [tuple(entry[m]/comb(REGION_LENGTH, m) for entry in result) for m in range(maximum+1)]


def arb_chronological_words(regions, length):
    """Little-endian count codes, with the first region multiplied on the left."""
    from flint import arb
    zero, one = arb(0), arb(1)
    words = [(one, zero, zero, one)]
    for region in range(length):
        words = [tuple(left[2*i]*right[j]+left[2*i+1]*right[2+j]
                       for i in range(2) for j in range(2))
                 for right in regions for left in words]
    return words


def arb_regions_and_words(z):
    """Directed R_0,R_1,R_2 and every eight-region chronological word."""
    regions = arb_region_matrices(z, 2)
    words = arb_chronological_words(regions, 8)
    return regions, words


def factored_column_vectors(regions, potential):
    """Cache four-region prefixes and suffixes applied to (1,potential)."""
    prefixes = arb_chronological_words(regions, 4)
    suffixes = [(matrix[0]+matrix[1]*potential, matrix[2]+matrix[3]*potential)
                for matrix in prefixes]
    return prefixes, suffixes


def factored_word_vector(code, prefixes, suffixes):
    """Evaluate the first-four matrix times the last-four potential vector."""
    right, left = divmod(int(code), len(prefixes))
    a, b, c, d = prefixes[left]
    x, y = suffixes[right]
    return a*x+b*y, c*x+d*y


def receipt(groups, z, log_potential, singleton_mark, double_mark, precision=256):
    """Replay all cases and the exact RS-overlap union with directed arithmetic."""
    from flint import arb, ctx
    if groups not in (1, 2) or precision < 128:
        raise ValueError('One or two groups and at least 128-bit precision required')
    z, log_potential, singleton_mark, double_mark = map(Fraction,
        (z, log_potential, singleton_mark, double_mark))
    if not 0 < z < 1 or log_potential > 0:
        raise ValueError('Invalid weight mark or state potential')
    def rational(value):
        return arb(value.numerator)/value.denominator
    start = perf_counter()
    with ctx.workprec(precision):
        zero, one = arb(0), arb(1)
        zz = rational(z)
        potential = arb(2)**rational(log_potential)
        f1, f2 = arb(2)**rational(singleton_mark), arb(2)**rational(double_mark)
        regions, words = arb_regions_and_words(z)
        checked = 0
        for category, plans in case_plans(groups).items():
            factor = f1**category[0]*f2**category[1]
            for first, second, law in plans:
                matrix = [sum((rational(probability)*words[code][entry] for code, probability in law), zero)
                          for entry in range(4)]
                contractions = (matrix[0]+matrix[1]*potential,
                                matrix[2]+matrix[3]*potential)
                if category == (0, 0):
                    valid = contractions[0] == one and contractions[1] < potential
                else:
                    valid = contractions[0] < factor and contractions[1] < factor*potential
                if not valid:
                    raise ArithmeticError(f'Directed potential failed for category {category}, sources {(first, second)}')
                checked += 1
        moment = sum((count*f1**single*f2**double for (single, double), count in overlap_counts()), zero)
        bound = (comb(REGION_LENGTH, groups)*28**groups*moment**groups
                 *zz**(-LIMIT)/(Q-1)/potential)
        logarithm = bound.log()/arb(2).log()
        if not logarithm < -40:
            raise ArithmeticError(f'Directed class bound did not reach 40 bits: {logarithm}')
        return dict(active_groups=groups, active_rows_per_group=2, mark=str(z),
                    log2_potential=str(log_potential), log2_singleton_mark=str(singleton_mark),
                    log2_double_mark=str(double_mark), failure_log2_upper=str(logarithm),
                    checked_support_cases=checked, source_row_pair_choices=28**groups,
                    exact_rs_overlap=True, cancellation_before_potential=True,
                    column_chronology=True, projective_normalizations=1,
                    alpha=1, precision_bits=precision, directed=True,
                    whole_code_certificate=False, elapsed_seconds=perf_counter()-start)


def propose_and_replay(groups=2, iterations=40, precision=256):
    proposed = proposal(groups, iterations)
    z = Fraction(round(proposed['mark']*10**9), 10**9)
    potential = Fraction(round(proposed['log2_potential']*64), 64)
    _, first, second, _, _ = potential_bound(groups, float(z), float(potential))
    first = Fraction(ceil(first*4096)+1, 4096)
    second = Fraction(ceil(second*4096)+1, 4096)
    return receipt(groups, z, potential, first, second, precision)


SAVED_WITNESSES = {
    1: ('499451537/500000000', '-3991/64', '-343125/4096', '-342223/4096'),
    2: ('199748791/200000000', '-499/8', '-12225/256', '-194985/4096'),
}


def replay(groups=2, precision=256):
    return receipt(groups, *SAVED_WITNESSES[groups], precision=precision)


@lru_cache(2)
def multiset_plans(groups):
    """Bounded group-symmetry reduction; region positions remain labelled.

    A singleton group chooses one missing output region. A double group
    chooses an edge on eight region vertices. Its cancellation either keeps
    that edge or occupies every region. Group permutations preserve this law.
    """
    if not 1 <= groups <= 4:
        raise ValueError('The multiset checker is capped at four groups')
    base = groups+1
    full = mask_code(FULL, base)
    pair_codes = tuple(mask_code(pair, base) for pair in PAIRS)
    result = {}
    for single in range(groups+1):
        single_rows = tuple(combinations_with_replacement(range(8), single))
        for double in range(groups-single+1):
            pair_rows = tuple(combinations_with_replacement(range(28), double))
            codes = np.empty((len(single_rows)*len(pair_rows), 1 << double), dtype=np.int32)
            weights = np.array([branch.bit_count()*log2(1/(Q-1))
                                +(double-branch.bit_count())*log2((Q-2)/(Q-1))
                                for branch in range(1 << double)])
            offset = 0
            for rows in single_rows:
                initial = (single+double)*full-sum(base**row for row in rows)
                for edges in pair_rows:
                    shifts = [pair_codes[edge]-full for edge in edges]
                    codes[offset] = [initial+sum(shifts[i] for i in range(double) if branch >> i & 1)
                                     for branch in range(1 << double)]
                    offset += 1
            result[single, double] = codes, weights
    if sum(len(codes) for codes, weights in result.values()) != comb(36+groups, groups):
        raise ArithmeticError('The multiset cases did not cover their combinatorial count')
    return result


def multiset_column_matrices(z, groups):
    words = chronological_words(region_operators(z, groups))
    result = {}
    for category, (indices, weights) in multiset_plans(groups).items():
        matrices = np.empty((len(indices), 2, 2))
        # Batches bound temporary storage even at four groups and 16 branches.
        for start in range(0, len(indices), 1024):
            selected = indices[start:start+1024]
            matrices[start:start+len(selected)] = np.logaddexp2.reduce(
                words[selected]+weights[None, :, None, None], axis=1)
        result[category] = matrices
    return result


def multiset_support_marks(costs):
    first = max(cost/single for (single, double), cost in costs.items() if single and not double)
    second = max(cost/double for (single, double), cost in costs.items() if double and not single)
    mixed = [(single, double, cost) for (single, double), cost in costs.items() if single and double]
    if not mixed:
        return first, second
    last = max([first]+[(cost-double*second)/single for single, double, cost in mixed])
    def second_for(single_mark):
        return max([second]+[(cost-single*single_mark)/double for single, double, cost in mixed])
    if last <= first:
        return first, second
    selected = minimize_scalar(lambda single_mark: group_log_moment(single_mark, second_for(single_mark)),
                               bounds=(first, last), method='bounded', options={'xatol': 1e-8})
    return float(selected.x), second_for(float(selected.x))


def multiset_proposal(groups, iterations=30):
    start = perf_counter()
    best = None
    def objective(ell):
        nonlocal best
        z = 2**(-ell)
        matrices = multiset_column_matrices(z, groups)
        def evaluate(potential):
            costs, worst = category_costs(matrices, potential)
            first, second = multiset_support_marks(costs)
            value = (log2(comb(REGION_LENGTH, groups))+groups*log2(28)-log2(Q-1)
                     -LIMIT*log2(z)-potential+groups*group_log_moment(first, second))
            return value, first, second, costs
        chosen = minimize_scalar(lambda p: evaluate(p)[0], bounds=(-150., 0.), method='bounded',
                                 options={'xatol': 1e-7})
        value, first, second, costs = evaluate(float(chosen.x))
        if best is None or value < best['log2_upper_screen']:
            best = dict(active_groups=groups, active_rows_per_group=2, mark=z,
                        log2_potential=float(chosen.x), log2_singleton_mark=first,
                        log2_double_mark=second, log2_upper_screen=value,
                        checked_multiset_cases=comb(36+groups, groups),
                        source_row_pair_choices=28**groups,
                        category_log2_costs={str(key): cost for key, cost in costs.items()},
                        exact_rs_overlap=True, cancellation_before_potential=True,
                        column_chronology=True, group_permutation_reduction=True,
                        directed=False, whole_code_certificate=False)
        return value
    minimize_scalar(objective, bounds=(1e-6, .03), method='bounded',
                    options={'maxiter': iterations, 'xatol': 1e-9})
    best['elapsed_seconds'] = perf_counter()-start
    return best


def multiset_receipt(groups, z, log_potential, singleton_mark, double_mark, precision=256):
    """Directed group-multiset checks, with no symmetry reduction of regions.

    Repeated edges keep independent cancellation coins: every labelled coin
    branch is averaged, including branches with duplicate output count codes.
    The cases certify universal inequalities, not additional union terms.
    """
    from flint import arb, ctx
    if groups not in (1, 2, 3, 4) or precision < 128:
        raise ValueError('One to four groups and at least 128-bit precision required')
    z, log_potential, singleton_mark, double_mark = map(Fraction,
        (z, log_potential, singleton_mark, double_mark))
    if not 0 < z < 1 or log_potential > 0:
        raise ValueError('Invalid weight mark or state potential')
    def rational(value):
        return arb(value.numerator)/value.denominator
    start = perf_counter()
    with ctx.workprec(precision):
        zero, one = arb(0), arb(1)
        zz = rational(z)
        potential = arb(2)**rational(log_potential)
        f1, f2 = arb(2)**rational(singleton_mark), arb(2)**rational(double_mark)
        regions = arb_region_matrices(z, groups)
        prefixes, suffixes = factored_column_vectors(regions, potential)
        checked = 0
        for (single, double), (indices, _) in multiset_plans(groups).items():
            factor = f1**single*f2**double
            terminal = factor*potential
            weights = [arb((Q-2)**(double-branch.bit_count()))/(Q-1)**double
                       for branch in range(1 << double)]
            for index, codes in enumerate(indices):
                first, second = zero, zero
                for code, weight in zip(codes, weights):
                    left, right = factored_word_vector(code, prefixes, suffixes)
                    first = first+weight*left
                    second = second+weight*right
                if (single, double) == (0, 0):
                    valid = first == one and second < potential
                else:
                    valid = first < factor and second < terminal
                if not valid:
                    raise ArithmeticError('Directed multiset potential failed for '
                                          f'category {(single, double)}, case {index}')
                checked += 1
        moment = sum((count*f1**single*f2**double for (single, double), count in overlap_counts()), zero)
        bound = (comb(REGION_LENGTH, groups)*28**groups*moment**groups
                 *zz**(-LIMIT)/(Q-1)/potential)
        logarithm = bound.log()/arb(2).log()
        if not logarithm < -40:
            raise ArithmeticError(f'Directed class bound did not reach 40 bits: {logarithm}')
        return dict(active_groups=groups, active_rows_per_group=2, mark=str(z),
                    log2_potential=str(log_potential), log2_singleton_mark=str(singleton_mark),
                    log2_double_mark=str(double_mark), failure_log2_upper=str(logarithm),
                    checked_support_cases=checked, source_row_pair_choices=28**groups,
                    cached_prefix_matrices=len(prefixes), cached_suffix_vectors=len(suffixes),
                    exact_rs_overlap=True, cancellation_before_potential=True,
                    column_chronology=True, group_permutation_reduction=True,
                    region_permutation_reduction=False, projective_normalizations=1,
                    alpha=1, precision_bits=precision, directed=True,
                    whole_code_certificate=False, elapsed_seconds=perf_counter()-start)


def multiset_propose_and_replay(groups, iterations=30, precision=256):
    proposed = multiset_proposal(groups, iterations)
    z = Fraction(round(proposed['mark']*10**9), 10**9)
    potential = Fraction(round(proposed['log2_potential']*64), 64)
    costs, _ = category_costs(multiset_column_matrices(float(z), groups), float(potential))
    first, second = multiset_support_marks(costs)
    first = Fraction(ceil(first*4096)+1, 4096)
    second = Fraction(ceil(second*4096)+1, 4096)
    return multiset_receipt(groups, z, potential, first, second, precision)


SAVED_MULTISET_WITNESSES = {
    **SAVED_WITNESSES,
    3: ('62411847/62500000', '-1997/32', '-145985/4096', '-72745/2048'),
    4: ('998438053/1000000000', '-3995/64', '-120989/4096', '-120555/4096'),
}


def multiset_replay(groups, precision=256):
    return multiset_receipt(groups, *SAVED_MULTISET_WITNESSES[groups], precision=precision)


def multiset_code_batches(groups, single, double, batch_size=1024):
    """Stream exact multiset count codes without a complete case or word table."""
    if not 1 <= groups <= 5 or min(single, double) < 0 or single+double > groups or batch_size < 1:
        raise ValueError('Invalid streaming multiset category or batch size')
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


def streaming_multiset_screen(groups, z, log_potential):
    """Fixed-witness FLOAT boundary screen, with O(batch_size*2^groups) storage."""
    if not 1 <= groups <= 5 or not 0 < z < 1 or log_potential > 0:
        raise ValueError('One to five groups and valid marks required')
    start = perf_counter()
    prefixes = chronological_words(region_operators(z, groups), 4)
    suffixes = np.logaddexp2(prefixes[:, :, 0], prefixes[:, :, 1]+log_potential)
    split = len(prefixes)
    costs, worst = {}, {}
    checked = 0
    for single in range(groups+1):
        for double in range(groups-single+1):
            weights = np.array([branch.bit_count()*log2(1/(Q-1))
                                +(double-branch.bit_count())*log2((Q-2)/(Q-1))
                                for branch in range(1 << double)])
            maximum, maximizer, offset = -np.inf, None, 0
            for codes in multiset_code_batches(groups, single, double):
                left, right = prefixes[codes % split], suffixes[codes // split]
                vectors = np.logaddexp2(left[:, :, :, 0]+right[:, :, 0, None],
                                       left[:, :, :, 1]+right[:, :, 1, None])
                contractions = np.logaddexp2.reduce(vectors+weights[None, :, None], axis=1)
                contractions[:, 1] -= log_potential
                index = np.unravel_index(np.argmax(contractions), contractions.shape)
                if contractions[index] > maximum:
                    maximum = float(contractions[index])
                    maximizer = (offset+int(index[0]), int(index[1]))
                offset += len(codes)
            costs[single, double] = maximum
            worst[single, double] = maximizer
            checked += offset
    if checked != comb(36+groups, groups):
        raise ArithmeticError('Streaming multiset count is incomplete')
    first, second = multiset_support_marks(costs)
    value = (log2(comb(REGION_LENGTH, groups))+groups*log2(28)-log2(Q-1)
             -LIMIT*log2(z)-log_potential+groups*group_log_moment(first, second))
    return dict(active_groups=groups, active_rows_per_group=2, mark=z,
                log2_potential=log_potential, log2_singleton_mark=first,
                log2_double_mark=second, log2_upper_screen=value,
                checked_multiset_cases=checked, cached_prefix_matrices=split,
                source_row_pair_choices=28**groups,
                category_log2_costs={str(key): cost for key, cost in costs.items()},
                worst_cases={str(key): case for key, case in worst.items()},
                exact_rs_overlap=True, cancellation_before_potential=True,
                column_chronology=True, group_permutation_reduction=True,
                directed=False, whole_code_certificate=False,
                elapsed_seconds=perf_counter()-start)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups', type=int, nargs='+', choices=(1, 2, 3, 4, 5), default=(1, 2))
    parser.add_argument('--iterations', type=int, default=40)
    parser.add_argument('--directed', action='store_true')
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--multiset', action='store_true')
    parser.add_argument('--stream-mark', type=float, help='Fixed-z FLOAT streaming screen (up to five groups)')
    parser.add_argument('--stream-potential', type=float, default=-62.4375)
    args = parser.parse_args()
    for groups in args.groups:
        if args.stream_mark is not None:
            result = streaming_multiset_screen(groups, args.stream_mark, args.stream_potential)
        elif args.multiset:
            result = (multiset_replay(groups, args.precision) if args.replay else
                      multiset_propose_and_replay(groups, args.iterations, args.precision) if args.directed
                      else multiset_proposal(groups, args.iterations))
        else:
            result = (replay(groups, args.precision) if args.replay else
                      propose_and_replay(groups, args.iterations, args.precision) if args.directed
                      else proposal(groups, args.iterations))
        print(json.dumps(result), flush=True)
