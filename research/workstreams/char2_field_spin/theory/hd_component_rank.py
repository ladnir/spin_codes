"""Component-class contractions and bounded gates; not a distance certificate.

Exact small kernels use the conditioned four-product post law. The positive
two-state control also covers hd_three_product.py, without equating ensembles.
"""

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial

from hd_additive_two_state import support_law
from hd_three_product import collapsed_majorant
from hd_two_state_regional_strides import ROWS, STRIDES


def identity(size):
    return [[int(i == j) for j in range(size)] for i in range(size)]


def multiply(left, right):
    size = len(left)
    return [[sum(left[i][k]*right[k][j] for k in range(size))
             for j in range(size)] for i in range(size)]


def power(matrix, exponent):
    result = identity(len(matrix))
    while exponent:
        if exponent & 1:
            result = multiply(result, matrix)
        matrix = multiply(matrix, matrix)
        exponent >>= 1
    return result


def exact_kernels(rows, q, mark, collapsed=False):
    """All normalized fixed-occupancy kernels, using exact polynomial DP."""
    mark = F(mark)
    if collapsed:
        transitions = collapsed_majorant(mark, q)
    else:
        transitions = []
        for active in (0, 1):
            transitions.append([
                [sum(p*mark**out for (new, out), p in support_law(old, active, q).items()
                     if new == target) for target in range(4)] for old in range(4)])
    size = len(transitions[0])
    polynomial = [identity(size)]
    for _ in range(rows):
        following = [[[F(0) for _ in range(size)] for _ in range(size)]
                     for _ in range(len(polynomial)+1)]
        for h, coefficient in enumerate(polynomial):
            for active in (0, 1):
                term = multiply(coefficient, transitions[active])
                for i in range(size):
                    for j in range(size):
                        following[h+active][i][j] += term[i][j]
        polynomial = following
    return [[[entry/comb(rows, h) for entry in row] for row in coefficient]
            for h, coefficient in enumerate(polynomial)]


def local_subset_polynomial(lengths, q):
    result = {0: 1}
    for length in lengths:
        following = dict(result)
        for h, count in result.items():
            following[h+length] = following.get(h+length, 0)+(q-1)*count
        result = following
    return result


def component_moment(cuts_by_region, rows, q, kernel):
    """Sum E[v^wt] over all nonzero field vectors in the contracted space.

    Each region has the same positive number of cuts. The component before
    its first cut is carried from the preceding region; the initial one is 0.
    """
    size = len(kernel(0))
    vector = [[0 for _ in range(size)] for _ in range(2)]
    vector[0][0] = 1
    for cuts in cuts_by_region:
        cuts = tuple(sorted(cuts))
        if not cuts or cuts[0] < 0 or cuts[-1] >= rows or len(set(cuts)) != len(cuts):
            raise ValueError("Need distinct in-range cuts in every region")
        interior = [right-left for left, right in zip(cuts, cuts[1:])]
        local = local_subset_polynomial(interior, q)
        following = [[0 for _ in range(size)] for _ in range(2)]
        for old_component, new_component in product((0, 1), repeat=2):
            shift = cuts[0]*old_component+(rows-cuts[-1])*new_component
            for inside, multiplicity in local.items():
                matrix = kernel(shift+inside)
                weight = multiplicity*(q-1)**new_component
                for old in range(size):
                    for new in range(size):
                        following[new_component][new] += weight*vector[old_component][old]*matrix[old][new]
        vector = following
    return sum(sum(row) for row in vector)-1


def direct_component_moment(cuts_by_region, rows, q, kernels):
    """Exhaustive component-flag oracle, independently of the frontier DP."""
    count = sum(map(len, cuts_by_region))
    size, total = len(kernels[0]), F(0)
    for flags in product((0, 1), repeat=count):
        if not any(flags):
            continue
        offset, carried = 0, 0
        vector = [F(1)]+[F(0)]*(size-1)
        for cuts in cuts_by_region:
            h, previous = carried*cuts[0], cuts[0]
            for index, cut in enumerate(cuts):
                end = cuts[index+1] if index+1 < len(cuts) else rows
                flag = flags[offset]
                offset += 1
                h += flag*(end-cut)
                carried = flag
            matrix = kernels[h]
            vector = [sum(vector[i]*matrix[i][j] for i in range(size)) for j in range(size)]
        total += (q-1)**sum(flags)*sum(vector)
    return total


def gap_moment_exact(rows, h, selected, base):
    """Moment of selected zero gaps in a uniform h-subset of [rows]."""
    missing, slots = rows-h, h+1
    if selected == 0:
        return F(1)
    if selected == slots:
        return base**(-missing)
    return sum(F(comb(g+selected-1, selected-1)
                 * comb(missing-g+slots-selected-1, slots-selected-1), comb(rows, h))
               * base**(-g) for g in range(missing+1))


def run_kernel_exact(rows, h, q, mark):
    """Closed positive run sum for the two-state control, exact small oracle."""
    mark, p = F(mark), F(1, q-1)
    base = mark+(1-mark)*p
    if h == 0:
        return [[F(1), F(0)], [F(0), base**rows]]
    result = [[F(0), F(0)], [F(0), F(0)]]
    for t in range((h-1)//2+1):
        result[0][1] += (comb(h-t-1, t)*mark**(2*t+1)*p**(2*t)
                         * base**(rows-2*t-1)*gap_moment_exact(rows, h, t+1, base))
    result[1][0] = p*p*result[0][1]
    for t in range(h//2+1):
        result[1][1] += (comb(h-t, t)*mark**(2*t)*p**(2*t)*base**(rows-2*t)
                         * gap_moment_exact(rows, h, t, base))
    for t in range(1, h//2+1):
        result[0][0] += (comb(h-t-1, t-1)*mark**(2*t)*p**(2*t)*base**(rows-2*t)
                         * gap_moment_exact(rows, h, t+1, base))
    return result


def arb_fraction(value):
    from flint import arb
    value = F(value)
    return arb(value.numerator)/value.denominator


def gap_moment_upper(rows, h, selected, base):
    """Directed upper via the positive beta-binomial factorial series."""
    from flint import arb
    missing, slots = rows-h, h+1
    if selected == 0 or missing == 0:
        return arb(1)
    if selected == slots:
        return base**(-missing)
    epsilon = 1/base-1
    term = total = arb(1)
    for ell in range(missing):
        # Bound every remaining term ratio by missing*epsilon/(ell+1).
        ratio = missing*epsilon/(ell+1)
        if ratio < 1:
            tail = term*ratio/(1-ratio)
            if tail < arb(2)**-110:
                return total+arb(tail.upper())
        term *= arb(missing-ell)/(ell+1)*arb(selected+ell)/(slots+ell)*epsilon
        total += term
    return total


def run_kernel_upper(rows, h, q, mark, returns=3):
    """Directed two-state kernel upper; return tail has a uniform majorant."""
    from flint import arb
    mark, p = arb_fraction(mark), arb(1)/(q-1)
    base = mark+(1-mark)*p
    if h == 0:
        return [[arb(1), arb(0)], [arb(0), base**rows]]
    result = [[arb(0), arb(0)], [arb(0), arb(0)]]
    for t in range(min((h-1)//2, returns)+1):
        result[0][1] += (comb(h-t-1, t)*mark**(2*t+1)*p**(2*t)*base**(rows-2*t-1)
                         * gap_moment_upper(rows, h, t+1, base))
    for t in range(min(h//2, returns)+1):
        result[1][1] += (comb(h-t, t)*mark**(2*t)*p**(2*t)*base**(rows-2*t)
                         * gap_moment_upper(rows, h, t, base))
    for t in range(1, min(h//2, returns+1)+1):
        result[0][0] += (comb(h-t-1, t-1)*mark**(2*t)*p**(2*t)*base**(rows-2*t)
                         * gap_moment_upper(rows, h, t+1, base))
    lam = h*p*p
    next_t = returns+1
    if not lam/(next_t+1) < 1:
        raise ValueError("Increase retained returns for this parameter range")
    tail = lam**next_t/factorial(next_t)/(1-lam/(next_t+1))
    result[0][1] += arb(tail.upper())
    result[1][1] += arb(tail.upper())
    result[0][0] += arb((p*p*tail).upper())
    result[1][0] = p*p*result[0][1]
    return result


def affine_cuts(active_rows, rows=ROWS, strides=STRIDES):
    return [tuple(sorted((pow(stride, -1, rows)*r) % rows for r in active_rows))
            for stride in strides]


def common_mark_obstruction():
    n, k, cutoff = 64*ROWS, 32*ROWS, (2*64*ROWS)//5
    return dict(small_mark_log2_lower_strict=cutoff-ROWS-512,
                large_mark_log2_lower=15*(k-1)-(n-cutoff),
                applies_to_common_mark_only=True, code_counterexample=False)


def singleton_receipts(denominator=2950, selected=(0, 1, 2, 8192, 16384, 32767), precision=128):
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        mark, q, cutoff = F(denominator-1, denominator), 65536, (2*64*ROWS)//5
        @lru_cache(None)
        def kernel(h):
            return run_kernel_upper(ROWS, h, q, mark)
        results = []
        for row in selected:
            moment = component_moment(affine_cuts((row,)), ROWS, q, kernel)
            # Dropping (1-z)^(R-1) increases this class contribution.
            upper = moment*arb(q-1)**-33*arb_fraction(mark)**(-cutoff)
            results.append(dict(row=row, log2_upper=str(upper.log()/arb(2).log()),
                                below_2_to_minus_40=bool(upper < arb(2)**-40)))
        return dict(mark=str(mark), precision=precision, kernels_evaluated=kernel.cache_info().currsize,
                    selected_profiles=results, all_singletons_bound=False,
                    whole_distance_certificate=False)
    finally:
        ctx.prec = previous


def singleton_uniform_receipt(denominator=2950, precision=128):
    """Uniform in all singleton cut positions; covers every regionwise placement."""
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        q, regions, rows = 65536, 64, ROWS
        mark_fraction = F(denominator-1, denominator)
        mark, p = arb_fraction(mark_fraction), arb(1)/(q-1)
        base = mark+(1-mark)*p
        base_power = base**rows
        eta = (1+p*p)**rows-1
        # No-return birth is maximized when h=1 (uniform first live position).
        birth = mark*(1-base_power)/(rows*(1-base))
        terminal_return = p*p*(1+p*p)**(rows-1)
        positive = [[terminal_return, birth+eta],
                    [p*p*(birth+eta), base_power+eta]]
        arbitrary = [[arb(1), birth+eta],
                     [p*p*(birth+eta), base_power+eta]]
        empty = [[arb(1), arb(0)], [arb(0), base_power]]
        full = run_kernel_upper(rows, rows, q, mark_fraction)
        matrix = [[arb(0) for _ in range(4)] for _ in range(4)]
        for old_component, new_component in product((0, 1), repeat=2):
            block = ((empty, positive), (arbitrary, full))[old_component][new_component]
            for old in range(2):
                for new in range(2):
                    matrix[2*old_component+old][2*new_component+new] = (q-1)**new_component*block[old][new]
        moment = sum(power(matrix, regions)[0])-1
        cutoff = (2*regions*rows)//5
        per_profile = moment*arb(q-1)**-33*mark**(-cutoff)
        class_upper = rows*per_profile
        if not class_upper < arb(2)**-40:
            raise ArithmeticError("The uniform singleton gate did not pass")
        return dict(mark=str(mark_fraction), precision=precision,
                    per_profile_log2_upper=str(per_profile.log()/arb(2).log()),
                    all_singletons_log2_upper=str(class_upper.log()/arb(2).log()),
                    all_singletons_below_2_to_minus_40=True,
                    any_fixed_regionwise_row_order=True,
                    exact_four_product_and_three_product_control=True,
                    whole_distance_certificate=False)
    finally:
        ctx.prec = previous


def pair_receipts(denominator=1475, precision=128):
    """A few a=2 feasibility profiles, not a bound for the full pair class."""
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        q, cutoff = 65536, (2*64*ROWS)//5
        mark = F(denominator-1, denominator)
        @lru_cache(None)
        def kernel(h):
            return run_kernel_upper(ROWS, h, q, mark)
        profiles = [(str(pair), affine_cuts(pair)) for pair in ((0, 1), (0, 16384), (1, 32767))]
        profiles.append(("common cuts (0,1), e.g. aligned group pair", [(0, 1)]*64))
        results = []
        for name, cuts in profiles:
            moment = component_moment(cuts, ROWS, q, kernel)
            upper = moment*arb(q-1)**-65*arb_fraction(mark)**(-cutoff)
            results.append(dict(profile=name, log2_upper=str(upper.log()/arb(2).log()),
                                below_2_to_minus_40=bool(upper < arb(2)**-40)))
        return dict(mark=str(mark), precision=precision, kernels_evaluated=kernel.cache_info().currsize,
                    selected_profiles=results, all_pairs_bound=False, whole_distance_certificate=False)
    finally:
        ctx.prec = previous


def main():
    print("common-mark obstruction:", common_mark_obstruction())
    print("singleton profile receipt:", singleton_receipts())
    print("all-singleton class receipt:", singleton_uniform_receipt())
    print("selected pair feasibility:", pair_receipts())


if __name__ == '__main__':
    main()
