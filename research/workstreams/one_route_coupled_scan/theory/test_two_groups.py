"""Exact bounded tests of the two-group regional route and insertion sums.

Enumerates images of one full physical pair or two labeled singleton pairs
for m=2,3,4 batches. All inactive pair labels have the same number of
permutation completions, so uniform distinct images give their exact route
marginal. Independent swap bits are enumerated explicitly.

Toy matrices are distinct and noncommuting; they are not state probabilities.
This checks route averaging and directed arithmetic, not a distance claim.
Writes no files and runs no work at import; invoke with Python -B.
"""

from collections import Counter
from decimal import Decimal as D, ROUND_CEILING, localcontext
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb

from check_endpoint_route import gap_convolve, regions
from check_two_groups import insertion_moment, regional_matrices


MASKS = (1, 2, 3, 4, 5, 6, 8, 9, 10, 12)
SINGLETONS = (1, 2, 4, 8)
FULL_PAIRS = (3, 12)
CROSS_PAIRS = (5, 6, 9, 10)
PRECISION = 80
RELATIVE_TOLERANCE = F(1, 10**70)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def zero():
    return [[F(0), F(0)], [F(0), F(0)]]


def identity():
    return [[F(1), F(0)], [F(0), F(1)]]


def exact_multiply(left, right):
    # Separate exact oracle; do not reuse the Decimal production multiplier.
    return [[sum((left[i][k] * right[k][j] for k in range(2)), F(0))
             for j in range(2)] for i in range(2)]


def add_in_place(total, value):
    for i in range(2):
        for j in range(2):
            total[i][j] += value[i][j]


def divide(matrix, denominator):
    return [[entry / denominator for entry in row] for row in matrix]


def average(matrices, masks):
    total = zero()
    for mask in masks:
        add_in_place(total, matrices[mask])
    return divide(total, len(masks))


def toy_matrices():
    matrices = {
        mask: [[F(0), F(mask + 2, 17)],
               [F(2 * mask + 1, 23), F(mask * mask + 3, 29)]]
        for mask in MASKS
    }
    for left, right in combinations(MASKS, 2):
        require(matrices[left] != matrices[right], "toy masks must differ")
        require(exact_multiply(matrices[left], matrices[right])
                != exact_multiply(matrices[right], matrices[left]),
                ("toy matrices must not commute", left, right))
    return matrices


def chronological_product(empty, matrices, masks):
    result = identity()
    for mask in masks:
        result = exact_multiply(result, matrices[mask] if mask else empty)
    return result


def exact_routes(batches, empty, matrices):
    pair_total, split_total = zero(), zero()
    pair_samples = split_samples = collisions = 0
    full_masks, cross_masks = Counter(), Counter()
    slots = 2 * batches

    # A complete pair occupies a uniform physical pair slot. Its swap bit
    # changes labels but not its support; count both values nonetheless.
    for slot, swap in product(range(slots), (0, 1)):
        masks = [0] * batches
        offset = 2 * (slot % 2)
        masks[slot // 2] = (1 << (offset + swap)) | (1 << (offset + 1 - swap))
        full_masks[masks[slot // 2]] += 1
        add_in_place(pair_total, chronological_product(empty, matrices, masks))
        pair_samples += 1

    # The two source pairs are labeled. Their images are ordered and distinct,
    # not independent samples with replacement. Swaps act on each singleton.
    for first, second in permutations(range(slots), 2):
        for first_swap, second_swap in product((0, 1), repeat=2):
            masks = [0] * batches
            masks[first // 2] |= 1 << (2 * (first % 2) + first_swap)
            masks[second // 2] |= 1 << (2 * (second % 2) + second_swap)
            if first // 2 == second // 2:
                collisions += 1
                cross_masks[masks[first // 2]] += 1
            add_in_place(split_total, chronological_product(empty, matrices, masks))
            split_samples += 1

    require(pair_samples == 2 * slots, "full-pair sample count")
    require(split_samples == 4 * slots * (slots - 1), "distinct-pair sample count")
    require(F(collisions, split_samples) == F(1, 2 * batches - 1),
            "same-batch collision probability")
    require(set(full_masks) == set(FULL_PAIRS)
            and all(count * 2 == pair_samples for count in full_masks.values()),
            "uniform full-pair mask")
    require(set(cross_masks) == set(CROSS_PAIRS)
            and all(count * 4 == collisions for count in cross_masks.values()),
            "uniform collision mask")
    return (divide(pair_total, pair_samples), divide(split_total, split_samples),
            pair_samples, split_samples)


def exact_insertion(empty, event, steps, degree):
    total = zero()
    samples = 0
    for positions in combinations(range(steps), degree):
        selected = set(positions)
        term = identity()
        for step in range(steps):
            term = exact_multiply(term, event if step in selected else empty)
        add_in_place(total, term)
        samples += 1
    require(samples == comb(steps, degree), "insertion coefficient count")
    return total


def decimal_upper(matrix):
    return [[D(entry.numerator) / D(entry.denominator) for entry in row]
            for row in matrix]


def require_upper_close(actual, exact, label):
    strict = 0
    for i in range(2):
        for j in range(2):
            # Fraction(Decimal) preserves the exact returned decimal value.
            observed, reference = F(actual[i][j]), exact[i][j]
            require(observed >= reference,
                    ("directed bound below exact value", label, i, j, observed, reference))
            require(observed - reference <= max(F(1), reference) * RELATIVE_TOLERANCE,
                    ("directed bound is not close", label, i, j, observed, reference))
            strict += observed > reference
    return strict


def check_routes(matrices):
    fixtures = pair_samples = split_samples = strict = 0
    # lambda=1 also checks the collision normalization without gap decay.
    for contraction in (F(2, 3), F(1)):
        empty = [[F(1), F(0)], [F(0), contraction]]
        for batches in (2, 3, 4):
            pair_exact, split_exact, pair_count, split_count = exact_routes(batches, empty, matrices)
            with localcontext() as context:
                context.prec = PRECISION
                context.rounding = ROUND_CEILING
                pair_upper, split_upper = regional_matrices(
                    contraction,
                    decimal_upper(average(matrices, SINGLETONS)),
                    decimal_upper(average(matrices, FULL_PAIRS)),
                    decimal_upper(average(matrices, CROSS_PAIRS)),
                    steps=batches,
                )
            strict += require_upper_close(pair_upper, pair_exact, ("pair", batches, contraction))
            strict += require_upper_close(split_upper, split_exact, ("split", batches, contraction))
            pair_samples += pair_count
            split_samples += split_count
            fixtures += 1
    require(strict > 0, "rounding fixtures must exercise nonexact upper rounding")
    return fixtures, pair_samples, split_samples


def check_insertions(matrices):
    fixtures = 0
    for contraction in (F(2, 3), F(1)):
        empty = [[F(1), F(0)], [F(0), contraction]]
        for event in (matrices[3], average(matrices, SINGLETONS)):
            for steps in (2, 3, 4):
                for degree in (1, 2):
                    exact = exact_insertion(empty, event, steps, degree)
                    with localcontext() as context:
                        context.prec = PRECISION
                        context.rounding = ROUND_CEILING
                        actual = insertion_moment(decimal_upper(empty), decimal_upper(event), steps, degree)
                    require_upper_close(actual, exact, ("insertion", steps, degree, contraction))
                    fixtures += 1
    return fixtures


def exact_tail(tails, gap, contraction):
    ell = len(tails) - 1
    if gap < ell:
        return tails[gap]
    # Geometric extension acts on rows, before the endpoint matrix.
    return exact_multiply([[F(1), F(0)], [F(0), contraction**(gap - ell)]], tails[ell])


def exact_endpoint_regions(contraction, singleton, paired, cross, steps):
    def suffix(length):
        return [[F(1), F(0)], [F(0), contraction**length]]

    pair, collision, distinct = zero(), zero(), zero()
    for prefix in range(steps):
        after = suffix(steps - 1 - prefix)
        add_in_place(pair, exact_multiply(exact_tail(paired, prefix, contraction), after))
        add_in_place(collision, exact_multiply(exact_tail(cross, prefix, contraction), after))
    # Ordered gaps, including zero gaps; no exchange of the two endpoint maps.
    for first in range(steps - 1):
        for middle in range(steps - 1 - first):
            after = steps - 2 - first - middle
            term = exact_multiply(exact_tail(singleton, first, contraction),
                                  exact_tail(singleton, middle, contraction))
            add_in_place(distinct, exact_multiply(term, suffix(after)))
    pair = divide(pair, steps)
    collision = divide(collision, steps)
    distinct = divide(distinct, comb(steps, 2))
    probability = F(1, 2 * steps - 1)
    split = [[probability * collision[i][j] + (1 - probability) * distinct[i][j]
              for j in range(2)] for i in range(2)]
    return pair, split


def check_endpoint_routes(matrices):
    # Arbitrary noncommuting rational tails test the convolution algebra, not
    # the separate assertion that a tail table bounds an actual empty span.
    tables = (
        [average(matrices, SINGLETONS), matrices[3], matrices[12]],
        [average(matrices, FULL_PAIRS), matrices[1], matrices[9]],
        [average(matrices, CROSS_PAIRS), matrices[2], matrices[10]],
    )
    for table in tables:
        for left, right in combinations(table, 2):
            require(exact_multiply(left, right) != exact_multiply(right, left),
                    "endpoint toy tail matrices must not commute")
    convolution_fixtures = regional_fixtures = old_fixtures = 0
    for contraction in (F(2, 3), F(1)):
        for ell in (0, 1, 2):
            tails = [table[:ell + 1] for table in tables]
            for steps in (2, 3, 4, 5):
                previous = [matrices[MASKS[(index + 4) % len(MASKS)]] for index in range(steps)]
                for table in tails:
                    exact = []
                    for degree in range(steps):
                        total = zero()
                        for gap in range(degree + 1):
                            add_in_place(total, exact_multiply(exact_tail(table, gap, contraction),
                                                              previous[degree - gap]))
                        exact.append(total)
                    with localcontext() as context:
                        context.prec = PRECISION
                        context.rounding = ROUND_CEILING
                        decay = D(contraction.numerator) / D(contraction.denominator)
                        actual = gap_convolve([decimal_upper(matrix) for matrix in table],
                                              [decimal_upper(matrix) for matrix in previous], decay)
                    require(len(actual) == steps, "gap convolution length")
                    for degree, (observed, reference) in enumerate(zip(actual, exact)):
                        require_upper_close(observed, reference, ("gap", ell, steps, degree, contraction))
                    convolution_fixtures += 1

                exact = exact_endpoint_regions(contraction, *tails, steps)
                with localcontext() as context:
                    context.prec = PRECISION
                    context.rounding = ROUND_CEILING
                    actual = regions(contraction,
                                     *[[decimal_upper(matrix) for matrix in table] for table in tails],
                                     steps=steps)
                    if ell == 0:
                        old = regional_matrices(contraction,
                                                *(decimal_upper(table[0]) for table in tails), steps=steps)
                for kind, (observed, reference) in enumerate(zip(actual, exact)):
                    require_upper_close(observed, reference, ("endpoint region", ell, steps, kind, contraction))
                    if ell == 0:
                        require_upper_close(old[kind], reference, ("old region", steps, kind, contraction))
                        for i in range(2):
                            for j in range(2):
                                require(abs(F(observed[i][j]) - F(old[kind][i][j]))
                                        <= max(F(1), reference[i][j]) * RELATIVE_TOLERANCE,
                                        ("ell0 must recover old regional formula", steps, kind, i, j))
                old_fixtures += ell == 0
                regional_fixtures += 1
    return convolution_fixtures, regional_fixtures, old_fixtures


def main():
    matrices = toy_matrices()
    route_fixtures, full_samples, singleton_samples = check_routes(matrices)
    insertion_fixtures = check_insertions(matrices)
    convolutions, endpoints, old_regions = check_endpoint_routes(matrices)
    print(f"regional fixtures={route_fixtures}; full-pair slot/swap samples={full_samples}; "
          f"ordered singleton-slot/swap samples={singleton_samples}")
    print(f"exact insertion fixtures={insertion_fixtures}; Decimal precision={PRECISION}")
    print(f"endpoint convolution fixtures={convolutions}; endpoint regional fixtures={endpoints}; "
          f"ell0 old-formula fixtures={old_regions}")
    print("PASS: exact tiny routes, ordered insertions and endpoint gaps; not a distance certificate")


if __name__ == "__main__":
    main()
