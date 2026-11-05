#!/usr/bin/env python3
"""Bounded exact GF16 checks of the one-diagonal forward kernel.

This is a local algebra/counting regression, not a distance certificate.
It uses only the standard library, writes no files, and should run with -B.
The root agent owns execution; intended cap: 30 seconds, 256 MiB.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb

from test_persistent_mds import INV, MUL, Q, add, inverse, matvec, rank, require


B = 2
D = tuple(tuple(INV[x ^ y] for y in (2, 3)) for x in (0, 1))
VECTORS = tuple(product(range(Q), repeat=B))
MARKERS = tuple(product((Fraction(1, 2), Fraction(3, 4), Fraction(1)),
                       (Fraction(1, 4), Fraction(1, 2), Fraction(1))))


def weight(vector):
    return sum(value != 0 for value in vector)


def step(state, input_vector, key):
    output = tuple(x ^ MUL[r][v] for x, r, v in zip(input_vector, key, state))
    return output, matvec(D, add(state, output))


def exact_support_count(length):
    if length == 0:
        return 1
    inclusion = sum((-1) ** j * comb(length, j) * Q ** max(0, length - B - j)
                    for j in range(length + 1))
    standard = (0 if length <= B else
                (Q - 1) * sum((-1) ** j * comb(length - 1, j)
                              * Q ** (length - B - 1 - j)
                              for j in range(length - B)))
    require(inclusion == standard >= 0, ("MDS coefficient", length, inclusion, standard))
    return inclusion


def split_coefficients(h):
    coefficients = Counter({(0, 0): 1})
    for m in range(1, h + 1):
        for outgoing in range(B + 1 - m, B + 1):
            coefficients[m, outgoing] = (comb(h, m) * comb(B, outgoing)
                                         * exact_support_count(m + outgoing))
    return coefficients


def evaluate(coefficients, z, s):
    return sum((count * z ** m * s ** outgoing
                for (m, outgoing), count in coefficients.items()), Fraction(0))


def sharper_count(m, zero_count):
    terms = [Fraction(comb(B, j), comb(zero_count, j)) * Q ** (m - j)
             for j in range(min(m, zero_count) + 1)]
    require(min(terms) == terms[-1], ("q>=b monotonicity", m, zero_count, terms))
    return min(terms)


def check_superregularity():
    cases = 0
    for size in range(1, B + 1):
        for rows in combinations(range(B), size):
            for columns in combinations(range(B), size):
                require(rank(tuple(tuple(D[i][j] for j in columns) for i in rows)) == size,
                        ("Cauchy minor", rows, columns))
                cases += 1
    return cases


def check_shells():
    fixtures = (
        ((1, 0), (0, 0)), ((1, 0), (1, 0)), ((1, 0), (0, 1)),
        ((0, 1), (0, 0)), ((0, 1), (0, 7)), ((0, 2), (8, 0)),
        ((1, 1), (0, 0)), ((1, 1), (1, 1)), ((2, 3), (4, 5)),
        ((7, 9), (3, 0)), ((15, 1), (0, 14)), ((6, 10), (11, 12)),
    )
    comparisons = potential_checks = 0
    for state, input_vector in fixtures:
        active = tuple(j for j, value in enumerate(state) if value)
        h = len(active)
        shells, weight_histogram, reset_weights = Counter(), Counter(), Counter()
        for key in VECTORS:
            output, outgoing = step(state, input_vector, key)
            m = sum(output[j] != 0 for j in active)
            h_next = weight(outgoing)
            shells[m, h_next] += 1
            weight_histogram[weight(output), h_next] += 1
            if h_next == 0:
                reset_weights[weight(output)] += 1
        for m in range(h + 1):
            for h_next in range(1, B + 1):
                count = shells[m, h_next]
                plain = (Q ** (B - h) * comb(h, m) * comb(B, h_next)
                         * Q ** max(0, m - B + h_next))
                sharp = Q ** (B - h) * comb(h, m) * sharper_count(m, B - h_next)
                require(count <= plain and count <= sharp,
                        ("joint shell", state, input_vector, m, h_next, count, plain, sharp))
                comparisons += 1
        feasible = all(input_vector[j] == 0 for j in range(B) if j not in active)
        expected_reset = Q ** (B - h) if feasible else 0
        require(sum(reset_weights.values()) == expected_reset,
                ("exact reset probability", state, input_vector, reset_weights))
        require(not reset_weights or set(reset_weights) == {h}, "reset output weight")
        for z, s in MARKERS:
            actual = sum((count * z ** w * s ** (h_next - B)
                          for (w, h_next), count in weight_histogram.items()), Fraction(0)) / Q ** B
            bound = Fraction(1, Q ** h) * s ** (-B) * evaluate(split_coefficients(h), z, s)
            require(actual <= bound, ("potential", state, input_vector, z, s, actual, bound))
            potential_checks += 1
    return len(fixtures), len(fixtures) * len(VECTORS), comparisons, potential_checks


def check_fourier_and_split():
    offsets = ((0, 0), (1, 0), (0, 1), (1, 1), (2, 3), (15, 7))
    split_checks = affine_checks = 0
    for h in range(1, B + 1):
        for active in combinations(range(B), h):
            inputs = tuple(product(range(Q), repeat=h))
            g = tuple(tuple(row[j] for j in active) for row in D)
            homogeneous = Counter((weight(u), weight(matvec(g, u))) for u in inputs)
            expected = +split_coefficients(h)
            require(homogeneous == expected, ("split enumerator", active, homogeneous, expected))
            split_checks += 1
            for offset in offsets:
                affine = Counter((weight(u), weight(add(matvec(g, u), offset))) for u in inputs)
                for z, s in MARKERS:
                    require(evaluate(affine, z, s) <= evaluate(homogeneous, z, s),
                            ("affine Fourier maximum", active, offset, z, s))
                    affine_checks += 1
    return split_checks, affine_checks


def check_activation_and_two_batch_reset():
    activation_checks = 0
    for first in VECTORS:
        output, outgoing = step((0, 0), first, (7, 11))
        require(output == first and outgoing == matvec(D, first), "quiet transition")
        if any(first):
            require(weight(outgoing) >= B + 1 - weight(first), ("activation MDS", first))
        activation_checks += 1

    # One sparse activation and one activation whose transformed state is sparse.
    first_inputs = ((1, 0), matvec(inverse(D), (1, 0)))
    require({weight(matvec(D, first)) for first in first_inputs} == {1, 2},
            "two-batch fixtures must exercise both incoming support weights")
    reset_checks = coefficient_assignments = 0
    for first in first_inputs:
        state = matvec(D, first)
        h = weight(state)
        for second in VECTORS:
            count = 0
            for key in VECTORS:
                output = tuple(x ^ MUL[r][v] for x, r, v in zip(second, key, state))
                if output == state:
                    require(weight(output) == h, "two-batch reset weight")
                    count += 1
                coefficient_assignments += 1
            feasible = all(x == 0 for x, v in zip(second, state) if v == 0)
            require(count == (Q ** (B - h) if feasible else 0),
                    ("two-batch exact reset", first, second, count))
            reset_checks += 1
    return activation_checks, reset_checks, coefficient_assignments


def main():
    minors = check_superregularity()
    fixtures, assignments, shell_checks, potential_checks = check_shells()
    split_checks, affine_checks = check_fourier_and_split()
    activations, reset_checks, reset_assignments = check_activation_and_two_batch_reset()
    print("PASS forward-shell toy interfaces:", dict(
        superregular_minors=minors, fixed_input_fixtures=fixtures,
        shell_key_assignments=assignments, joint_shell_checks=shell_checks,
        potential_checks=potential_checks, split_enumerators=split_checks,
        affine_fourier_checks=affine_checks, activation_checks=activations,
        two_batch_reset_inputs=reset_checks, two_batch_key_assignments=reset_assignments))


if __name__ == "__main__":
    main()
