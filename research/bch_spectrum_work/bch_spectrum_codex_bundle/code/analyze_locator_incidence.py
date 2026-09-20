#!/usr/bin/env python3
"""Probe the six-equation locator incidence variety at known minimum words.

For a 24-element set T in GF(256)^*, let R_T be the degree-at-most-23
remainder of u^146 modulo prod_{t in T}(u+t).  The set T is admissible exactly
when R_T has constant coefficient one and coefficients 19 through 23 equal
zero.  This script evaluates the formal Jacobian of those six coordinates at
subsets of every normalized locator discovered by the radius-five search.

The rank experiment is evidence about observed components.  It is not a point
count and is not a proof that undiscovered components are nonsingular.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import Counter
from pathlib import Path

from affine_wambach import ALPHA, MODULUS, gf_mul, gf_pow
from verify_locator_reduction import evaluate, locator_coefficients


ROOT = Path(__file__).resolve().parents[1]


def multiplication_table() -> list[bytes]:
    return [bytes(gf_mul(a, b) for b in range(256)) for a in range(256)]


def field_tables() -> tuple[list[int], dict[int, int]]:
    powers = [gf_pow(ALPHA, exponent) for exponent in range(255)]
    logs = {value: exponent for exponent, value in enumerate(powers)}
    assert len(logs) == 255
    return powers, logs


def support_from_word(word: int, powers: list[int]) -> list[int]:
    support = []
    for coordinate in range(255):
        if (word >> coordinate) & 1:
            support.append(powers[coordinate])
    if (word >> 255) & 1:
        support.append(0)
    assert len(support) == 38
    return support


def normalize_punctured_support(
    punctured: list[int], powers: list[int], logs: dict[int, int]
) -> list[int]:
    locator = locator_coefficients(punctured)
    leading_log = logs[locator[37]]
    # Scaling every support point by alpha^s multiplies the leading locator
    # coefficient by alpha^(37s).  Since gcd(37,255)=1, the shift is unique.
    shift = (-leading_log * pow(37, -1, 255)) % 255
    scale = powers[shift]
    normalized = [gf_mul(scale, x) for x in punctured]
    assert locator_coefficients(normalized)[37] == 1
    return normalized


def polynomial_product_roots(points: list[int], mul: list[bytes]) -> list[int]:
    coefficients = [1]
    for point in points:
        updated = [0] * (len(coefficients) + 1)
        row = mul[point]
        for degree, coefficient in enumerate(coefficients):
            updated[degree] ^= row[coefficient]
            updated[degree + 1] ^= coefficient
        coefficients = updated
    return coefficients


def synthetic_quotient(polynomial: list[int], root: int, mul: list[bytes]) -> list[int]:
    degree = len(polynomial) - 1
    quotient = [0] * degree
    quotient[-1] = polynomial[-1]
    row = mul[root]
    for index in range(degree - 2, -1, -1):
        quotient[index] = polynomial[index + 1] ^ row[quotient[index + 1]]
    assert polynomial[0] == row[quotient[0]]
    return quotient


def field_inverse(value: int, logs: dict[int, int], powers: list[int]) -> int:
    if value == 0:
        raise ZeroDivisionError
    return powers[-logs[value] % 255]


def matrix_rank(matrix: list[list[int]], mul: list[bytes], powers, logs) -> int:
    rows = [row[:] for row in matrix]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next(
            (index for index in range(rank, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = field_inverse(rows[rank][column], logs, powers)
        inverse_row = mul[inverse]
        rows[rank] = [inverse_row[value] for value in rows[rank]]
        for index in range(len(rows)):
            if index == rank or rows[index][column] == 0:
                continue
            factor_row = mul[rows[index][column]]
            rows[index] = [
                left ^ factor_row[right]
                for left, right in zip(rows[index], rows[rank])
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def derivative_value(g: list[int], point: int, mul, powers, logs) -> int:
    # Formal derivatives retain only odd-degree terms in characteristic two.
    value = 0
    for degree in range(1, len(g), 2):
        value ^= mul[g[degree]][gf_pow(point, degree - 1)]
    return value


def incidence_jacobian_rank(points: list[int], g: list[int], mul, powers, logs) -> int:
    root_polynomial = polynomial_product_roots(points, mul)
    selected_degrees = (0, 19, 20, 21, 22, 23)
    columns = []
    for point in points:
        lagrange_numerator = synthetic_quotient(root_polynomial, point, mul)
        denominator = evaluate(lagrange_numerator, point)
        scale = mul[derivative_value(g, point, mul, powers, logs)][
            field_inverse(denominator, logs, powers)
        ]
        scale_row = mul[scale]
        columns.append(
            [scale_row[lagrange_numerator[degree]] for degree in selected_degrees]
        )
    matrix = [
        [columns[column][row] for column in range(len(columns))]
        for row in range(6)
    ]
    return matrix_rank(matrix, mul, powers, logs)


def normalized_locators(orbit_path: Path, powers, logs) -> list[tuple[int, ...]]:
    data = json.loads(orbit_path.read_text(encoding="utf-8"))
    orbits = data["shells"]["38"]["orbits"]
    return normalized_locators_from_words(
        [int(orbit["representative_hex"], 16) for orbit in orbits], powers, logs
    )


def normalized_locators_from_words(
    words: list[int], powers, logs
) -> list[tuple[int, ...]]:
    locators = set()
    for word in words:
        support = support_from_word(word, powers)
        for origin in support:
            punctured = [point ^ origin for point in support if point != origin]
            normalized = normalize_punctured_support(punctured, powers, logs)
            locator = locator_coefficients(normalized)
            assert locator[0] == locator[37] == 1
            assert all(locator[degree] == 0 for degree in range(1, 37, 2))
            g = tuple(locator[2 * degree] for degree in range(19))
            locators.add(g)
    return sorted(locators)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subsets-per-locator", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0xB37)
    args = parser.parse_args()

    mul = multiplication_table()
    powers, logs = field_tables()
    orbit_path = ROOT / "generated" / "wambach_p_orbits_r5_w38.json"
    locators = set(normalized_locators(orbit_path, powers, logs))
    receipt_paths = [
        ROOT / "generated" / "frobenius_invariant_incidence.json",
        ROOT / "generated" / "frobenius4_invariant_incidence.json",
        ROOT / "generated" / "frobenius16_locator_extension_monte_carlo.json",
    ]
    receipts = [json.loads(path.read_text()) for path in receipt_paths]
    new_words = [
        int(record["canonical_word_hex"], 16)
        for receipt in receipts
        for record in receipt["new_affine_orbit_records"]
    ]
    locators.update(normalized_locators_from_words(new_words, powers, logs))
    locators = sorted(locators)
    assert len(locators) == 138 * 38

    rng = random.Random(args.seed)
    rank_histogram: Counter[int] = Counter()
    agreement_histogram: Counter[int] = Counter()
    derivative_zero_histogram: Counter[int] = Counter()
    degree_histogram: Counter[int] = Counter()
    for g_tuple in locators:
        g = list(g_tuple)
        degree_histogram[
            max(index for index, coefficient in enumerate(g) if coefficient)
        ] += 1
        agreements = [
            point
            for point in range(1, 256)
            if evaluate(g, point) == gf_pow(point, 146)
        ]
        # With z=u^128, these agreements are roots of
        # z^37 + g(z^2), a nonzero polynomial of degree 37.
        assert len(agreements) <= 37
        agreement_histogram[len(agreements)] += 1
        derivative_zero_histogram[
            sum(derivative_value(g, point, mul, powers, logs) == 0 for point in agreements)
        ] += 1
        for _ in range(args.subsets_per_locator):
            subset = rng.sample(agreements, 24)
            rank_histogram[
                incidence_jacobian_rank(subset, g, mul, powers, logs)
            ] += 1

    samples = len(locators) * args.subsets_per_locator
    result = {
        "classification": "diagnostic finite-field Jacobian experiment; not a point-count proof",
        "field_modulus": hex(MODULUS),
        "source_affine_orbits": 138,
        "source_radius_five_affine_orbits": 115,
        "source_new_frobenius_affine_orbits": 4,
        "source_new_frobenius4_affine_orbits": 10,
        "source_new_frobenius16_sample_affine_orbits": 9,
        "normalized_locator_polynomials": len(locators),
        "expected_locators_per_affine_orbit": 38,
        "subsets_per_locator": args.subsets_per_locator,
        "sampled_admissible_24_sets": samples,
        "incidence_equations": [
            "constant coefficient of R_T minus 1",
            "coefficients 19,20,21,22,23 of R_T",
        ],
        "R_T_definition": "u^146 mod product_{t in T}(u+t)",
        "jacobian_rank_histogram": {
            str(rank): count for rank, count in sorted(rank_histogram.items())
        },
        "agreement_count_histogram": {
            str(count): frequency
            for count, frequency in sorted(agreement_histogram.items())
        },
        "interpolant_degree_histogram": {
            str(degree): frequency
            for degree, frequency in sorted(degree_histogram.items())
        },
        "zero_derivatives_on_37_agreements_histogram": {
            str(count): frequency
            for count, frequency in sorted(derivative_zero_histogram.items())
        },
        "scope": (
            "Full rank at sampled known points shows that the six constraints "
            "are locally independent there. It does not exclude singular or "
            "high-point components elsewhere."
        ),
    }
    output = ROOT / "generated" / "locator_incidence_jacobian.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
