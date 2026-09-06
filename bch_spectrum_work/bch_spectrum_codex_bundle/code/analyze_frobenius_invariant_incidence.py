#!/usr/bin/env python3
"""Enumerate all 24-sets fixed by x -> x^2 and their incidences."""

from __future__ import annotations

import itertools
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from affine_wambach import gf_pow
from analyze_locator_incidence import (
    multiplication_table,
    polynomial_product_roots,
)
from bch_quotient import binary_poly_divmod
from prepare_wambach_bases import (
    affine_canonical,
    generator_polynomial_generic,
    gf_pow_table,
)
from verify_locator_reduction import evaluate
from verify_symmetric_incidence import complete_homogeneous, monomial_remainder


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "generated" / "frobenius_invariant_incidence.json"


def main() -> None:
    multiply = multiplication_table()
    seen: set[int] = set()
    orbits: list[tuple[int, ...]] = []
    for value in range(1, 256):
        if value in seen:
            continue
        orbit = []
        current = value
        while current not in orbit:
            orbit.append(current)
            seen.add(current)
            current = multiply[current][current]
        orbits.append(tuple(orbit))
    orbit_histogram = Counter(map(len, orbits))
    assert orbit_histogram == {1: 1, 2: 1, 4: 3, 8: 30}

    degree_four = [orbit for orbit in orbits if len(orbit) == 4]
    degree_eight = [orbit for orbit in orbits if len(orbit) == 8]
    candidates: list[tuple[str, tuple[int, ...]]] = []
    candidates.extend(
        ("8+8+8", sum(selected, ()))
        for selected in itertools.combinations(degree_eight, 3)
    )
    candidates.extend(
        ("8+8+4+4", sum(eight + four, ()))
        for eight in itertools.combinations(degree_eight, 2)
        for four in itertools.combinations(degree_four, 2)
    )
    assert len(candidates) == 5365

    admissible_by_type: Counter[str] = Counter()
    set_agreement_histogram: Counter[int] = Counter()
    polynomial_subset_multiplicity: defaultdict[tuple[int, ...], int] = defaultdict(int)
    polynomial_agreements: dict[tuple[int, ...], int] = {}
    for orbit_type, points_tuple in candidates:
        points = list(points_tuple)
        homogeneous = complete_homogeneous(points, 127, multiply)
        product = 1
        for point in points:
            product = multiply[product][point]
        if multiply[product][homogeneous[122]] != 1 or any(
            homogeneous[index] for index in range(123, 128)
        ):
            continue

        modulus = polynomial_product_roots(points, multiply)
        remainder = monomial_remainder(146, modulus, multiply)
        assert remainder[0] == 1 and not any(remainder[19:])
        polynomial = tuple(remainder[:19])
        if polynomial not in polynomial_agreements:
            polynomial_agreements[polynomial] = sum(
                evaluate(list(polynomial), point) == gf_pow(point, 146)
                for point in range(1, 256)
            )
        agreements = polynomial_agreements[polynomial]
        assert agreements >= 24
        admissible_by_type[orbit_type] += 1
        set_agreement_histogram[agreements] += 1
        polynomial_subset_multiplicity[polynomial] += 1

    distinct_polynomial_agreements = Counter(polynomial_agreements.values())
    distinct_polynomial_degrees = Counter(
        max(index for index, coefficient in enumerate(polynomial) if coefficient)
        for polynomial in polynomial_agreements
    )
    locator_polynomials = {
        polynomial
        for polynomial, agreements in polynomial_agreements.items()
        if agreements == 37
    }
    assert len(locator_polynomials) == 10

    powers, logs = gf_pow_table(8, 0x14D)
    inverses = [0] * 256
    for value in range(1, 256):
        inverses[value] = gf_pow(value, 254)
    p_generator = generator_polynomial_generic(8, 37, 0x14D)
    locator_orbits: defaultdict[int, list[tuple[int, ...]]] = defaultdict(list)
    locator_orbit_stabilizers: dict[int, int] = {}
    for polynomial in locator_polynomials:
        agreements = [
            point
            for point in range(1, 256)
            if evaluate(list(polynomial), point) == gf_pow(point, 146)
        ]
        support = [gf_pow(point, 127) for point in agreements]
        word = (1 << 255) | sum(1 << logs[point] for point in support)
        assert word.bit_count() == 38
        assert binary_poly_divmod(word & ((1 << 255) - 1), p_generator)[1] == 0
        canonical, stabilizer = affine_canonical(
            word, powers, logs, 8, 0x14D, multiply, inverses
        )
        locator_orbits[canonical].append(polynomial)
        locator_orbit_stabilizers[canonical] = stabilizer

    known_orbit_data = json.loads(
        (ROOT / "generated" / "wambach_p_orbits_r5_w38.json").read_text()
    )
    known_canonical_orbits = {
        int(record["canonical_hex"], 16)
        for record in known_orbit_data["shells"]["38"]["orbits"]
    }
    new_canonical_orbits = set(locator_orbits) - known_canonical_orbits
    new_p_words = sum(
        256 * 255 // locator_orbit_stabilizers[canonical]
        for canonical in new_canonical_orbits
    )

    reference = Fraction(math.comb(255, 24), 256**6)
    contribution_ratio = Fraction(sum(admissible_by_type.values()), 1) / reference
    result = {
        "classification": "exact exhaustive Frobenius-orbit computation",
        "action": "x -> x^2 on GF(256)^*",
        "field_orbit_size_histogram": {
            str(size): count for size, count in sorted(orbit_histogram.items())
        },
        "invariant_24_sets": len(candidates),
        "candidate_type_histogram": {
            "8+8+8": math.comb(30, 3),
            "8+8+4+4": math.comb(30, 2) * math.comb(3, 2),
        },
        "admissible_24_sets": sum(admissible_by_type.values()),
        "admissible_type_histogram": dict(sorted(admissible_by_type.items())),
        "set_agreement_count_histogram": {
            str(count): frequency
            for count, frequency in sorted(set_agreement_histogram.items())
        },
        "distinct_admissible_polynomials": len(polynomial_agreements),
        "distinct_polynomial_agreement_count_histogram": {
            str(count): frequency
            for count, frequency in sorted(distinct_polynomial_agreements.items())
        },
        "distinct_polynomial_degree_histogram": {
            str(degree): frequency
            for degree, frequency in sorted(distinct_polynomial_degrees.items())
        },
        "invariant_subset_multiplicity_per_polynomial_histogram": {
            str(count): frequency
            for count, frequency in sorted(
                Counter(polynomial_subset_multiplicity.values()).items()
            )
        },
        "locator_polynomials_with_37_agreements": len(locator_polynomials),
        "distinct_affine_orbits_represented_by_locators": len(locator_orbits),
        "locator_polynomials_per_affine_orbit_histogram": {
            str(count): frequency
            for count, frequency in sorted(
                Counter(map(len, locator_orbits.values())).items()
            )
        },
        "affine_orbits_already_known": len(set(locator_orbits) & known_canonical_orbits),
        "new_affine_orbits": len(new_canonical_orbits),
        "new_affine_orbit_sizes": sorted(
            256 * 255 // locator_orbit_stabilizers[canonical]
            for canonical in new_canonical_orbits
        ),
        "new_affine_orbit_records": [
            {
                "canonical_word_hex": f"{canonical:064x}",
                "stabilizer_size": locator_orbit_stabilizers[canonical],
                "orbit_size": 256 * 255 // locator_orbit_stabilizers[canonical],
                "frobenius_invariant_normalized_locators_in_orbit": len(
                    locator_orbits[canonical]
                ),
                "one_g_coefficient_vector_hex_ascending": [
                    f"{coefficient:02x}"
                    for coefficient in locator_orbits[canonical][0]
                ],
            }
            for canonical in sorted(new_canonical_orbits)
        ],
        "new_rigorous_shell_consequences": {
            "additional_A38_P_words": new_p_words,
            "improved_A38_P_lower": int(known_orbit_data["shells"]["38"]["rigorous_lower_bound"])
            + new_p_words,
            "improved_h38_lower": (
                int(known_orbit_data["shells"]["38"]["rigorous_lower_bound"])
                + new_p_words
            )
            // 255,
            "improved_A38_C_lower": 31
            * (
                int(known_orbit_data["shells"]["38"]["rigorous_lower_bound"])
                + new_p_words
            )
            // 255,
        },
        "R6_contribution": {
            "numerator": str(contribution_ratio.numerator),
            "denominator": str(contribution_ratio.denominator),
            "decimal": float(contribution_ratio),
            "log2": math.log2(float(contribution_ratio)),
        },
        "interpretation": (
            "Frobenius-invariant sets are strongly enriched relative to an "
            "independent rank model, but their entire absolute contribution is "
            "negligible on the global R6 scale."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
