#!/usr/bin/env python3
"""Check and summarize the exhaustive x -> x^4 invariant 24-set run."""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from affine_wambach import gf_pow
from analyze_locator_incidence import multiplication_table
from bch_quotient import binary_poly_divmod
from prepare_wambach_bases import (
    affine_canonical,
    generator_polynomial_generic,
    gf_pow_table,
)
from verify_locator_reduction import evaluate


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "generated" / "frobenius4_invariant_incidence_exact_raw.json"
OUTPUT = ROOT / "generated" / "frobenius4_invariant_incidence.json"


def parse_polynomial(encoded: str) -> tuple[int, ...]:
    assert len(encoded) == 38
    return tuple(int(encoded[index : index + 2], 16) for index in range(0, 38, 2))


def main() -> None:
    raw = json.loads(RAW.read_text(encoding="utf-8"))
    assert raw["exhaustive"] is True

    cases = []
    total_population = 0
    for singles in range(4):
        for pairs in range(7):
            remaining = 24 - singles - 2 * pairs
            if remaining < 0 or remaining % 4:
                continue
            quads = remaining // 4
            if quads > 60:
                continue
            population = (
                math.comb(3, singles)
                * math.comb(6, pairs)
                * math.comb(60, quads)
            )
            cases.append((singles, pairs, quads, population))
            total_population += population
    assert cases == [
        (0, 0, 6, 50_063_860),
        (0, 2, 5, 81_922_680),
        (0, 4, 4, 7_314_525),
        (0, 6, 3, 34_220),
        (2, 1, 5, 98_307_216),
        (2, 3, 4, 29_258_100),
        (2, 5, 3, 615_960),
    ]
    assert total_population == raw["samples"] == raw["invariant_24_set_population"]
    assert all(record["samples"] == record["population"] for record in raw["candidate_types"])

    set_histogram = Counter(
        {int(agreements): count for agreements, count in raw["agreement_count_histogram"].items()}
    )
    assert sum(set_histogram.values()) == raw["admissible_hits"] == 67_900
    polynomial_records = raw["distinct_polynomials"]
    polynomial_histogram = Counter(record["agreements"] for record in polynomial_records)
    multiplicity_histogram = Counter(record["sample_multiplicity"] for record in polynomial_records)
    assert sum(record["sample_multiplicity"] for record in polynomial_records) == 67_900
    assert all(
        sum(
            record["sample_multiplicity"]
            for record in polynomial_records
            if record["agreements"] == agreements
        )
        == count
        for agreements, count in set_histogram.items()
    )
    # Every x -> x^4 invariant set of size at least 24 contains an invariant
    # 24-subset: check the finite orbit-count statement directly. Hence the
    # set enumeration discovers every GF(4)-coefficient polynomial with N>=24.
    assert all(
        any(
            selected_singles + 2 * selected_pairs + 4 * selected_quads == 24
            for selected_singles in range(singles + 1)
            for selected_pairs in range(pairs + 1)
            for selected_quads in range(quads + 1)
        )
        for singles in range(4)
        for pairs in range(7)
        for quads in range(61)
        if singles + 2 * pairs + 4 * quads >= 24
    )
    coefficient_family_i24 = sum(
        math.comb(record["agreements"], 24) for record in polynomial_records
    )
    endpoint_i24 = len(
        [record for record in polynomial_records if record["agreements"] == 37]
    ) * math.comb(37, 24)

    multiply = multiplication_table()
    powers, logs = gf_pow_table(8, 0x14D)
    inverses = [0] + [gf_pow(value, 254) for value in range(1, 256)]
    p_generator = generator_polynomial_generic(8, 37, 0x14D)
    endpoint_records = [record for record in polynomial_records if record["agreements"] == 37]
    assert len(endpoint_records) == 34
    endpoint_orbits: defaultdict[int, list[tuple[tuple[int, ...], int, int]]] = defaultdict(list)
    endpoint_subset_multiplicity = Counter()
    for record in endpoint_records:
        polynomial = parse_polynomial(record["coefficients_hex_ascending"])
        agreements = [
            point
            for point in range(1, 256)
            if evaluate(list(polynomial), point) == gf_pow(point, 146)
        ]
        assert len(agreements) == 37
        assert all(gf_pow(point, 4) in agreements for point in agreements)
        support = [gf_pow(point, 127) for point in agreements]
        word = (1 << 255) | sum(1 << logs[point] for point in support)
        assert word.bit_count() == 38
        assert binary_poly_divmod(word & ((1 << 255) - 1), p_generator)[1] == 0
        canonical, stabilizer = affine_canonical(
            word, powers, logs, 8, 0x14D, multiply, inverses
        )
        endpoint_orbits[canonical].append(
            (polynomial, stabilizer, record["sample_multiplicity"])
        )
        endpoint_subset_multiplicity[record["sample_multiplicity"]] += 1
    assert len(endpoint_orbits) == 15
    assert Counter(map(len, endpoint_orbits.values())) == {2: 13, 4: 2}
    assert endpoint_subset_multiplicity == {84: 26, 364: 8}
    assert {records[0][1] for records in endpoint_orbits.values()} == {1}

    wambach = json.loads(
        (ROOT / "generated" / "wambach_p_orbits_r5_w38.json").read_text(encoding="utf-8")
    )
    known_canonical = {
        int(record["canonical_hex"], 16)
        for record in wambach["shells"]["38"]["orbits"]
    }
    binary = json.loads(
        (ROOT / "generated" / "frobenius_invariant_incidence.json").read_text(
            encoding="utf-8"
        )
    )
    known_canonical.update(
        int(record["canonical_word_hex"], 16)
        for record in binary["new_affine_orbit_records"]
    )
    new_canonical = set(endpoint_orbits) - known_canonical
    assert len(set(endpoint_orbits) & known_canonical) == 5
    assert len(new_canonical) == 10
    new_p_words = sum(
        256 * 255 // endpoint_orbits[canonical][0][1]
        for canonical in new_canonical
    )
    prior_p_lower = binary["new_rigorous_shell_consequences"]["improved_A38_P_lower"]
    improved_p_lower = prior_p_lower + new_p_words
    assert improved_p_lower % 255 == 0

    independent_reference = Fraction(total_population, 4**6)
    conditional_ratio = Fraction(raw["admissible_hits"], 1) / independent_reference
    global_reference = Fraction(math.comb(255, 24), 256**6)
    contribution_ratio = Fraction(raw["admissible_hits"], 1) / global_reference
    coefficient_family_ratio = Fraction(coefficient_family_i24, 1) / global_reference
    endpoint_family_ratio = Fraction(endpoint_i24, 1) / global_reference
    result = {
        "classification": "exact exhaustive Frobenius-orbit computation",
        "action": "x -> x^4 on GF(256)^*",
        "coefficient_subfield": "GF(4)",
        "field_orbit_size_histogram": {"1": 3, "2": 6, "4": 60},
        "invariant_24_sets": total_population,
        "candidate_type_records": [
            {
                "single_orbits": singles,
                "pair_orbits": pairs,
                "quad_orbits": quads,
                "sets": population,
                "admissible_sets": raw_record["hits"],
            }
            for (singles, pairs, quads, population), raw_record in zip(
                cases, raw["candidate_types"], strict=True
            )
        ],
        "admissible_24_sets": raw["admissible_hits"],
        "set_agreement_count_histogram": {
            str(agreements): count for agreements, count in sorted(set_histogram.items())
        },
        "distinct_admissible_polynomials": len(polynomial_records),
        "GF4_coefficient_family_I24_contribution": coefficient_family_i24,
        "GF4_coefficient_family_R6_contribution": {
            "numerator": str(coefficient_family_ratio.numerator),
            "denominator": str(coefficient_family_ratio.denominator),
            "decimal": float(coefficient_family_ratio),
            "log2": math.log2(float(coefficient_family_ratio)),
        },
        "GF4_endpoint_R6_contribution": {
            "numerator": str(endpoint_family_ratio.numerator),
            "denominator": str(endpoint_family_ratio.denominator),
            "decimal": float(endpoint_family_ratio),
            "fraction_of_GF4_family_I24": endpoint_i24 / coefficient_family_i24,
        },
        "distinct_polynomial_agreement_count_histogram": {
            str(agreements): count
            for agreements, count in sorted(polynomial_histogram.items())
        },
        "invariant_subset_multiplicity_per_polynomial_histogram": {
            str(multiplicity): count
            for multiplicity, count in sorted(multiplicity_histogram.items())
        },
        "conditional_GF4_rank_model_ratio": {
            "numerator": str(conditional_ratio.numerator),
            "denominator": str(conditional_ratio.denominator),
            "decimal": float(conditional_ratio),
        },
        "invariant_set_R6_contribution": {
            "numerator": str(contribution_ratio.numerator),
            "denominator": str(contribution_ratio.denominator),
            "decimal": float(contribution_ratio),
            "log2": math.log2(float(contribution_ratio)),
        },
        "locator_polynomials_with_37_agreements": len(endpoint_records),
        "locator_invariant_subset_multiplicity_histogram": {
            str(multiplicity): count
            for multiplicity, count in sorted(endpoint_subset_multiplicity.items())
        },
        "distinct_affine_orbits_represented_by_locators": len(endpoint_orbits),
        "locator_polynomials_per_affine_orbit_histogram": {
            str(count): frequency
            for count, frequency in sorted(Counter(map(len, endpoint_orbits.values())).items())
        },
        "affine_orbits_already_known": len(set(endpoint_orbits) & known_canonical),
        "new_affine_orbits": len(new_canonical),
        "new_affine_orbit_sizes": sorted(
            256 * 255 // endpoint_orbits[canonical][0][1]
            for canonical in new_canonical
        ),
        "new_affine_orbit_records": [
            {
                "canonical_word_hex": f"{canonical:064x}",
                "stabilizer_size": endpoint_orbits[canonical][0][1],
                "orbit_size": 256 * 255 // endpoint_orbits[canonical][0][1],
                "frobenius4_invariant_normalized_locators_in_orbit": len(
                    endpoint_orbits[canonical]
                ),
                "one_g_coefficient_vector_hex_ascending": [
                    f"{coefficient:02x}"
                    for coefficient in endpoint_orbits[canonical][0][0]
                ],
            }
            for canonical in sorted(new_canonical)
        ],
        "new_rigorous_shell_consequences": {
            "additional_A38_P_words": new_p_words,
            "improved_A38_P_lower": improved_p_lower,
            "improved_h38_lower": improved_p_lower // 255,
            "improved_A38_C_lower": 31 * improved_p_lower // 255,
        },
        "kernel": {
            "identity": (
                "h_n=sum_{0<=i<=24, i congruent n mod 2} "
                "e_i*h_((n-i)/2)^2"
            ),
            "needed_h_indices": "1..31, 49..63, and 122..127",
            "batch_width": 32,
            "raw_runtime_seconds": raw["seconds"],
            "checks": (
                "opening batches compared with a scalar GF(4) recurrence; every hit "
                "was independently reconstructed as x^146 mod P_T and evaluated"
            ),
        },
        "interpretation": (
            "The invariant-set equations are only mildly enriched over six "
            "independent GF(4) equations. The full GF(4)-coefficient family's "
            "contribution is negligible in the global R6 sum and is 99.79 percent "
            "endpoint-generated. Its locators nevertheless yield ten new exact "
            "affine orbits."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
