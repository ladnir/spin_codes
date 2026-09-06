#!/usr/bin/env python3
"""Check the GF(16) invariant-base sample and certify its endpoint witnesses."""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from affine_wambach import gf_pow
from analyze_locator_incidence import (
    multiplication_table,
    normalized_locators,
    normalized_locators_from_words,
)
from bch_quotient import binary_poly_divmod
from prepare_wambach_bases import (
    affine_canonical,
    generator_polynomial_generic,
    gf_pow_table,
)
from verify_locator_reduction import evaluate


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "generated" / "frobenius16_locator_extension_monte_carlo_raw.json"
OUTPUT = ROOT / "generated" / "frobenius16_locator_extension_monte_carlo.json"


def parse_polynomial(encoded: str) -> tuple[int, ...]:
    assert len(encoded) == 38
    return tuple(int(encoded[index : index + 2], 16) for index in range(0, 38, 2))


def in_subfield(polynomial: tuple[int, ...], size: int) -> bool:
    return all(gf_pow(coefficient, size) == coefficient for coefficient in polynomial)


def main() -> None:
    raw = json.loads(RAW.read_text(encoding="utf-8"))
    assert raw["samples"] == 10_000_000
    assert raw["field_orbit_size_histogram"] == {"1": 15, "2": 120}

    base_types = []
    base_population = 0
    for singletons in range(0, 16, 2):
        pairs = (18 - singletons) // 2
        population = math.comb(15, singletons) * math.comb(120, pairs)
        base_types.append((singletons, pairs, population))
        base_population += population
    target_population = sum(
        math.comb(15, singletons) * math.comb(120, (24 - singletons) // 2)
        for singletons in range(0, 16, 2)
    )
    assert base_population == raw["invariant_18_set_population"] == 199_417_781_755_185
    assert target_population == raw["invariant_24_set_population"] == 348_770_148_711_057_905
    assert sum(record["samples"] for record in raw["base_type_records"]) == raw["samples"]
    assert [
        (record["singletons"], record["pairs"], record["population"])
        for record in raw["base_type_records"]
    ] == base_types

    estimator = raw["weighted_invariant_extension_estimator"]
    assert 0.99 < estimator["ratio"] < 1.02
    assert estimator["normal_95_upper_ratio"] < 1.02
    endpoint_records = raw["distinct_endpoint_polynomials"]
    assert raw["endpoint_hits"] == sum(
        record["sample_multiplicity"] for record in endpoint_records
    ) == 10
    assert len(endpoint_records) == 10

    multiply = multiplication_table()
    powers, logs = gf_pow_table(8, 0x14D)
    inverses = [0] + [gf_pow(value, 254) for value in range(1, 256)]
    generator = generator_polynomial_generic(8, 37, 0x14D)
    sampled_orbits: defaultdict[int, list[tuple[tuple[int, ...], int, dict]]] = defaultdict(list)
    for record in endpoint_records:
        polynomial = parse_polynomial(record["coefficients_hex_ascending"])
        assert in_subfield(polynomial, 16)
        agreements = [
            point
            for point in range(1, 256)
            if evaluate(list(polynomial), point) == gf_pow(point, 146)
        ]
        assert len(agreements) == 37
        assert all(gf_pow(point, 16) in agreements for point in agreements)
        fixed_points = sum(gf_pow(point, 16) == point for point in agreements)
        assert fixed_points == record["agreement_fixed_points_in_GF16_star"]
        support = [gf_pow(point, 127) for point in agreements]
        word = (1 << 255) | sum(1 << logs[point] for point in support)
        assert word.bit_count() == 38
        assert binary_poly_divmod(word & ((1 << 255) - 1), generator)[1] == 0
        canonical, stabilizer = affine_canonical(
            word, powers, logs, 8, 0x14D, multiply, inverses
        )
        sampled_orbits[canonical].append((polynomial, stabilizer, record))
    assert len(sampled_orbits) == 9
    assert {records[0][1] for records in sampled_orbits.values()} == {1}

    wambach_path = ROOT / "generated" / "wambach_p_orbits_r5_w38.json"
    wambach = json.loads(wambach_path.read_text(encoding="utf-8"))
    known_words = {
        int(record["canonical_hex"], 16)
        for record in wambach["shells"]["38"]["orbits"]
    }
    for filename in (
        "frobenius_invariant_incidence.json",
        "frobenius4_invariant_incidence.json",
    ):
        receipt = json.loads((ROOT / "generated" / filename).read_text(encoding="utf-8"))
        known_words.update(
            int(record["canonical_word_hex"], 16)
            for record in receipt["new_affine_orbit_records"]
        )
    assert not set(sampled_orbits) & known_words

    prior_locators = set(normalized_locators(wambach_path, powers, logs))
    prior_locators.update(normalized_locators_from_words(list(known_words - {
        int(record["canonical_hex"], 16)
        for record in wambach["shells"]["38"]["orbits"]
    }), powers, logs))
    assert len(prior_locators) == 4_902
    prior_gf16_locators = {polynomial for polynomial in prior_locators if in_subfield(polynomial, 16)}
    assert len(prior_gf16_locators) == 126

    new_words = list(sampled_orbits)
    new_locators = set(normalized_locators_from_words(new_words, powers, logs))
    assert len(new_locators) == 9 * 38
    new_gf16_locators = {polynomial for polynomial in new_locators if in_subfield(polynomial, 16)}
    assert len(new_gf16_locators) == 64
    assert not prior_gf16_locators & new_gf16_locators
    per_orbit_gf16_count = {}
    for canonical in sampled_orbits:
        closure = normalized_locators_from_words([canonical], powers, logs)
        per_orbit_gf16_count[canonical] = sum(
            in_subfield(polynomial, 16) for polynomial in closure
        )
    assert sum(per_orbit_gf16_count.values()) == 64

    previous = json.loads(
        (ROOT / "generated" / "frobenius4_invariant_incidence.json").read_text(
            encoding="utf-8"
        )
    )["new_rigorous_shell_consequences"]
    additional_p_words = 9 * 65_280
    improved_p_lower = previous["improved_A38_P_lower"] + additional_p_words
    assert improved_p_lower == 9_008_640 and improved_p_lower % 255 == 0

    result = {
        "classification": {
            "incidence_count": "diagnostic Monte Carlo; not a proof or upper bound",
            "endpoint_records": "exact finite-field witnesses and affine closures",
        },
        "action": "x -> x^16 on GF(256)^*",
        "coefficient_subfield": "GF(16)",
        "field_orbit_size_histogram": {"1": 15, "2": 120},
        "samples": raw["samples"],
        "seed": raw["seed"],
        "runtime_seconds": raw["seconds"],
        "invariant_18_set_population": base_population,
        "invariant_24_set_population": target_population,
        "extra_agreement_histogram": raw["extra_agreement_histogram"],
        "weighted_invariant_extension_estimator": estimator,
        "estimator_definition": (
            "For an invariant 18-set B, sum 1/m(T) over invariant agreeing "
            "6-extensions T\\B, where m(T) is the number of invariant 18-subsets "
            "of T. Multiplication by the number of invariant 18-sets is unbiased "
            "for the number of admissible invariant 24-sets."
        ),
        "previously_certified_GF16_endpoint_inventory": {
            "affine_orbits_containing_GF16_locators": 21,
            "normalized_GF16_locators": len(prior_gf16_locators),
        },
        "sampled_endpoint_hits": raw["endpoint_hits"],
        "distinct_sampled_endpoint_polynomials": len(endpoint_records),
        "new_affine_orbits": len(sampled_orbits),
        "new_affine_orbit_sizes": [65_280] * len(sampled_orbits),
        "new_normalized_locators_from_affine_closure": len(new_locators),
        "new_GF16_endpoint_locators_from_affine_closure": len(new_gf16_locators),
        "certified_GF16_endpoint_inventory_after_closure": {
            "affine_orbits_containing_GF16_locators": 30,
            "normalized_GF16_locators": len(prior_gf16_locators | new_gf16_locators),
        },
        "new_affine_orbit_records": [
            {
                "canonical_word_hex": f"{canonical:064x}",
                "stabilizer_size": sampled_orbits[canonical][0][1],
                "orbit_size": 65_280,
                "sampled_GF16_polynomials_in_orbit": len(sampled_orbits[canonical]),
                "GF16_normalized_locators_in_affine_closure": per_orbit_gf16_count[
                    canonical
                ],
                "one_g_coefficient_vector_hex_ascending": [
                    f"{coefficient:02x}"
                    for coefficient in sampled_orbits[canonical][0][0]
                ],
            }
            for canonical in sorted(sampled_orbits)
        ],
        "new_rigorous_shell_consequences": {
            "additional_A38_P_words": additional_p_words,
            "improved_A38_P_lower": improved_p_lower,
            "improved_h38_lower": improved_p_lower // 255,
            "improved_A38_C_lower": 31 * improved_p_lower // 255,
        },
        "scope": (
            "The ratio estimates the GF(16)-invariant 24-set stratum only. "
            "The nine new affine orbits and their shell consequences are exact, "
            "but the sample does not exhaust the GF(16)-coefficient family."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
