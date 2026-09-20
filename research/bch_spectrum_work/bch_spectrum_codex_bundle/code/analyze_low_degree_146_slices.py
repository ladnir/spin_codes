#!/usr/bin/env python3
"""Exactly enumerate low-degree perturbations of x^146 over GF(256)."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np

from affine_wambach import MODULUS, gf_mul, gf_pow


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "generated" / "low_degree_146_slices.json"
Q = 256


def absolute_trace(value: int, multiply: np.ndarray) -> int:
    result = 0
    conjugate = value
    for _ in range(8):
        result ^= conjugate
        conjugate = int(multiply[conjugate, conjugate])
    assert result in (0, 1)
    return result


def decode_coefficients(index: int, degree: int) -> list[str]:
    return [f"0x{(index >> (8 * offset)) & 255:02x}" for offset in range(degree)]


def enumerate_degree(
    degree: int, multiply: np.ndarray, target: np.ndarray
) -> dict[str, object]:
    """Count roots of x^146 + 1 + sum_j c_j x^j for all coefficients."""

    family_size = Q**degree
    prefix_size = Q ** (degree - 1)
    counts = np.zeros(family_size, dtype=np.uint8)
    prefixes = np.arange(prefix_size, dtype=np.uint64)
    prefix_coefficients = [
        ((prefixes >> (8 * offset)) & 255).astype(np.uint8)
        for offset in range(degree - 1)
    ]

    for value in range(1, Q):
        residual = np.full(prefix_size, target[value], dtype=np.uint8)
        power = value
        for coefficient in prefix_coefficients:
            residual ^= multiply[coefficient, power]
            power = gf_mul(power, value)
        inverse_last_power = gf_pow(power, 254)
        last_coefficient = multiply[inverse_last_power, residual]
        indices = prefixes + (
            last_coefficient.astype(np.uint64) << (8 * (degree - 1))
        )
        # The low bytes of each index equal its prefix, so indices are distinct.
        counts[indices] += 1

    histogram = np.bincount(counts, minlength=Q)
    maximum = int(np.max(counts))
    maximizing = np.flatnonzero(counts == maximum)
    assert int(np.sum(histogram)) == family_size
    assert sum(root_count * int(number) for root_count, number in enumerate(histogram)) == (
        255 * prefix_size
    )
    return {
        "degree_at_most": degree,
        "family_size": family_size,
        "maximum_number_of_roots": maximum,
        "number_attaining_maximum": int(len(maximizing)),
        "maximizing_coefficient_vectors_c1_up": [
            decode_coefficients(int(index), degree) for index in maximizing
        ],
        "root_count_histogram": {
            str(root_count): int(number)
            for root_count, number in enumerate(histogram)
            if number
        },
    }


def main() -> None:
    multiply = np.asarray(
        [[gf_mul(left, right) for right in range(Q)] for left in range(Q)],
        dtype=np.uint8,
    )
    target = np.asarray([gf_pow(value, 146) ^ 1 for value in range(Q)], dtype=np.uint8)
    slices = [enumerate_degree(degree, multiply, target) for degree in (1, 2, 3)]

    exceptional_by_linear_coefficient: Counter[int] = Counter()
    exceptional_total = 0
    for linear_coefficient in range(Q):
        signed_trace_sum = sum(
            1 - 2 * absolute_trace(output_mask, multiply)
            for output_mask in range(1, Q)
            if gf_pow(output_mask, 130) == linear_coefficient
        )
        contribution = 64 * signed_trace_sum
        exceptional_by_linear_coefficient[contribution] += 1
        exceptional_total += contribution
    assert exceptional_total == -64

    result = {
        "classification": "exact exhaustive finite-field computation",
        "field_modulus": hex(MODULUS),
        "monomial_exponent": 146,
        "polynomial_family": "x^146 + 1 + sum_{j=1}^r c_j x^j",
        "slices": slices,
        "locator_consequence": (
            "Every perturbation of degree at most 3 has fewer than 18 roots, "
            "so these slices contribute zero to the 18-base and 24-point incidences."
        ),
        "exceptional_walsh_ridge": {
            "condition": "b=a^131; for b=a*c this is c=a^130",
            "signed_contribution_histogram_over_c": {
                str(contribution): count
                for contribution, count in sorted(exceptional_by_linear_coefficient.items())
            },
            "signed_total_over_all_256_linear_coefficients": exceptional_total,
            "explanation": (
                "The exceptional numerator contribution for c is "
                "64*sum_{a:a^130=c}(-1)^Tr(a). Summing over c counts every "
                "nonzero a once and therefore gives 64*sum_{a!=0} psi(a)=-64."
            ),
        },
        "scope": (
            "The ordinary Walsh ridge controls the affine slice. Higher-degree "
            "perturbations require generalized exponential sums, so this receipt "
            "does not bound the full degree-18 locator family."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "slice_maxima": {
            str(row["degree_at_most"]): row["maximum_number_of_roots"]
            for row in slices
        },
        "exceptional_ridge_signed_total": exceptional_total,
    }, indent=2))


if __name__ == "__main__":
    main()
