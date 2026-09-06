#!/usr/bin/env python3
"""Enumerate all GF(2)-coefficient perturbations of x^146 over GF(256)."""

from __future__ import annotations

import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

from affine_wambach import MODULUS, gf_pow


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "generated" / "binary_coefficient_146_family.json"


def main() -> None:
    points = np.arange(1, 256, dtype=np.uint8)
    target = np.asarray([gf_pow(int(point), 146) for point in points], dtype=np.uint8)
    monomials = np.asarray(
        [
            [gf_pow(int(point), degree) for point in points]
            for degree in range(1, 19)
        ],
        dtype=np.uint8,
    )

    values = np.ones(255, dtype=np.uint8)
    previous_gray = 0
    histogram: Counter[int] = Counter()
    locator_masks: list[int] = []
    for index in range(1 << 18):
        gray = index ^ (index >> 1)
        if index:
            toggled_degree = (gray ^ previous_gray).bit_length() - 1
            values ^= monomials[toggled_degree]
        agreements = int(np.count_nonzero(values == target))
        histogram[agreements] += 1
        if agreements == 37:
            locator_masks.append(gray)
        previous_gray = gray

    assert sum(histogram.values()) == 1 << 18
    assert max(histogram) == 37
    assert len(locator_masks) == 10
    assert not any(histogram[count] for count in range(29, 37))

    incidence = sum(
        frequency * math.comb(agreements, 24)
        for agreements, frequency in histogram.items()
    )
    random_reference = Fraction(math.comb(255, 24), 256**6)
    ratio = Fraction(incidence, 1) / random_reference
    endpoint_ratio = Fraction(
        len(locator_masks) * math.comb(37, 24), 1
    ) / random_reference
    result = {
        "classification": "exact exhaustive finite-field computation",
        "field_modulus": hex(MODULUS),
        "family": (
            "g(x)=1+sum_{j=1}^{18} c_j*x^j with every c_j in GF(2), "
            "tested against x^146 on GF(256)^*"
        ),
        "family_size": 1 << 18,
        "agreement_count_histogram": {
            str(agreements): frequency
            for agreements, frequency in sorted(histogram.items())
        },
        "maximum_agreements": max(histogram),
        "polynomials_with_at_least_24_agreements": sum(
            frequency for agreements, frequency in histogram.items() if agreements >= 24
        ),
        "locator_polynomials_with_37_agreements": len(locator_masks),
        "locator_coefficient_masks_c1_through_c18_hex": [
            f"0x{mask:05x}" for mask in sorted(locator_masks)
        ],
        "gap_below_endpoint": "no polynomial has 29 through 36 agreements",
        "I24_contribution": incidence,
        "R6_contribution": {
            "numerator": str(ratio.numerator),
            "denominator": str(ratio.denominator),
            "decimal": float(ratio),
            "log2": math.log2(float(ratio)),
        },
        "endpoint_R6_contribution": {
            "numerator": str(endpoint_ratio.numerator),
            "denominator": str(endpoint_ratio.denominator),
            "decimal": float(endpoint_ratio),
        },
        "interpretation": (
            "This Frobenius-fixed coefficient family is visibly non-random and "
            "contains endpoint locators, but its complete contribution to R6 is "
            "negligible relative to the certified factor-6 target."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "maximum_agreements": max(histogram),
        "endpoint_locators": len(locator_masks),
        "R6_contribution": float(ratio),
    }, indent=2))


if __name__ == "__main__":
    main()
