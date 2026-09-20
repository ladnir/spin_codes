#!/usr/bin/env python3
"""Verify the Kasami--Tokura enumerator below twice the RM distance."""

from __future__ import annotations

import csv
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
REPOSITORY = WORKSTREAM.parent.parent
RM49_SPECTRUM = REPOSITORY / "scripts" / "rm512_256_spectrum.csv"
OUTPUT = HERE / "rm511_low_weight_exact.json"


def gaussian_binomial(n: int, k: int) -> int:
    if not 0 <= k <= n:
        return 0
    numerator = 1
    denominator = 1
    for index in range(k):
        numerator *= (1 << (n - index)) - 1
        denominator *= (1 << (k - index)) - 1
    if numerator % denominator:
        raise AssertionError("nonintegral Gaussian binomial")
    return numerator // denominator


def family_two(m: int, r: int, h: int) -> int:
    odd_product = 1
    for exponent in range(3, 2 * h, 2):
        odd_product *= (1 << exponent) - 1
    return (
        (1 << (r + h * h + h - 2))
        * gaussian_binomial(m, r - 2)
        * gaussian_binomial(m - r + 2, 2 * h)
        * odd_product
    )


def family_three(m: int, r: int, h: int) -> int:
    return (
        (1 << (r + h * h + h - 1))
        * gaussian_binomial(m, r - h)
        * gaussian_binomial(m - r, h)
        * gaussian_binomial(m - r + h, h)
    )


def canonical_counts(m: int, r: int) -> dict[int, int]:
    """Apply the four cases of the Kasami--Tokura formula."""
    distance = 1 << (m - r)
    a = min(r, m - r)
    b = 1 + (m - r) / 2
    result = {distance: (1 << r) * gaussian_binomial(m, m - r)}
    for h in range(2, int(max(a, b)) + 1):
        weight = 2 * distance - (1 << (m - r + 1 - h))
        if h == 2 or max(a, 2) < h <= b:
            count = family_two(m, r, h)
        elif max(b, 2) < h <= a:
            count = family_three(m, r, h)
        elif 2 < h <= min(a, b):
            count = family_two(m, r, h) + family_three(m, r, h)
        else:
            raise AssertionError(f"uncovered formula case h={h}")
        result[weight] = count
    return result


def main() -> int:
    with RM49_SPECTRUM.open(newline="", encoding="utf-8") as source:
        rm49 = {int(row["weight"]): int(row["count"]) for row in csv.DictReader(source)}
    expected_rm49 = {weight: count for weight, count in rm49.items() if 0 < weight < 64}
    computed_rm49 = canonical_counts(9, 4)
    if computed_rm49 != expected_rm49:
        raise AssertionError("formula does not reproduce the exact RM(4,9) low spectrum")

    rm511 = canonical_counts(11, 5)
    payload = {
        "schema": "rm511-low-weight-enumerator-v1",
        "status": "EXACT_INTEGER_FORMULA_CHECK",
        "code": {"name": "RM(5,11)", "length": 2048, "dimension": 1024, "distance": 64},
        "scope": "all nonzero coefficients below twice the minimum distance",
        "coefficients": {str(weight): str(count) for weight, count in sorted(rm511.items())},
        "regression": {
            "code": "RM(4,9)",
            "matched_weights": sorted(computed_rm49),
            "matches_authenticated_full_spectrum": True,
        },
        "source": {
            "result": "Kasami--Tokura formula for weights below twice the minimum distance",
            "doi": "10.1109/TIT.1970.1054545",
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print(f"output={OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
