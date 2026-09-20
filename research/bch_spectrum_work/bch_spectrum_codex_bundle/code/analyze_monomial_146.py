#!/usr/bin/env python3
"""Compute exact differential and Walsh spectra of x -> x^146 on GF(256)."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from affine_wambach import MODULUS, gf_mul, gf_pow


ROOT = Path(__file__).resolve().parents[1]


def absolute_trace(value: int, multiply: list[bytes]) -> int:
    result = 0
    conjugate = value
    for _ in range(8):
        result ^= conjugate
        conjugate = multiply[conjugate][conjugate]
    assert result in (0, 1)
    return result


def main() -> None:
    multiply = [bytes(gf_mul(a, b) for b in range(256)) for a in range(256)]
    trace = bytes(absolute_trace(value, multiply) for value in range(256))
    monomial = bytes(gf_pow(value, 146) for value in range(256))

    differential_histogram: Counter[int] = Counter()
    differential_uniformity = 0
    for difference in range(1, 256):
        fibers = Counter(
            monomial[value ^ difference] ^ monomial[value]
            for value in range(256)
        )
        differential_histogram.update(fibers.values())
        differential_uniformity = max(differential_uniformity, max(fibers.values()))

    walsh_histogram: Counter[int] = Counter()
    maximum_pairs: list[tuple[int, int]] = []
    maximum_absolute_walsh = 0
    for output_mask in range(1, 256):
        masked_output = bytes(
            multiply[output_mask][monomial[value]] for value in range(256)
        )
        for input_mask in range(256):
            input_row = multiply[input_mask]
            coefficient = sum(
                1 - 2 * trace[masked_output[value] ^ input_row[value]]
                for value in range(256)
            )
            walsh_histogram[coefficient] += 1
            absolute = abs(coefficient)
            if absolute > maximum_absolute_walsh:
                maximum_absolute_walsh = absolute
                maximum_pairs = [(output_mask, input_mask)]
            elif absolute == maximum_absolute_walsh:
                maximum_pairs.append((output_mask, input_mask))

    inverse_exponent = pow(146, -1, 255)
    expected_maximum_pairs = [
        (output_mask, gf_pow(output_mask, inverse_exponent))
        for output_mask in range(1, 256)
    ]
    assert sorted(maximum_pairs) == expected_maximum_pairs
    assert differential_uniformity == 6
    assert maximum_absolute_walsh == 64

    result = {
        "classification": "exact exhaustive finite-field computation",
        "field_modulus": hex(MODULUS),
        "monomial_exponent": 146,
        "binary_weight_of_exponent": 3,
        "inverse_exponent_mod_255": inverse_exponent,
        "differential_uniformity": differential_uniformity,
        "nonzero_differential_fiber_size_histogram": {
            str(size): count
            for size, count in sorted(differential_histogram.items())
        },
        "maximum_absolute_walsh_coefficient": maximum_absolute_walsh,
        "walsh_coefficient_histogram": {
            str(coefficient): count
            for coefficient, count in sorted(walsh_histogram.items())
        },
        "maximum_walsh_pairs": (
            "exactly (a,b)=(a,a^131) for every a in GF(256)^*"
        ),
        "scope": (
            "This is a building-block spectrum for a character-sum proof. "
            "It does not by itself bound the 24-point incidence variety."
        ),
    }
    output = ROOT / "generated" / "monomial_146_spectrum.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
