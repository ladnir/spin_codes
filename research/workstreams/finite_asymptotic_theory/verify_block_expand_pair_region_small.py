#!/usr/bin/env python3
"""Verify the exact four-symbol kernel for one Block Expand region.

The proof formula is an exponential-generating-function coefficient.  This
script checks it against direct enumeration for small region lengths and
small input pair types.  It is a formula check, not a length-512 certificate.
"""

from __future__ import annotations

import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "block_expand_pair_region_small_exact.json"
TYPE4 = tuple[int, int, int, int]
EXP3 = tuple[int, int, int]
POLY = dict[EXP3, Fraction]


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first,) + tail


def multiply(first: POLY, second: POLY, caps: EXP3) -> POLY:
    result: POLY = {}
    for (i1, j1, k1), c1 in first.items():
        for (i2, j2, k2), c2 in second.items():
            exponent = (i1 + i2, j1 + j2, k1 + k2)
            if all(exponent[index] <= caps[index] for index in range(3)):
                result[exponent] = result.get(exponent, Fraction(0)) + c1 * c2
    return result


def power(base: POLY, exponent: int, caps: EXP3) -> POLY:
    result: POLY = {(0, 0, 0): Fraction(1)}
    for _ in range(exponent):
        result = multiply(result, base, caps)
    return result


def local_egf(parity: tuple[int, int], caps: EXP3) -> POLY:
    """Return the truncated EGF for one bin with the requested output parity."""
    result: POLY = {}
    for i in range(caps[0] + 1):
        for j in range(caps[1] + 1):
            for k in range(caps[2] + 1):
                if ((i + k) & 1, (j + k) & 1) == parity:
                    result[(i, j, k)] = Fraction(
                        1, math.factorial(i) * math.factorial(j) * math.factorial(k)
                    )
    return result


def egf_count(length: int, inputs: EXP3, outputs: TYPE4) -> int:
    """Count labeled edge assignments using the four local EGFs."""
    polynomial: POLY = {(0, 0, 0): Fraction(1)}
    parities = ((0, 0), (0, 1), (1, 0), (1, 1))
    for parity, multiplicity in zip(parities, outputs):
        polynomial = multiply(
            polynomial,
            power(local_egf(parity, inputs), multiplicity, inputs),
            inputs,
        )

    choose_bins = math.factorial(length)
    for multiplicity in outputs:
        choose_bins //= math.factorial(multiplicity)
    labeled = polynomial.get(inputs, Fraction(0)) * choose_bins
    for count in inputs:
        labeled *= math.factorial(count)
    if labeled.denominator != 1:
        raise AssertionError("EGF assignment count is not integral")
    return labeled.numerator


def brute_counts(length: int, inputs: EXP3) -> Counter[TYPE4]:
    symbols = [(1, 0)] * inputs[0] + [(0, 1)] * inputs[1] + [(1, 1)] * inputs[2]
    result: Counter[TYPE4] = Counter()
    for locations in itertools.product(range(length), repeat=len(symbols)):
        bins = [[0, 0] for _ in range(length)]
        for symbol, location in zip(symbols, locations):
            bins[location][0] ^= symbol[0]
            bins[location][1] ^= symbol[1]
        counts = Counter(map(tuple, bins))
        result[(counts[(0, 0)], counts[(0, 1)], counts[(1, 0)], counts[(1, 1)])] += 1
    return result


def main() -> None:
    cases = []
    tested_entries = 0
    for length in range(1, 5):
        for inputs in itertools.product(range(4), repeat=3):
            if sum(inputs) > 5:
                continue
            brute = brute_counts(length, inputs)
            formula = {
                outputs: egf_count(length, inputs, outputs)
                for outputs in compositions(length, 4)
            }
            formula = Counter({key: value for key, value in formula.items() if value})
            if formula != brute:
                raise AssertionError(
                    f"regional pair kernel mismatch: length={length}, inputs={inputs}"
                )
            total = length ** sum(inputs)
            if sum(formula.values()) != total:
                raise AssertionError("regional assignment counts do not sum correctly")
            cases.append(
                {
                    "region_length": length,
                    "input_nonzero_type": list(inputs),
                    "assignment_count": total,
                    "nonzero_output_types": len(formula),
                }
            )
            tested_entries += math.comb(length + 3, 3)

    receipt = {
        "schema": "block-expand-pair-region-small-exact-v1",
        "status": "pass",
        "symbol_order": ["00", "01", "10", "11"],
        "input_nonzero_type_order": ["10", "01", "11"],
        "formula": (
            "p!q!r! * multinomial(ell;m00,m10,m01,m11) * "
            "[X^p Y^q Z^r] product_ab F_ab(X,Y,Z)^m_ab"
        ),
        "local_egf": (
            "F_ab=(1/4) sum_{s,t in {+1,-1}} "
            "s^a t^b exp(sX+tY+stZ)"
        ),
        "region_lengths_tested": [1, 2, 3, 4],
        "maximum_nonzero_input_coordinates": 5,
        "cases_tested": len(cases),
        "output_type_entries_tested": tested_entries,
        "notes": [
            "Every equality is exact over Python integers and rational numbers.",
            "The denominator ell^(p+q+r) converts assignment counts to probabilities.",
            "This receipt verifies the regional formula only; it is not a length-512 shell-cap certificate.",
        ],
        "cases": cases,
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: receipt[key] for key in ("status", "cases_tested", "output_type_entries_tested")}, indent=2))


if __name__ == "__main__":
    main()
