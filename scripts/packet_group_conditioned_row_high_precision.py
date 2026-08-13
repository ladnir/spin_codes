"""Independent high-precision conditioned-row outer reference.

This module intentionally does not import the outward evaluator. It rebuilds
the packet moments, split-spectrum sums, and graph average with Decimal
arithmetic. The result is a numerical cross-check, not a directed interval.
"""

from __future__ import annotations

import csv
import math
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


ROWS = 64
BAND_SIZES = (42, 43, 43)
NORMAL_TILES = 128
HOLE_TILES = 128
GRAPH_DIMENSION = 24


def decimal_fraction(value: Fraction) -> Decimal:
    return Decimal(value.numerator) / Decimal(value.denominator)


def polynomial_multiply(left: list[Decimal], right: list[Decimal]) -> list[Decimal]:
    result = [Decimal(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i + j] += x * y
    return result


def polynomial_power(base: list[Decimal], exponent: int) -> list[Decimal]:
    result = [Decimal(1)]
    factor = base
    power = exponent
    while power:
        if power & 1:
            result = polynomial_multiply(result, factor)
        power >>= 1
        if power:
            factor = polynomial_multiply(factor, factor)
    return result


def load_spectrum(path: Path) -> dict[tuple[int, int], int]:
    result = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["band0_weight"]), int(row["band1_weight"]))
            count = int(row["count"])
            if count:
                result[key] = count
    return result


def load_graph_spectrum(path: Path) -> dict[int, int]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            int(row["weight"]): int(row["count"])
            for row in csv.DictReader(handle)
            if int(row["count"])
        }


def log2_sum_exp(terms: list[Decimal], ln2: Decimal) -> Decimal:
    maximum = max(terms)
    total = sum(((value - maximum) * ln2).exp() for value in terms)
    return maximum + total.ln() / ln2


def evaluate_conditioned_row_high_precision(
    group_bits: int,
    profile: list[int],
    log_variables: list[float],
    coefficient_values: float | list[float],
    theta_value: float,
    spectrum01_path: Path,
    spectrum12_path: Path,
    graph_spectrum_path: Path,
    *,
    precision: int = 180,
) -> tuple[Decimal, dict[str, int]]:
    """Evaluate the exact-dyadic frozen parameters with high precision."""

    if group_bits not in (1, 2, 4, 8) or 64 % group_bits:
        raise ValueError("high-precision reference: unsupported packet width")
    if len(profile) != group_bits + 1 or len(log_variables) != group_bits + 1:
        raise ValueError("high-precision reference: dimension mismatch")

    spectrum01 = load_spectrum(spectrum01_path)
    spectrum12 = load_spectrum(spectrum12_path)
    graph_spectrum = load_graph_spectrum(graph_spectrum_path)
    punctured01 = {
        (a, b): (42 - a) * spectrum01.get((a, b), 0)
        + (a + 1) * spectrum01.get((a + 1, b), 0)
        for a in range(42)
        for b in range(44)
    }
    punctured01 = {key: count for key, count in punctured01.items() if count}
    masses = {
        "spectrum01": sum(spectrum01.values()),
        "spectrum12": sum(spectrum12.values()),
        "punctured01": sum(punctured01.values()),
        "graph": sum(graph_spectrum.values()),
    }
    expected = {
        "spectrum01": 1 << 64,
        "spectrum12": 1 << 64,
        "punctured01": 42 * (1 << 64),
        "graph": 1 << GRAPH_DIMENSION,
    }
    if masses != expected or spectrum01.get((0, 0)) != 1 or spectrum12.get((0, 0)) != 1:
        raise ValueError("high-precision reference: source spectrum mass changed")

    with localcontext() as context:
        context.prec = precision
        ln2 = Decimal(2).ln()
        variables = [
            decimal_fraction(Fraction.from_float(math.exp(float(value))))
            for value in log_variables
        ]
        atom = [
            Decimal(math.comb(group_bits, weight)) * variables[weight]
            for weight in range(group_bits + 1)
        ]
        coefficients = polynomial_power(atom, 64 // group_bits)
        moments = [
            coefficient / Decimal(math.comb(64, weight))
            for weight, coefficient in enumerate(coefficients)
        ]

        if isinstance(coefficient_values, list):
            band_coefficients = [
                decimal_fraction(Fraction.from_float(float(value)))
                for value in coefficient_values
            ]
        else:
            band1 = decimal_fraction(Fraction.from_float(float(coefficient_values)))
            band_coefficients = [Decimal(1) - band1, band1, band1]
        theta = decimal_fraction(Fraction.from_float(float(theta_value)))

        def conditioned_norm(coefficient: Decimal, bit: int) -> Decimal:
            terms = [
                Decimal(math.comb(63, weight)).ln() / ln2
                - Decimal(63)
                + moments[weight + bit].ln() / (coefficient * ln2)
                for weight in range(64)
            ]
            return coefficient * log2_sum_exp(terms, ln2)

        factors = [
            [conditioned_norm(coefficient, bit) for bit in (0, 1)]
            for coefficient in band_coefficients
        ]
        zeros = [pair[0] for pair in factors]
        ratios = [pair[1] - pair[0] for pair in factors]

        def pair_log_enumerator(
            table: dict[tuple[int, int], int],
            left: Decimal,
            right: Decimal,
            divisor: int,
        ) -> Decimal:
            return log2_sum_exp(
                [
                    Decimal(Fraction(count, divisor).numerator).ln() / ln2
                    - Decimal(Fraction(count, divisor).denominator).ln() / ln2
                    + Decimal(a) * left
                    + Decimal(b) * right
                    for (a, b), count in table.items()
                ],
                ln2,
            )

        left01 = Decimal(2) * ratios[0]
        middle01 = Decimal(2) * theta * ratios[1]
        middle12 = Decimal(2) * (Decimal(1) - theta) * ratios[1]
        right12 = Decimal(2) * ratios[2]

        def cauchy(first: dict[tuple[int, int], int], divisor: int) -> Decimal:
            return (
                pair_log_enumerator(first, left01, middle01, divisor)
                + pair_log_enumerator(spectrum12, middle12, right12, 1)
            ) / Decimal(2)

        normal_pair = cauchy(spectrum01, 1)
        punctured_pair = cauchy(punctured01, 42)
        free_message = Decimal((ROWS - 1) * ROWS)
        normal_tile = (
            free_message
            + sum(Decimal(size) * zero for size, zero in zip(BAND_SIZES, zeros))
            + normal_pair
        )
        punctured_common = (
            free_message
            + Decimal(41) * zeros[0]
            + Decimal(43) * zeros[1]
            + Decimal(43) * zeros[2]
            + punctured_pair
        )
        hole_tiles = [punctured_common + factors[0][bit] for bit in (0, 1)]
        graph = log2_sum_exp(
            [
                Decimal(count).ln() / ln2
                - Decimal(GRAPH_DIMENSION)
                + Decimal(HOLE_TILES - weight) * hole_tiles[0]
                + Decimal(weight) * hole_tiles[1]
                for weight, count in graph_spectrum.items()
            ],
            ln2,
        )
        constant = Decimal(NORMAL_TILES) * normal_tile + graph
        charges = [value.ln() / ln2 for value in variables]
        value = constant - sum(
            Decimal(count) * charge for count, charge in zip(profile, charges)
        )
        return +value, masses
