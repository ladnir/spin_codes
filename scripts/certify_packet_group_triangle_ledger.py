#!/usr/bin/env python3
"""Outward verifier for a finite ``g=2`` triangular profile ledger.

The triangle-cover program is a discovery/combinatorial component.  This
module consumes only its selected fixed witnesses and its actually used
witness-at-vertex assignments.  Each distinct witness is recomputed once with
the existing exact/outward shared-drive machinery; only referenced affine
profile inequalities are then evaluated.

The accepted ledger schema is deliberately small, while allowing a few aliases
so the discovery artifact can settle independently::

    {
      "group_bits": 2,
      "complete_cover": true,
      "triangles": [
        {"vertices": [[a0,a1,a2], ...], "witness": "row-name"},
        {"vertices": [[...], ...],
         "vertex_witnesses": ["same-name", "same-name", "same-name"]},
        {"vertices": [[...], ...],
         "mixture_witness_references": ["source.json:0", "source.json:3"],
         "mixture_weights": ["1/3", "2/3"]}
      ],
      "terminal_triangles": [
        {"vertices": [[a0,a1,a2], ...],
         "terminal_integer_lattice_count": 3,
         "terminal_integer_profiles_sha256": "..."}
      ],
      "discrete_point_assignments": [
        {"profile": [a0,a1,a2], "witness_reference": "source.json:7"}
      ],
      "complete_global_integer_cover": true
    }

An already flattened ``used_inequalities``/``inequalities``/
``vertex_assignments`` array of ``{"profile": [...], "witness": "..."}``
records is also accepted.  The verifier intentionally does not rerun witness
selection or edit/reinterpret the triangulation.  In the hybrid form,
``triangles + terminal_triangles`` must partition the exact hull.  Covered
triangles use the affine vertex argument; every lattice profile of every
terminal triangle must instead occur in ``discrete_point_assignments``.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import itertools
import json
import math
import os
import re
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from certify_packet8_hard_face_drive_inner import exact_float, outward_transport_image
from certify_packet_group_profile import outward_normalization
from outward_log2 import LN2, Interval, log2_fraction, log2_int, self_check
from packet_group_drive_stratified import (
    accumulator_pattern_counts,
    block_histograms,
    compatible_pair_polynomials,
    compatible_pair_polynomials_outward,
    point_caps,
    profile_count,
)
from packet_group_outer_profile import (
    BLOCK_BITS,
    D,
    GROUPS,
    K,
    N,
    atom_count,
    full_spectrum,
)
from packet_group_profile_bound import GRAPH_REPLACEMENTS, split_cap_table
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_repeated_value_state_transfer import best_witness
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS
from certify_three_band_exact_length import load_graph_spectrum


GROUP_BITS = 2
BLOCK_ATOMS = BLOCK_BITS // GROUP_BITS
DECIMAL_GUARD = 150
TRANSCENDENTAL_PAD = Decimal("1e-120")
EDGE_NAME = re.compile(r"^edge_(\d+)_(\d+)_")


_WORKER_SPLIT_CAPS = None
_OUTWARD_NORMALIZATION_CACHE: dict[tuple[int, int, int], Interval] = {}


def _up_add(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    raw = left + right
    return np.where(right > 0.0, np.nextafter(raw, np.inf), left)


def _up_mul(left: np.ndarray | float, right: np.ndarray | float) -> np.ndarray:
    left_array = np.asarray(left)
    right_array = np.asarray(right)
    raw = left_array * right_array
    positive_underflow = (left_array > 0.0) & (right_array > 0.0) & (raw == 0.0)
    if np.any(positive_underflow):
        raise FloatingPointError(
            "outward multiplication underflowed a structurally positive product"
        )
    return np.where(raw > 0.0, np.nextafter(raw, np.inf), raw)


def _clean_outward_pairs(fugacities: np.ndarray) -> np.ndarray:
    exact_support = compatible_pair_polynomials(GROUP_BITS, fugacities)
    upper = compatible_pair_polynomials_outward(GROUP_BITS, fugacities)
    upper[exact_support == 0.0] = 0.0
    return upper


def _block_histograms_outward(fugacities: np.ndarray) -> np.ndarray:
    """Existing outward block DP with exact structural zeros retained."""

    atoms = BLOCK_ATOMS
    patterns = accumulator_pattern_counts(GROUP_BITS)
    pairs = _clean_outward_pairs(fugacities)
    entries = []
    for incoming in range(2):
        for outgoing in range(2):
            for selected in range(GROUP_BITS + 1):
                for drive in range(GROUP_BITS + 1):
                    for emitted in range(GROUP_BITS + 1):
                        count = int(patterns[incoming, outgoing, drive, emitted])
                        coefficient = float(_up_mul(float(count), pairs[drive, selected]))
                        if coefficient:
                            entries.append(
                                (incoming, outgoing, selected, drive, emitted, coefficient)
                            )
    dp = np.zeros((2, BLOCK_BITS + 1, BLOCK_BITS + 1, BLOCK_BITS + 1))
    dp[0, 0, 0, 0] = 1.0
    maximum = 0
    for _ in range(atoms):
        next_dp = np.zeros_like(dp)
        for incoming, outgoing, selected, drive, emitted, coefficient in entries:
            source = dp[incoming, : maximum + 1, : maximum + 1, : maximum + 1]
            contribution = _up_mul(coefficient, source)
            destination = next_dp[
                outgoing,
                selected : selected + maximum + 1,
                drive : drive + maximum + 1,
                emitted : emitted + maximum + 1,
            ]
            destination[...] = _up_add(destination, contribution)
        dp = next_dp
        maximum += GROUP_BITS
    return _up_add(dp[0], dp[1])


def _convolve_outward(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.zeros(len(left) + len(right) - 1)
    for index, value in enumerate(left):
        if value == 0.0:
            continue
        contribution = _up_mul(value, right)
        target = result[index : index + len(right)]
        target[...] = _up_add(target, contribution)
    return result


def _point_caps_outward(fugacities: np.ndarray) -> np.ndarray:
    atoms = BLOCK_ATOMS
    local = _clean_outward_pairs(fugacities)
    maxima = np.zeros((BLOCK_BITS + 1, BLOCK_BITS + 1))
    for weights in itertools.combinations_with_replacement(
        range(GROUP_BITS + 1), atoms
    ):
        polynomial = np.array([1.0])
        for weight in weights:
            polynomial = _convolve_outward(polynomial, local[weight])
        total = sum(weights)
        maxima[:, total] = np.maximum(maxima[:, total], polynomial)
    for state_weight in range(BLOCK_BITS + 1):
        positive = maxima[state_weight] > 0.0
        quotient = maxima[state_weight, positive] / math.comb(BLOCK_BITS, state_weight)
        maxima[state_weight, positive] = np.nextafter(quotient, np.inf)
    return maxima


def _decimal_unary_bounds(value: Decimal, operation: str) -> Interval:
    """Guarded bounds for Decimal exp/ln, matching outward_log2's policy."""

    if operation == "ln" and value <= 0:
        raise ValueError("logarithm requires a positive value")
    with localcontext() as context:
        context.prec = DECIMAL_GUARD
        middle = value.exp() if operation == "exp" else value.ln()
        pad = TRANSCENDENTAL_PAD * max(Decimal(1), abs(middle))
        return Interval(middle - pad, middle + pad)


def _exp_interval(value: Interval) -> Interval:
    return Interval(
        _decimal_unary_bounds(value.lo, "exp").lo,
        _decimal_unary_bounds(value.hi, "exp").hi,
    )


def _ln_interval(value: Interval) -> Interval:
    if value.lo <= 0:
        raise ValueError("interval logarithm requires a positive lower endpoint")
    return Interval(
        _decimal_unary_bounds(value.lo, "ln").lo,
        _decimal_unary_bounds(value.hi, "ln").hi,
    )


def _ln_sum_exp(terms: Iterable[Interval]) -> Interval:
    rows = tuple(terms)
    if not rows:
        raise ValueError("log-sum-exp requires a nonempty input")
    shift = max(row.hi for row in rows)
    total = Interval.exact(0)
    shift_interval = Interval.exact(shift)
    for row in rows:
        total = total + _exp_interval(row - shift_interval)
    return _ln_interval(total) + shift_interval


def _log2_sum_exp(terms: Iterable[Interval]) -> Interval:
    """Outward ``log2(sum(2**term))`` for base-two log intervals."""

    rows = tuple(terms)
    if not rows:
        raise ValueError("base-two log-sum-exp requires a nonempty input")
    return _ln_sum_exp(row * LN2 for row in rows) / LN2


def _fraction_polynomial_power(base: tuple[Fraction, ...], power: int) -> list[Fraction]:
    result = [Fraction(1)]
    for _ in range(power):
        product = [Fraction(0)] * (len(result) + len(base) - 1)
        for left, left_value in enumerate(result):
            for right, right_value in enumerate(base):
                product[left + right] += left_value * right_value
        result = product
    return result


def _fraction_polynomial_multiply(
    left: Iterable[Fraction], right: Iterable[Fraction]
) -> list[Fraction]:
    left_rows = tuple(left)
    right_rows = tuple(right)
    result = [Fraction(0)] * (len(left_rows) + len(right_rows) - 1)
    for left_index, left_value in enumerate(left_rows):
        for right_index, right_value in enumerate(right_rows):
            result[left_index + right_index] += left_value * right_value
    return result


def _linear_outer(
    log_variables: list[float], band1_coefficient: float
) -> tuple[Interval, tuple[Interval, ...], dict[str, Any]]:
    """Outward affine linear-BL outer witness from frozen binary64 parameters."""

    if len(log_variables) != GROUP_BITS + 1:
        raise ValueError("linear-BL witness has the wrong variable dimension")
    variables = tuple(exact_float(math.exp(float(value))) for value in log_variables)
    if any(value <= 0 for value in variables):
        raise ValueError("linear-BL variables must be positive")
    band1 = exact_float(float(band1_coefficient))
    band0 = Fraction(1) - band1
    if not Fraction(0) < band0 <= Fraction(1, 2) or not Fraction(1, 2) <= band1 < 1:
        raise ValueError("linear-BL band coefficient is outside [1/2,1)")

    local = tuple(
        Fraction(math.comb(GROUP_BITS, weight)) * variables[weight]
        for weight in range(GROUP_BITS + 1)
    )
    coefficients = _fraction_polynomial_power(local, BLOCK_ATOMS)
    moments = tuple(
        coefficient / math.comb(BLOCK_BITS, weight)
        for weight, coefficient in enumerate(coefficients)
    )

    def band_norm(coefficient: Fraction) -> Interval:
        terms = []
        coefficient_interval = Interval.exact(coefficient.numerator) / Interval.exact(
            coefficient.denominator
        )
        for weight, moment in enumerate(moments):
            probability = Fraction(math.comb(BLOCK_BITS, weight), 1 << BLOCK_BITS)
            terms.append(
                (log2_fraction(probability) * LN2)
                + (log2_fraction(moment) * LN2) / coefficient_interval
            )
        return coefficient_interval * _ln_sum_exp(terms) / LN2

    norm0 = band_norm(band0)
    norm1 = band_norm(band1)
    constant = (
        Interval.exact(K)
        + norm0.times_int(256 * 42)
        + norm1.times_int(256 * 86)
    )
    replacement_ratio = max(
        variables[new] / variables[old]
        for old in range(GROUP_BITS + 1)
        for new in range(GROUP_BITS + 1)
        if abs(new - old) <= 1
    )
    graph = log2_fraction(replacement_ratio).times_int(GRAPH_REPLACEMENTS)
    constant = constant + graph
    charges = tuple(log2_fraction(value) for value in variables)
    return constant, charges, {
        "type": "linear_bl",
        "variables": [f"{value.numerator}/{value.denominator}" for value in variables],
        "band1": f"{band1.numerator}/{band1.denominator}",
        "graph_log2_interval": [str(graph.lo), str(graph.hi)],
    }


def _exact_graph_linear_outer(
    log_variables: list[float], band1_coefficient: float
) -> tuple[Interval, tuple[Interval, ...], dict[str, Any]]:
    """Outward linear-BL branch with exact graph-hole group moments."""

    if GROUP_BITS != 4 or BLOCK_ATOMS != 16:
        raise ValueError("exact graph linear-BL outer is currently g=4 only")
    variables = tuple(exact_float(math.exp(float(value))) for value in log_variables)
    if len(variables) != 5 or any(value <= 0 for value in variables):
        raise ValueError("exact graph linear-BL variables are invalid")
    band1 = exact_float(float(band1_coefficient))
    band0 = Fraction(1) - band1
    if not Fraction(0) < band0 <= Fraction(1, 2) or not Fraction(1, 2) <= band1 < 1:
        raise ValueError("exact graph linear-BL coefficient is invalid")

    atom = tuple(
        Fraction(math.comb(4, weight)) * variables[weight] for weight in range(5)
    )
    normal_coefficients = _fraction_polynomial_power(atom, 16)
    normal_moments = tuple(
        coefficient / math.comb(64, weight)
        for weight, coefficient in enumerate(normal_coefficients)
    )
    base15 = _fraction_polynomial_power(atom, 15)
    hole_moments = []
    for graph_bit in (0, 1):
        hole_atom = tuple(
            Fraction(math.comb(3, other)) * variables[other + graph_bit]
            for other in range(4)
        )
        coefficients = _fraction_polynomial_multiply(base15, hole_atom)
        hole_moments.append(
            tuple(
                coefficient / math.comb(63, weight)
                for weight, coefficient in enumerate(coefficients)
            )
        )

    def band_norm(coefficient: Fraction, moments, bits: int) -> Interval:
        coefficient_interval = Interval.exact(coefficient.numerator) / Interval.exact(
            coefficient.denominator
        )
        terms = [
            log2_fraction(Fraction(math.comb(bits, weight), 1 << bits)) * LN2
            + (log2_fraction(moment) * LN2) / coefficient_interval
            for weight, moment in enumerate(moments)
        ]
        return coefficient_interval * _ln_sum_exp(terms) / LN2

    normal0 = band_norm(band0, normal_moments, 64)
    normal1 = band_norm(band1, normal_moments, 64)
    hole0 = band_norm(band0, hole_moments[0], 63)
    hole1 = band_norm(band0, hole_moments[1], 63)
    graph = _log2_sum_exp(
        log2_fraction(Fraction(count, 1 << 24))
        + hole0.times_int(128 - weight)
        + hole1.times_int(weight)
        for weight, count in enumerate(load_graph_spectrum())
        if count
    )
    constant = (
        Interval.exact(K)
        + normal0.times_int(256 * 42 - 128)
        + normal1.times_int(256 * 86)
        + graph
    )
    charges = tuple(log2_fraction(value) for value in variables)
    return constant, charges, {
        "type": "exact_graph_linear_bl",
        "variables": [f"{value.numerator}/{value.denominator}" for value in variables],
        "band1": f"{band1.numerator}/{band1.denominator}",
        "normal_band0_log2_interval": [str(normal0.lo), str(normal0.hi)],
        "normal_band1_log2_interval": [str(normal1.lo), str(normal1.hi)],
        "hole0_band0_log2_interval": [str(hole0.lo), str(hole0.hi)],
        "hole1_band0_log2_interval": [str(hole1.lo), str(hole1.hi)],
        "graph_factor_log2_interval": [str(graph.lo), str(graph.hi)],
        "graph_spectrum_sha256": _file_digest(
            Path(__file__).with_name("ebch128_graph24_spectrum.csv")
        ),
        "lemma": "distinct band-zero hole groups; P^15 P_b; exact graph24 average",
    }


@lru_cache(maxsize=1)
def _conditioned_row_split_spectra():
    """Load and validate the exact pair spectra used by the conditioned-row branch."""

    root = Path(__file__).resolve().parent.parent / "out"

    def load(path: Path) -> dict[tuple[int, int], int]:
        table: dict[tuple[int, int], int] = {}
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                key = (int(row["band0_weight"]), int(row["band1_weight"]))
                count = int(row["count"])
                if count <= 0 or key in table:
                    raise ValueError("conditioned-row split spectrum is malformed")
                table[key] = count
        if sum(table.values()) != 1 << 64 or table.get((0, 0)) != 1:
            raise ValueError("conditioned-row split spectrum has wrong mass")
        return table

    path01 = root / "ebch85_band01_split_spectrum.csv"
    path12 = root / "ebch86_band12_split_spectrum.csv"
    spectrum01 = load(path01)
    spectrum12 = load(path12)
    punctured01 = {
        (a, b): (42 - a) * spectrum01.get((a, b), 0)
        + (a + 1) * spectrum01.get((a + 1, b), 0)
        for a in range(42)
        for b in range(44)
    }
    punctured01 = {key: count for key, count in punctured01.items() if count}
    if sum(punctured01.values()) != 42 * (1 << 64):
        raise ValueError("conditioned-row punctured spectrum has wrong mass")
    return spectrum01, punctured01, spectrum12, {
        "spectrum01_sha256": _file_digest(path01),
        "spectrum12_sha256": _file_digest(path12),
    }


def _conditioned_row_exact_graph_outer(
    log_variables: list[float], band1_coefficient: float, theta_value: float
) -> tuple[Interval, tuple[Interval, ...], dict[str, Any]]:
    """Outward one-conditioned-row BL branch with exact graph averaging."""

    if GROUP_BITS != 4 or BLOCK_ATOMS != 16:
        raise ValueError("conditioned-row outer is currently g=4 only")
    variables = tuple(exact_float(math.exp(float(value))) for value in log_variables)
    if len(variables) != 5 or any(value <= 0 for value in variables):
        raise ValueError("conditioned-row variables are invalid")
    band1 = exact_float(float(band1_coefficient))
    band0 = Fraction(1) - band1
    theta = exact_float(float(theta_value))
    if not Fraction(0) < band0 <= Fraction(1, 2):
        raise ValueError("conditioned-row band-zero coefficient is invalid")
    if not Fraction(1, 2) <= band1 < 1:
        raise ValueError("conditioned-row band-one coefficient is invalid")
    if not Fraction(0) <= theta <= 1:
        raise ValueError("conditioned-row Cauchy split is invalid")

    atom = tuple(
        Fraction(math.comb(4, weight)) * variables[weight] for weight in range(5)
    )
    coefficients = _fraction_polynomial_power(atom, 16)
    moments = tuple(
        coefficient / math.comb(64, weight)
        for weight, coefficient in enumerate(coefficients)
    )

    def scale(value: Interval, factor: Fraction) -> Interval:
        return value * (
            Interval.exact(factor.numerator) / Interval.exact(factor.denominator)
        )

    def conditioned_norm(coefficient: Fraction, bit: int) -> Interval:
        coefficient_interval = (
            Interval.exact(coefficient.numerator)
            / Interval.exact(coefficient.denominator)
        )
        terms = [
            log2_fraction(Fraction(math.comb(63, weight), 1 << 63))
            + log2_fraction(moments[weight + bit]) / coefficient_interval
            for weight in range(64)
        ]
        return coefficient_interval * _log2_sum_exp(terms)

    factors = (
        (conditioned_norm(band0, 0), conditioned_norm(band0, 1)),
        (conditioned_norm(band1, 0), conditioned_norm(band1, 1)),
        (conditioned_norm(band1, 0), conditioned_norm(band1, 1)),
    )
    zeros = tuple(pair[0] for pair in factors)
    ratios = tuple(pair[1] - pair[0] for pair in factors)
    spectrum01, punctured01, spectrum12, source_report = (
        _conditioned_row_split_spectra()
    )

    def pair_log_enumerator(table, left: Interval, right: Interval, divisor: int):
        return _log2_sum_exp(
            log2_fraction(Fraction(count, divisor))
            + left.times_int(a)
            + right.times_int(b)
            for (a, b), count in table.items()
        )

    left01 = ratios[0].times_int(2)
    middle01 = scale(ratios[1].times_int(2), theta)
    middle12 = scale(ratios[1].times_int(2), Fraction(1) - theta)
    right12 = ratios[2].times_int(2)

    def cauchy(first, divisor: int) -> Interval:
        return (
            pair_log_enumerator(first, left01, middle01, divisor)
            + pair_log_enumerator(spectrum12, middle12, right12, 1)
        ) / Interval.exact(2)

    normal_pair = cauchy(spectrum01, 1)
    punctured_pair = cauchy(punctured01, 42)
    free_message = Interval.exact(63 * 64)
    normal_tile = (
        free_message
        + zeros[0].times_int(42)
        + zeros[1].times_int(43)
        + zeros[2].times_int(43)
        + normal_pair
    )
    punctured_common = (
        free_message
        + zeros[0].times_int(41)
        + zeros[1].times_int(43)
        + zeros[2].times_int(43)
        + punctured_pair
    )
    hole0_tile = punctured_common + factors[0][0]
    hole1_tile = punctured_common + factors[0][1]
    graph = _log2_sum_exp(
        log2_fraction(Fraction(count, 1 << 24))
        + hole0_tile.times_int(128 - weight)
        + hole1_tile.times_int(weight)
        for weight, count in enumerate(load_graph_spectrum())
        if count
    )
    constant = normal_tile.times_int(128) + graph
    charges = tuple(log2_fraction(value) for value in variables)
    return constant, charges, {
        "type": "conditioned_row_exact_graph",
        "variables": [f"{value.numerator}/{value.denominator}" for value in variables],
        "band1": f"{band1.numerator}/{band1.denominator}",
        "theta": f"{theta.numerator}/{theta.denominator}",
        "normal_tile_log2_interval": [str(normal_tile.lo), str(normal_tile.hi)],
        "hole0_tile_log2_interval": [str(hole0_tile.lo), str(hole0_tile.hi)],
        "hole1_tile_log2_interval": [str(hole1_tile.lo), str(hole1_tile.hi)],
        "graph_factor_log2_interval": [str(graph.lo), str(graph.hi)],
        "graph_spectrum_sha256": _file_digest(
            Path(__file__).with_name("ebch128_graph24_spectrum.csv")
        ),
        **source_report,
        "lemma": (
            "condition one BCH row; apply three-band BL to 63 rows; "
            "pair-spectrum Cauchy; condition punctured row in hole tiles"
        ),
    }


def _exact_graph_total_spectrum_outer(
    log_variables: list[float], log_beta: float
) -> tuple[Interval, tuple[Interval, ...], dict[str, Any]]:
    """Outward total-spectrum branch with exact graph-hole alpha factors."""

    if GROUP_BITS != 4 or BLOCK_ATOMS != 16:
        raise ValueError("exact graph total-spectrum outer is currently g=4 only")
    variables = tuple(exact_float(math.exp(float(value))) for value in log_variables)
    beta = exact_float(math.exp(float(log_beta)))
    if len(variables) != 5 or any(value <= 0 for value in variables) or beta <= 0:
        raise ValueError("exact graph total-spectrum parameters are invalid")
    atom = tuple(
        Fraction(math.comb(4, weight)) * variables[weight] for weight in range(5)
    )
    normal_coefficients = _fraction_polynomial_power(atom, 16)
    normal_alpha = max(
        coefficient / math.comb(64, weight) / beta**weight
        for weight, coefficient in enumerate(normal_coefficients)
    )
    base15 = _fraction_polynomial_power(atom, 15)
    hole_alphas = []
    for graph_bit in (0, 1):
        hole_atom = [Fraction(0)] * 5
        for old_bit in (0, 1):
            for other in range(4):
                hole_atom[old_bit + other] += (
                    Fraction(math.comb(3, other)) * variables[other + graph_bit]
                )
        coefficients = _fraction_polynomial_multiply(base15, hole_atom)
        hole_alphas.append(
            max(
                coefficient / math.comb(64, weight) / beta**weight
                for weight, coefficient in enumerate(coefficients)
            )
        )
    enumerator = sum(
        (
            Fraction(count) * beta**weight
            for weight, count in enumerate(full_spectrum())
            if count
        ),
        Fraction(0),
    )
    graph_factor = sum(
        (
            Fraction(count, 1 << 24)
            * hole_alphas[0] ** (128 - weight)
            * hole_alphas[1] ** weight
            for weight, count in enumerate(load_graph_spectrum())
            if count
        ),
        Fraction(0),
    )
    constant = (
        log2_fraction(normal_alpha).times_int(GROUPS - 128)
        + log2_fraction(enumerator).times_int(K // 64)
        + log2_fraction(graph_factor)
    )
    charges = tuple(log2_fraction(value) for value in variables)
    return constant, charges, {
        "type": "exact_graph_total_spectrum",
        "variables": [f"{value.numerator}/{value.denominator}" for value in variables],
        "beta": f"{beta.numerator}/{beta.denominator}",
        "normal_alpha_log2_interval": [
            str(log2_fraction(normal_alpha).lo),
            str(log2_fraction(normal_alpha).hi),
        ],
        "hole0_alpha_log2_interval": [
            str(log2_fraction(hole_alphas[0]).lo),
            str(log2_fraction(hole_alphas[0]).hi),
        ],
        "hole1_alpha_log2_interval": [
            str(log2_fraction(hole_alphas[1]).lo),
            str(log2_fraction(hole_alphas[1]).hi),
        ],
        "graph_spectrum_sha256": _file_digest(
            Path(__file__).with_name("ebch128_graph24_spectrum.csv")
        ),
        "lemma": "distinct graph holes; punctured-bit alpha; exact graph24 average",
    }


def _total_weight_outer(log_pole: float) -> tuple[Interval, tuple[Interval, ...], dict[str, Any]]:
    """Exact-polynomial/outward-log total-EBCH affine witness."""

    pole = exact_float(math.exp(float(log_pole)))
    if pole <= 0:
        raise ValueError("total-weight pole must be positive")
    spectrum_sum = sum(
        (Fraction(count) * pole**weight for weight, count in enumerate(full_spectrum()) if count),
        Fraction(0),
    )
    replacement_sum = sum(
        (pole ** (-delta) for delta in range(-GRAPH_REPLACEMENTS, GRAPH_REPLACEMENTS + 1)),
        Fraction(0),
    )
    constant = (
        log2_fraction(spectrum_sum).times_int(K // 64)
        + log2_fraction(replacement_sum)
    )
    per_bit = log2_fraction(pole)
    charges = tuple(per_bit.times_int(weight) for weight in range(GROUP_BITS + 1))
    return constant, charges, {
        "type": "total_weight",
        "pole": f"{pole.numerator}/{pole.denominator}",
    }


def _hypergeometric_band0_rows(weight: int) -> tuple[tuple[int, Fraction], ...]:
    selected_min = max(0, 42 - (128 - weight))
    selected_max = min(42, weight)
    return tuple(
        (
            selected,
            Fraction(
                math.comb(weight, selected)
                * math.comb(128 - weight, 42 - selected),
                math.comb(128, 42),
            ),
        )
        for selected in range(selected_min, selected_max + 1)
    )


def _exact_graph_puncture_total_weight_outer(
    pole: Fraction,
) -> tuple[Interval, tuple[Interval, ...], dict[str, Any]]:
    """Exact graph-spectrum/puncture-averaged affine total-weight branch."""

    if not 0 < pole < 1:
        raise ValueError("exact graph-puncture outer pole must lie in (0,1)")
    log_q = log2_fraction(pole)
    # Maclaurin over 128 distinct selected tiles followed by
    # (1+x)^128 <= exp(128x) leaves exp(lambda*J), where J is the
    # hypergeometric band-zero occupancy of one permuted EBCH word.
    lam = (Fraction(1, 1) / pole - 1) / 5376
    lam_interval = Interval.exact(lam.numerator) / Interval.exact(lam.denominator)
    block_terms = []
    for weight, count in enumerate(full_spectrum()):
        if not count:
            continue
        occupancy = _ln_sum_exp(
            log2_fraction(probability) * LN2 + lam_interval.times_int(selected)
            for selected, probability in _hypergeometric_band0_rows(weight)
        )
        block_terms.append(
            log2_fraction(Fraction(count)) * LN2
            + log_q.times_int(weight) * LN2
            + occupancy
        )
    block_log2 = _ln_sum_exp(block_terms) / LN2
    graph_terms = []
    for weight, count in enumerate(load_graph_spectrum()):
        if count:
            graph_terms.append(
                log2_fraction(Fraction(count, 1 << 24)) * LN2
                + log_q.times_int(weight) * LN2
            )
    graph_log2 = _ln_sum_exp(graph_terms) / LN2
    constant = block_log2.times_int(K // 64) + graph_log2
    charges = tuple(log_q.times_int(weight) for weight in range(GROUP_BITS + 1))
    return constant, charges, {
        "type": "exact_graph_puncture_total_weight",
        "pole": f"{pole.numerator}/{pole.denominator}",
        "graph_spectrum_sha256": _file_digest(
            Path(__file__).with_name("ebch128_graph24_spectrum.csv")
        ),
        "block_log2_interval": [str(block_log2.lo), str(block_log2.hi)],
        "graph_factor_log2_interval": [str(graph_log2.lo), str(graph_log2.hi)],
        "puncture_lemma": (
            "distinct-tile Maclaurin; independent coordinate bijections; "
            "hypergeometric(128,w,42); (1+x)^128<=exp(128x)"
        ),
    }


def _load_source_row(row: dict[str, Any], artifact: Path) -> dict[str, Any] | None:
    source_name = row.get("source")
    if source_name is None:
        return None
    # Discovery artifacts may be produced on Windows and verified on Linux.
    # Stored source names are logical relative paths, so accept either native
    # separator without changing the source artifact or its digest.
    source_path = Path(str(source_name).replace("\\", "/"))
    candidates = [source_path]
    if not source_path.is_absolute():
        candidates.extend((artifact.parent / source_path, Path.cwd() / source_path))
    source_path = next((path for path in candidates if path.exists()), None)
    if source_path is None:
        raise ValueError(f"cannot resolve witness source {source_name!r}")
    data = json.loads(source_path.read_text(encoding="utf-8"))
    rows = data if isinstance(data, list) else data.get("rows", [data])
    index = int(row.get("source_index", 0))
    if not 0 <= index < len(rows):
        raise ValueError(f"witness source index {index} is outside {source_path}")
    return rows[index]


def _parameters(row: dict[str, Any], artifact: Path) -> dict[str, Any]:
    source = _load_source_row(row, artifact)
    merged = dict(source or {})
    merged.update(row)

    if "fugacities" in merged:
        fugacities = [float(value) for value in merged["fugacities"]]
    elif "ratio" in merged:
        match = EDGE_NAME.match(str(merged.get("name", "")))
        if match is None:
            raise ValueError("edge witness lacks parseable base/active classes")
        base, active = (int(value) for value in match.groups())
        fugacities = [0.0] * (GROUP_BITS + 1)
        fugacities[base] = 1.0
        fugacities[active] = float(merged["ratio"])
    else:
        raise ValueError("witness lacks fugacities or an edge ratio")
    if len(fugacities) != GROUP_BITS + 1 or any(value < 0 for value in fugacities):
        raise ValueError("invalid witness fugacities")

    outer_type = str(merged.get("outer_type", "linear_bl"))
    outer_details = merged.get("outer_details", {})
    linear_details = merged.get("linear_bl_outer_details", outer_details)
    if outer_type == "linear_bl":
        log_variables = linear_details.get("log_variables", merged.get("outer_log_variables"))
        band1 = linear_details.get("band1_coefficient", 0.5)
        if log_variables is None:
            raise ValueError("linear-BL witness lacks frozen log variables")
        outer = ("linear_bl", [float(value) for value in log_variables], float(band1))
    elif outer_type == "exact_graph_linear_bl":
        log_variables = outer_details.get("log_variables", merged.get("outer_log_variables"))
        band1 = outer_details.get("band1_coefficient", 0.5)
        if log_variables is None:
            raise ValueError("exact graph linear-BL witness lacks frozen log variables")
        outer = (
            "exact_graph_linear_bl",
            [float(value) for value in log_variables],
            float(band1),
        )
    elif outer_type == "conditioned_row_exact_graph":
        log_variables = outer_details.get(
            "log_variables", merged.get("outer_log_variables")
        )
        band1 = outer_details.get("band1_coefficient")
        theta = outer_details.get("pair_cauchy_theta")
        if log_variables is None or band1 is None or theta is None:
            raise ValueError("conditioned-row witness lacks frozen parameters")
        outer = (
            "conditioned_row_exact_graph",
            [float(value) for value in log_variables],
            float(band1),
            float(theta),
        )
    elif outer_type == "exact_graph_total_spectrum":
        log_variables = outer_details.get("log_variables", merged.get("outer_log_variables"))
        outer_point = outer_details.get("outer_point")
        log_beta = outer_details.get("log_beta")
        if log_beta is None and isinstance(outer_point, list) and len(outer_point) > GROUP_BITS:
            log_beta = outer_point[GROUP_BITS]
        if log_variables is None or log_beta is None:
            raise ValueError("exact graph total-spectrum witness lacks frozen parameters")
        outer = (
            "exact_graph_total_spectrum",
            [float(value) for value in log_variables],
            float(log_beta),
        )
    elif outer_type == "total_weight":
        log_pole = outer_details.get("log_pole", merged.get("outer_log_pole"))
        if log_pole is None and source is not None:
            log_pole = source.get("outer_log_pole")
        if log_pole is None:
            raise ValueError("total-weight witness lacks its frozen log pole")
        outer = ("total_weight", float(log_pole))
    elif outer_type == "exact_graph_puncture_total_weight":
        expected_graph_digest = outer_details.get("source_graph_spectrum_sha256")
        actual_graph_digest = _file_digest(
            Path(__file__).with_name("ebch128_graph24_spectrum.csv")
        )
        if (
            expected_graph_digest is not None
            and str(expected_graph_digest) != actual_graph_digest
        ):
            raise ValueError("exact graph-puncture spectrum digest mismatch")
        pole_exact = outer_details.get("pole_exact", merged.get("outer_pole_exact"))
        log_pole = outer_details.get("log_pole", merged.get("outer_log_pole"))
        if pole_exact is not None:
            exact_pole = Fraction(str(pole_exact))
        elif log_pole is not None:
            exact_pole = exact_float(math.exp(float(log_pole)))
        else:
            raise ValueError("exact graph-puncture witness lacks its frozen log pole")
        outer = ("exact_graph_puncture_total_weight", exact_pole)
    else:
        raise ValueError(f"unsupported outer witness type {outer_type!r}")
    return {
        "fugacities": fugacities,
        "pole": float(merged["pole"]),
        "outer": outer,
        "collatz_vector": merged.get(
            "collatz_vector", merged.get("inner_details", {}).get("collatz_vector")
        ),
    }


def harden_witness(
    name: str,
    row: dict[str, Any],
    artifact: Path,
    split_caps,
    iterations: int,
) -> dict[str, Any]:
    # Outer-only atlas rows are already one frozen affine BL branch and do not
    # consume the profile normalization or an inner witness.
    if "fugacities" not in row and "ratio" not in row:
        details = row.get("outer_details", {})
        log_variables = details.get("log_variables", row.get("outer_log_variables"))
        if log_variables is None:
            raise ValueError("outer-only witness lacks frozen log variables")
        outer_constant, outer_charges, outer_report = _linear_outer(
            [float(value) for value in log_variables],
            float(details.get("band1_coefficient", 0.5)),
        )
        return {
            "name": name,
            "constant": outer_constant,
            "charges": outer_charges,
            "subtract_normalization": False,
            "report": {
                "name": name,
                "kind": "outer",
                "outer": outer_report,
                "constant_log2_interval": [str(outer_constant.lo), str(outer_constant.hi)],
                "charge_log2_intervals": [
                    [str(charge.lo), str(charge.hi)] for charge in outer_charges
                ],
            },
        }
    parameters = _parameters(row, artifact)
    fugacities = np.asarray(parameters["fugacities"], dtype=np.float64)
    pole = float(parameters["pole"])
    if not 0 < pole < 1:
        raise ValueError("inner pole must lie in (0,1)")

    diagnostic_histograms = block_histograms(GROUP_BITS, fugacities)
    diagnostic_caps = point_caps(GROUP_BITS, fugacities)
    kernel = SharedDriveStratifiedKernel(
        diagnostic_histograms, diagnostic_caps, split_caps, pole
    )
    frozen_vector = parameters.get("collatz_vector")
    if frozen_vector is None:
        (
            _eigenvalue,
            _domination,
            values,
            _worst,
            selected_iteration,
            _score,
        ) = best_witness(kernel, iterations, INNER_BLOCKS)
        vector_source = "best_finite_block_bound_on_power_trajectory"
    else:
        values = np.asarray(frozen_vector, dtype=np.float64)
        if (
            values.shape != (BLOCK_BITS + 1,)
            or not np.all(np.isfinite(values))
            or np.any(values <= 0.0)
        ):
            raise ValueError(
                "frozen Collatz vector must contain 65 positive finite values"
            )
        selected_iteration = None
        vector_source = "artifact_frozen_exact_dyadics"
    histograms_upper = _block_histograms_outward(fugacities)
    caps_upper = _point_caps_outward(fugacities)
    if not np.all(histograms_upper >= diagnostic_histograms):
        raise ValueError("outward histogram does not enclose diagnostic histogram")
    if not np.all(caps_upper >= diagnostic_caps):
        raise ValueError("outward cap does not enclose diagnostic cap")
    image, transport_stats = outward_transport_image(
        values, histograms_upper, caps_upper, split_caps, pole
    )
    value_fractions = [exact_float(float(value)) for value in values]
    eigenvalue_upper = max(entry / value for entry, value in zip(image, value_fractions))
    domination = max(Fraction(1) / value for value in value_fractions)
    inner_mgf = (
        log2_fraction(domination)
        + log2_fraction(eigenvalue_upper).times_int(INNER_BLOCKS)
        + log2_fraction(value_fractions[0])
    )
    inner_constant = inner_mgf - log2_fraction(exact_float(pole)).times_int(D)
    inner_charges = tuple(
        None if value == 0 else log2_fraction(exact_float(float(value)))
        for value in fugacities
    )

    outer = parameters["outer"]
    if outer[0] == "linear_bl":
        outer_constant, outer_charges, outer_report = _linear_outer(outer[1], outer[2])
    elif outer[0] == "exact_graph_linear_bl":
        outer_constant, outer_charges, outer_report = _exact_graph_linear_outer(
            outer[1], outer[2]
        )
    elif outer[0] == "conditioned_row_exact_graph":
        outer_constant, outer_charges, outer_report = (
            _conditioned_row_exact_graph_outer(outer[1], outer[2], outer[3])
        )
    elif outer[0] == "exact_graph_total_spectrum":
        outer_constant, outer_charges, outer_report = _exact_graph_total_spectrum_outer(
            outer[1], outer[2]
        )
    elif outer[0] == "total_weight":
        outer_constant, outer_charges, outer_report = _total_weight_outer(outer[1])
    else:
        outer_constant, outer_charges, outer_report = (
            _exact_graph_puncture_total_weight_outer(outer[1])
        )
    charges = tuple(
        None if inner is None else outer_charge + inner
        for outer_charge, inner in zip(outer_charges, inner_charges)
    )
    constant = outer_constant + inner_constant
    return {
        "name": name,
        "constant": constant,
        "charges": charges,
        "subtract_normalization": True,
        "report": {
            "name": name,
            "outer": outer_report,
            "fugacities": fugacities.tolist(),
            "pole": pole,
            "collatz_vector_source": vector_source,
            "collatz_selected_iteration": selected_iteration,
            "collatz_vector": [float(value) for value in values],
            "outward_float_invariant": (
                "positive dyadic multiplication must remain positive in "
                "binary64; any underflow aborts certification"
            ),
            "inner_constant_log2_interval": [
                str(inner_constant.lo),
                str(inner_constant.hi),
            ],
            "inner_charge_log2_intervals": [
                None if charge is None else [str(charge.lo), str(charge.hi)]
                for charge in inner_charges
            ],
            "constant_log2_interval": [str(constant.lo), str(constant.hi)],
            "charge_log2_intervals": [
                None if charge is None else [str(charge.lo), str(charge.hi)]
                for charge in charges
            ],
            **transport_stats,
        },
    }


def _serialize_interval(value: Interval) -> list[str]:
    return [str(value.lo), str(value.hi)]


def _deserialize_interval(value: list[str]) -> Interval:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError("invalid cached interval")
    return Interval(Decimal(value[0]), Decimal(value[1]))


def _serialize_hardened(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": value["name"],
        "constant": _serialize_interval(value["constant"]),
        "charges": [
            None if charge is None else _serialize_interval(charge)
            for charge in value["charges"]
        ],
        "subtract_normalization": bool(value["subtract_normalization"]),
        "report": value["report"],
    }


def _deserialize_hardened(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": str(value["name"]),
        "constant": _deserialize_interval(value["constant"]),
        "charges": tuple(
            None if charge is None else _deserialize_interval(charge)
            for charge in value["charges"]
        ),
        "subtract_normalization": bool(value["subtract_normalization"]),
        "report": value["report"],
    }


@lru_cache(maxsize=None)
def _file_digest(path: Path) -> str:
    """Hash an immutable certificate input once per verifier process."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def _resolved_source_digest(row: dict[str, Any], artifact: Path) -> str | None:
    source_name = row.get("source")
    if source_name is None:
        return None
    source_path = Path(str(source_name).replace("\\", "/"))
    candidates = [source_path]
    if not source_path.is_absolute():
        candidates.extend((artifact.parent / source_path, Path.cwd() / source_path))
    resolved = next((path.resolve() for path in candidates if path.exists()), None)
    return None if resolved is None else _file_digest(resolved)


def _cache_key(
    name: str,
    row: dict[str, Any],
    artifact: Path,
    iterations: int,
) -> str:
    script_root = Path(__file__).resolve().parent
    dependency_names = (
        "certify_packet_group_triangle_ledger.py",
        "EBCH128_64.wd",
        "ebch128_systematic_split_slices.csv",
        "ebch128_64_spectrum.csv",
        "ebch128_graph24_spectrum.csv",
    )
    dependencies = {
        item: (
            os.environ.get("G2_HARDENING_CODE_SHA256", _file_digest(script_root / item))
            if item == "certify_packet_group_triangle_ledger.py"
            else _file_digest(script_root / item)
        )
        for item in dependency_names
        if (script_root / item).exists()
    }
    payload = {
        "schema": "g2-outward-hardened-witness-cache-v1",
        "name": name,
        "row": row,
        "artifact_sha256": _file_digest(artifact.resolve()),
        "nested_source_sha256": _resolved_source_digest(row, artifact),
        "iterations": iterations,
        "dependencies": dependencies,
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _cache_path(directory: Path, name: str, key: str) -> Path:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)[:80]
    return directory / f"{safe_name}.{key}.json"


def _read_cached_hardened(path: Path, expected_key: str) -> dict[str, Any]:
    cached = json.loads(path.read_text(encoding="utf-8"))
    if cached.get("schema") != "g2-outward-hardened-witness-cache-v1":
        raise ValueError(f"invalid hardened witness cache schema: {path}")
    if cached.get("cache_key") != expected_key:
        raise ValueError(f"hardened witness cache key mismatch: {path}")
    return _deserialize_hardened(cached["hardened"])


def _write_cached_hardened(
    path: Path,
    cache_key: str,
    hardened: dict[str, Any],
) -> None:
    payload = {
        "schema": "g2-outward-hardened-witness-cache-v1",
        "cache_key": cache_key,
        "hardened": _serialize_hardened(hardened),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _initialize_hardening_worker() -> None:
    global _WORKER_SPLIT_CAPS
    _WORKER_SPLIT_CAPS = split_cap_table()


def _harden_worker(task) -> tuple[str, dict[str, Any]]:
    name, row, artifact_text, iterations = task
    if _WORKER_SPLIT_CAPS is None:
        raise RuntimeError("outward hardening worker was not initialized")
    hardened = harden_witness(
        name,
        row,
        Path(artifact_text),
        _WORKER_SPLIT_CAPS,
        iterations,
    )
    return name, hardened


def harden_selected_witnesses(
    used_names: list[str],
    witnesses: dict[str, tuple[dict[str, Any], Path]],
    iterations: int,
    workers: int,
    checkpoint_dir: Path | None,
) -> dict[str, dict[str, Any]]:
    hardened: dict[str, dict[str, Any]] = {}
    tasks = []
    cache_records = {}
    for name in used_names:
        if name == "full_bijection":
            hardened[name] = harden_full_bijection()
            print("hardened witness=full_bijection source=builtin", flush=True)
            continue
        row, artifact = witnesses[name]
        key = _cache_key(name, row, artifact, iterations)
        cache_path = None if checkpoint_dir is None else _cache_path(checkpoint_dir, name, key)
        if cache_path is not None and cache_path.exists():
            hardened[name] = _read_cached_hardened(cache_path, key)
            print(f"hardened witness={name} source=checkpoint", flush=True)
            continue
        tasks.append((name, row, str(artifact), iterations))
        cache_records[name] = (key, cache_path)

    if not tasks:
        return {name: hardened[name] for name in used_names}

    completed = len(hardened)
    if workers == 1:
        split_caps = split_cap_table()
        for name, row, artifact_text, local_iterations in tasks:
            value = harden_witness(
                name,
                row,
                Path(artifact_text),
                split_caps,
                local_iterations,
            )
            hardened[name] = value
            key, cache_path = cache_records[name]
            if cache_path is not None:
                _write_cached_hardened(cache_path, key, value)
            completed += 1
            print(
                f"hardened witness={name} progress={completed}/{len(used_names)}",
                flush=True,
            )
    else:
        with ProcessPoolExecutor(
            max_workers=min(workers, len(tasks)),
            initializer=_initialize_hardening_worker,
        ) as executor:
            future_to_name = {
                executor.submit(_harden_worker, task): task[0] for task in tasks
            }
            for future in as_completed(future_to_name):
                name = future_to_name[future]
                try:
                    returned_name, value = future.result()
                except BaseException as error:
                    raise RuntimeError(f"outward hardening failed for witness {name}") from error
                if returned_name != name:
                    raise RuntimeError("outward hardening worker returned the wrong witness")
                hardened[name] = value
                key, cache_path = cache_records[name]
                if cache_path is not None:
                    _write_cached_hardened(cache_path, key, value)
                completed += 1
                print(
                    f"hardened witness={name} progress={completed}/{len(used_names)}",
                    flush=True,
                )
    return {name: hardened[name] for name in used_names}


def _rows_from_artifact(path: Path) -> list[dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = (
        data
        if isinstance(data, list)
        else data.get(
            "rows",
            data.get(
                "witnesses",
                [data] if isinstance(data, dict) and "fugacities" in data else [],
            ),
        )
    )
    if isinstance(rows, dict):
        rows = list(rows.values())
    if not isinstance(rows, list):
        raise ValueError(f"witness artifact {path} has no row array")
    return rows


def load_witnesses(paths: list[Path]) -> dict[str, tuple[dict[str, Any], Path]]:
    records: list[tuple[str, dict[str, Any], Path, int]] = []
    for path in paths:
        for index, row in enumerate(_rows_from_artifact(path)):
            name = str(row.get("name", f"{path.name}:{index}"))
            records.append((name, row, path, index))
    name_counts = Counter(name for name, _row, _path, _index in records)
    result: dict[str, tuple[dict[str, Any], Path]] = {}
    for name, row, path, index in records:
        # Canonical indexed identifiers remain available even when discovery
        # artifacts contain several refinements with the same display name.
        aliases = (f"{path.name}:{index}", f"{path}:{index}")
        for alias in aliases:
            if alias in result:
                raise ValueError(f"duplicate canonical witness identifier {alias!r}")
            result[alias] = (row, path)
        if name_counts[name] == 1:
            result[name] = (row, path)
    return result


def _witness_name(reference: Any, selected: list[str]) -> str:
    if isinstance(reference, dict):
        reference = reference.get("witness", reference.get("witness_name", reference.get("name")))
    if isinstance(reference, int):
        if not 0 <= reference < len(selected):
            raise ValueError(f"witness index {reference} is out of range")
        return selected[reference]
    if reference is None:
        raise ValueError("vertex assignment lacks a witness")
    return str(reference)


def _parse_weight(value: Any) -> Fraction:
    """Parse a serialized rational, preserving a JSON float's exact dyadic value."""

    if isinstance(value, bool):
        raise ValueError("mixture weights must be numeric, not booleans")
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("mixture weights must be finite")
        return Fraction.from_float(value)
    if isinstance(value, str):
        item = value.strip()
        try:
            if "/" in item:
                return Fraction(item)
            if item.lower().startswith(("0x", "+0x", "-0x")):
                parsed = float.fromhex(item)
                if not math.isfinite(parsed):
                    raise ValueError
                return Fraction.from_float(parsed)
            return Fraction(Decimal(item))
        except (ArithmeticError, ValueError) as error:
            raise ValueError(f"invalid exact mixture weight {value!r}") from error
    if isinstance(value, dict) and set(value) >= {"numerator", "denominator"}:
        try:
            return Fraction(int(value["numerator"]), int(value["denominator"]))
        except (ArithmeticError, TypeError, ValueError) as error:
            raise ValueError(f"invalid rational mixture weight {value!r}") from error
    raise ValueError(f"unsupported mixture weight encoding {value!r}")


def _weight_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _mixture_from_row(
    row: dict[str, Any], selected: list[str]
) -> tuple[tuple[str, Fraction], ...] | None:
    """Return a canonical fixed mixture, or ``None`` for the legacy path."""

    forbidden = (
        "vertex_mixture_weights",
        "vertex_mixture_witness_references",
        "vertex_mixtures",
    )
    if any(key in row for key in forbidden):
        raise ValueError("vertex-specific mixture weights are not certifiable")
    references = row.get("mixture_witness_references", row.get("mixture_witnesses"))
    exact_weights = row.get("mixture_weights_exact")
    weights = exact_weights if exact_weights is not None else row.get("mixture_weights")
    if references is None and weights is None:
        return None
    if references is None or weights is None:
        raise ValueError("fixed mixture requires both witness references and weights")
    if not isinstance(references, list) or not isinstance(weights, list):
        raise ValueError("fixed mixture witness references and weights must be flat arrays")
    if len(references) != len(weights) or not references:
        raise ValueError("fixed mixture arrays must have the same nonzero length")

    combined: dict[str, Fraction] = {}
    for reference, raw_weight in zip(references, weights):
        if isinstance(reference, (list, tuple)):
            raise ValueError("mixture witness references must not vary by vertex")
        weight = _parse_weight(raw_weight)
        if weight < 0:
            raise ValueError("mixture weights must be nonnegative")
        name = _witness_name(reference, selected)
        combined[name] = combined.get(name, Fraction(0)) + weight
    total_weight = sum(combined.values(), Fraction(0))
    if total_weight != 1:
        raise ValueError(
            f"fixed mixture weights must sum exactly to 1 (got {_weight_text(total_weight)})"
        )
    # Zero components are mathematically absent.  Positive components are all
    # evaluated at every vertex, including components that cannot support it.
    return tuple(sorted((name, weight) for name, weight in combined.items() if weight))


def extract_assignments(
    ledger: dict[str, Any],
) -> list[tuple[tuple[int, int, int], tuple[tuple[str, Fraction], ...]]]:
    """Flatten supported ledger schemas to unique profile/fixed-mixture pairs."""

    selected = [
        str(item.get("name", item.get("witness", item))) if isinstance(item, dict) else str(item)
        for item in ledger.get("selected_witnesses", [])
    ]
    pairs: list[tuple[tuple[int, int, int], tuple[tuple[str, Fraction], ...]]] = []

    def add(
        profile: Any,
        reference: Any = None,
        mixture: tuple[tuple[str, Fraction], ...] | None = None,
    ) -> None:
        if isinstance(profile, dict):
            if reference is None and mixture is None:
                reference = profile.get("witness", profile.get("witness_name"))
            profile = profile.get("profile", profile.get("vertex"))
        if not isinstance(profile, (list, tuple)) or len(profile) != 3:
            raise ValueError("g=2 vertex profile must contain three counts")
        vertex = tuple(int(value) for value in profile)
        if any(value < 0 for value in vertex) or sum(vertex) != atom_count(GROUP_BITS):
            raise ValueError(f"invalid g=2 profile vertex {vertex}")
        if mixture is None:
            mixture = ((_witness_name(reference, selected), Fraction(1)),)
        pairs.append((vertex, mixture))

    flat = next(
        (ledger[key] for key in ("used_inequalities", "inequalities", "vertex_assignments") if key in ledger),
        None,
    )
    if flat is not None:
        for row in flat:
            mixture = _mixture_from_row(row, selected)
            reference = row.get(
                "witness", row.get("witness_reference", row.get("witness_name"))
            )
            add(row.get("profile", row.get("vertex")), reference, mixture)
    else:
        triangles = ledger.get("triangles") or ledger.get("triangle_to_witness_ledger") or ledger.get("cells", [])
        for triangle in triangles:
            vertices = triangle.get(
                "vertices", triangle.get("vertex_profiles", triangle.get("profiles"))
            )
            if not isinstance(vertices, list) or len(vertices) != 3:
                raise ValueError("each triangular cell must contain three vertices")
            mixture = _mixture_from_row(triangle, selected)
            if mixture is not None:
                for vertex in vertices:
                    add(vertex, mixture=mixture)
                continue
            references = triangle.get("vertex_witnesses", triangle.get("witnesses"))
            if references is None:
                reference = triangle.get(
                    "witness",
                    triangle.get(
                        "witness_reference",
                        triangle.get("witness_name", triangle.get("witness_index")),
                    ),
                )
                references = [reference] * 3
            if len(references) != 3:
                raise ValueError("triangle vertex-witness array must have length three")
            resolved = [_witness_name(reference, selected) for reference in references]
            if len(set(resolved)) != 1:
                raise ValueError(
                    "one fixed witness must be used at all three vertices of a triangle"
                )
            for vertex, reference in zip(vertices, references):
                add(vertex, reference)
        for segment in ledger.get("segments", ledger.get("edge_segments", [])):
            vertices = segment.get("vertices", segment.get("profiles"))
            if not isinstance(vertices, list) or len(vertices) != 2:
                raise ValueError("each edge segment must contain two vertices")
            mixture = _mixture_from_row(segment, selected)
            reference = segment.get(
                "witness",
                segment.get(
                    "witness_reference", segment.get("witness_name", segment.get("witness_index"))
                ),
            )
            for vertex in vertices:
                add(vertex, reference, mixture)
    if not pairs:
        raise ValueError("triangle ledger contains no witness-at-vertex assignments")
    return sorted(set(pairs))


def extract_discrete_assignments(
    ledger: dict[str, Any],
) -> list[tuple[tuple[int, int, int], tuple[tuple[str, Fraction], ...]]]:
    """Parse the explicitly assigned lattice profiles in hybrid terminal cells."""

    rows = ledger.get("discrete_point_assignments", [])
    if not isinstance(rows, list):
        raise ValueError("discrete point assignments must be an array")
    if not rows:
        return []
    return extract_assignments(
        {
            "selected_witnesses": ledger.get("selected_witnesses", []),
            "used_inequalities": rows,
        }
    )


def _profile_vertex(value: Any) -> tuple[int, int, int]:
    if isinstance(value, dict):
        value = value.get("profile", value.get("vertex"))
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError("geometry vertex must be a g=2 profile")
    result = tuple(int(item) for item in value)
    if any(item < 0 for item in result) or sum(result) != atom_count(GROUP_BITS):
        raise ValueError(f"invalid geometry profile {result}")
    return result


def _point(profile: tuple[int, int, int]) -> tuple[int, int]:
    return profile[1], profile[2]


def _cross(left: tuple[int, int], middle: tuple[int, int], right: tuple[int, int]) -> int:
    return (middle[0] - left[0]) * (right[1] - left[1]) - (
        middle[1] - left[1]
    ) * (right[0] - left[0])


def _closed_triangle_lattice_count(
    vertices: tuple[tuple[int, int, int], ...],
) -> tuple[int, int, int]:
    """Return exact ``(area2, boundary, total)`` for a closed lattice triangle."""

    if len(vertices) != 3:
        raise ValueError("a triangular cell must contain exactly three vertices")
    points = tuple(_point(vertex) for vertex in vertices)
    area2 = abs(_cross(points[0], points[1], points[2]))
    if area2 == 0:
        raise ValueError("triangular aggregation cell is degenerate")
    boundary = sum(
        math.gcd(
            abs(points[(index + 1) % 3][0] - points[index][0]),
            abs(points[(index + 1) % 3][1] - points[index][1]),
        )
        for index in range(3)
    )
    numerator = area2 + boundary + 2
    if numerator % 2:
        raise ValueError("Pick lattice count is not integral")
    return area2, boundary, numerator // 2


def extract_triangle_cells(
    ledger: dict[str, Any],
) -> list[
    tuple[
        int,
        tuple[tuple[int, int, int], ...],
        tuple[tuple[str, Fraction], ...],
        int,
    ]
]:
    """Extract fixed-witness cells with independently checked lattice counts."""

    selected = [
        str(item.get("name", item.get("witness", item))) if isinstance(item, dict) else str(item)
        for item in ledger.get("selected_witnesses", [])
    ]
    rows = ledger.get("triangles") or ledger.get("triangle_to_witness_ledger") or ledger.get("cells", [])
    result = []
    for index, row in enumerate(rows):
        raw_vertices = row.get("vertices", row.get("vertex_profiles", row.get("profiles")))
        if not isinstance(raw_vertices, list) or len(raw_vertices) != 3:
            raise ValueError("each aggregation cell must contain three profile vertices")
        vertices = tuple(_profile_vertex(vertex) for vertex in raw_vertices)
        mixture = _mixture_from_row(row, selected)
        if mixture is None:
            references = row.get("vertex_witnesses", row.get("witnesses"))
            if references is None:
                reference = row.get(
                    "witness",
                    row.get(
                        "witness_reference",
                        row.get("witness_name", row.get("witness_index")),
                    ),
                )
                references = [reference] * 3
            if not isinstance(references, list) or len(references) != 3:
                raise ValueError("triangle vertex-witness array must have length three")
            resolved = tuple(_witness_name(reference, selected) for reference in references)
            if len(set(resolved)) != 1:
                raise ValueError("cell-local aggregation requires one fixed witness per triangle")
            mixture = ((resolved[0], Fraction(1)),)

        area2, boundary, total = _closed_triangle_lattice_count(vertices)
        for key, actual in (
            ("area2", area2),
            ("boundary_lattice_points", boundary),
            ("total_lattice_points", total),
        ):
            reported = row.get(key)
            if reported is not None and int(reported) != actual:
                raise ValueError(
                    f"triangle {index} {key} mismatch: reported {reported}, exact {actual}"
                )
        cell_id = int(row.get("triangle_id", index))
        result.append((cell_id, vertices, mixture, total))
    return result


def _polygon_area2(polygon: tuple[tuple[int, int], ...]) -> int:
    return sum(
        left[0] * right[1] - left[1] * right[0]
        for left, right in zip(polygon, polygon[1:] + polygon[:1])
    )


def _ccw_polygon(points: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    polygon = tuple(points)
    area = _polygon_area2(polygon)
    if not area:
        raise ValueError("geometry polygon is degenerate")
    return polygon if area > 0 else tuple(reversed(polygon))


def _on_segment(point, left, right) -> bool:
    return _cross(left, right, point) == 0 and all(
        min(a, b) <= value <= max(a, b)
        for value, a, b in zip(point, left, right)
    )


def _inside_convex(point, polygon) -> bool:
    return all(
        _cross(left, right, point) >= 0
        for left, right in zip(polygon, polygon[1:] + polygon[:1])
    )


def _mesh_covers_polygon(
    triangle_rows: list[dict[str, Any]], polygon_profiles: tuple[tuple[int, int, int], ...]
) -> dict[str, int]:
    """Exact oriented-edge and area audit for a conforming triangular mesh."""

    polygon = _ccw_polygon(_point(profile) for profile in polygon_profiles)
    raw_edges: list[tuple[tuple[int, int], tuple[int, int]]] = []
    area_sum = 0
    for row in triangle_rows:
        raw_vertices = row.get("vertices", row.get("vertex_profiles", row.get("profiles")))
        profiles = tuple(_profile_vertex(value) for value in raw_vertices)
        points = tuple(_point(profile) for profile in profiles)
        area = _cross(points[0], points[1], points[2])
        if not area:
            raise ValueError("triangular cover contains a degenerate triangle")
        if area < 0:
            points = (points[0], points[2], points[1])
            area = -area
        if not all(_inside_convex(point, polygon) for point in points):
            raise ValueError("triangle vertex lies outside the certified hull")
        area_sum += area
        for left, right in zip(points, points[1:] + points[:1]):
            raw_edges.append((left, right))

    # Adaptive neighboring triangles may create a conforming T-junction: one
    # side retains a long edge while the other side contains two subedges.
    # Group edges by their exact primitive integer line, sort that line's
    # existing endpoints once, and split each edge by binary search.  This is
    # O(E log E + P), where P is the number of resulting edge pieces; the old
    # all-mesh-vertices-per-edge scan was quadratic and dominated large ledgers.
    mesh_points = {point for edge in raw_edges for point in edge}
    line_groups: dict[
        tuple[int, int, int],
        list[tuple[tuple[int, int], tuple[int, int], int]],
    ] = defaultdict(list)
    line_points: dict[
        tuple[int, int, int], dict[int, tuple[int, int]]
    ] = defaultdict(dict)
    for left, right in raw_edges:
        dx = right[0] - left[0]
        dy = right[1] - left[1]
        divisor = math.gcd(abs(dx), abs(dy))
        if not divisor:
            raise ValueError("triangular mesh contains a zero-length edge")
        ux = dx // divisor
        uy = dy // divisor
        if ux < 0 or (ux == 0 and uy < 0):
            ux = -ux
            uy = -uy
        # (uy,-ux) is a primitive normal.  Its dot product is constant on
        # precisely this infinite integer line.
        line = (ux, uy, uy * left[0] - ux * left[1])
        axis = 0 if ux else 1
        line_groups[line].append((left, right, axis))
        line_points[line][left[axis]] = left
        line_points[line][right[axis]] = right

    edges: Counter[tuple[tuple[int, int], tuple[int, int]]] = Counter()
    directed_edges: Counter[tuple[tuple[int, int], tuple[int, int]]] = Counter()
    for line, rows in line_groups.items():
        coordinate_to_point = line_points[line]
        coordinates = sorted(coordinate_to_point)
        for left, right, axis in rows:
            left_coordinate = left[axis]
            right_coordinate = right[axis]
            low = min(left_coordinate, right_coordinate)
            high = max(left_coordinate, right_coordinate)
            begin_index = bisect.bisect_left(coordinates, low)
            end_index = bisect.bisect_right(coordinates, high)
            local_coordinates = coordinates[begin_index:end_index]
            if left_coordinate > right_coordinate:
                local_coordinates = list(reversed(local_coordinates))
            points = [coordinate_to_point[value] for value in local_coordinates]
            if not points or points[0] != left or points[-1] != right:
                raise ValueError("exact line sweep lost a mesh edge endpoint")
            for begin, end in zip(points, points[1:]):
                directed_edges[(begin, end)] += 1
                edge = (begin, end) if begin < end else (end, begin)
                edges[edge] += 1

    boundary_intervals: list[list[tuple[Fraction, Fraction]]] = [
        [] for _ in polygon
    ]
    interior_edges = 0
    for (left, right), multiplicity in edges.items():
        if multiplicity == 2:
            if directed_edges[(left, right)] != 1 or directed_edges[(right, left)] != 1:
                raise ValueError("paired mesh edge does not have opposite orientations")
            interior_edges += 1
            continue
        if multiplicity != 1:
            raise ValueError("triangular mesh edge has multiplicity greater than two")
        matches = [
            index
            for index, (begin, end) in enumerate(
                zip(polygon, polygon[1:] + polygon[:1])
            )
            if _on_segment(left, begin, end) and _on_segment(right, begin, end)
        ]
        if len(matches) != 1:
            raise ValueError("unpaired triangle edge is not on exactly one hull edge")
        index = matches[0]
        begin, end = polygon[index], polygon[(index + 1) % len(polygon)]
        axis = 0 if begin[0] != end[0] else 1
        denominator = end[axis] - begin[axis]
        first = Fraction(left[axis] - begin[axis], denominator)
        second = Fraction(right[axis] - begin[axis], denominator)
        boundary_intervals[index].append((min(first, second), max(first, second)))

    for intervals in boundary_intervals:
        cursor = Fraction(0)
        for low, high in sorted(intervals):
            if low != cursor or high <= low:
                raise ValueError("triangular mesh does not cover a hull edge exactly")
            cursor = high
        if cursor != 1:
            raise ValueError("triangular mesh leaves a hull-edge gap")
    polygon_area = _polygon_area2(polygon)
    if area_sum != polygon_area:
        raise ValueError(
            f"triangle doubled area {area_sum} does not equal hull area {polygon_area}"
        )
    return {
        "triangles": len(triangle_rows),
        "mesh_vertices": len(mesh_points),
        "interior_edges": interior_edges,
        "doubled_area": area_sum,
    }


def _edge_segments_cover(segments: list[dict[str, Any]]) -> dict[str, int]:
    """Check exact coverage of the three support-two integer edge segments."""

    expected = {
        2: (21, atom_count(GROUP_BITS)),  # a2=0, parameter a1
        1: (11, atom_count(GROUP_BITS)),  # a1=0, parameter a2
        0: (0, atom_count(GROUP_BITS)),   # a0=0, parameter a2
    }
    intervals: dict[int, list[tuple[int, int]]] = {index: [] for index in expected}
    for row in segments:
        raw_vertices = row.get("vertices", row.get("vertex_profiles", row.get("profiles")))
        vertices = tuple(_profile_vertex(value) for value in raw_vertices)
        zero = [index for index in range(3) if vertices[0][index] == vertices[1][index] == 0]
        if len(zero) != 1:
            raise ValueError("edge segment is not contained in one support-two face")
        index = zero[0]
        parameter = 1 if index == 2 else 2
        values = sorted((vertices[0][parameter], vertices[1][parameter]))
        if values[0] == values[1]:
            raise ValueError("edge segment is degenerate")
        intervals[index].append((values[0], values[1]))
    for index, (begin, end) in expected.items():
        cursor = begin
        for low, high in sorted(intervals[index]):
            if low > cursor or high <= cursor:
                raise ValueError(f"support face a{index}=0 is not exactly covered")
            cursor = max(cursor, high)
        if cursor != end:
            raise ValueError(f"support face a{index}=0 has an uncovered tail")
    return {"edge_segments": len(segments)}


def _terminal_triangles(ledger: dict[str, Any]) -> list[dict[str, Any]]:
    rows = ledger.get("terminal_triangles", ledger.get("terminal_triangle_cells", []))
    if not isinstance(rows, list):
        raise ValueError("terminal triangle cells must be an array")
    return rows


def _ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def _triangle_lattice_profiles(row: dict[str, Any]) -> set[tuple[int, int, int]]:
    """Enumerate a triangle's lattice profiles in scanline-linear time."""

    raw = row.get("vertices", row.get("vertex_profiles", row.get("profiles")))
    if not isinstance(raw, list) or len(raw) != 3:
        raise ValueError("each terminal triangle must contain three vertices")
    profiles = tuple(_profile_vertex(value) for value in raw)
    points = tuple(_point(profile) for profile in profiles)
    area = _cross(points[0], points[1], points[2])
    if not area:
        raise ValueError("terminal triangle is degenerate")
    if area < 0:
        points = (points[0], points[2], points[1])

    # Scan along the narrower coordinate.  Each exact edge half-plane gives
    # a rational lower or upper bound for the other coordinate, avoiding a
    # potentially quadratic bounding-box walk for long thin terminal cells.
    swap = max(point[0] for point in points) - min(point[0] for point in points) > (
        max(point[1] for point in points) - min(point[1] for point in points)
    )
    scan_points = tuple((y, x) if swap else (x, y) for x, y in points)
    if _cross(scan_points[0], scan_points[1], scan_points[2]) < 0:
        scan_points = (scan_points[0], scan_points[2], scan_points[1])
    low_scan = min(point[0] for point in scan_points)
    high_scan = max(point[0] for point in scan_points)
    low_other = min(point[1] for point in scan_points)
    high_other = max(point[1] for point in scan_points)
    result: set[tuple[int, int, int]] = set()
    for scan in range(low_scan, high_scan + 1):
        lower = Fraction(low_other)
        upper = Fraction(high_other)
        feasible = True
        for left, right in zip(scan_points, scan_points[1:] + scan_points[:1]):
            dx = right[0] - left[0]
            dy = right[1] - left[1]
            if dx == 0:
                if -dy * (scan - left[0]) < 0:
                    feasible = False
                    break
                continue
            boundary = Fraction(dy * (scan - left[0]), dx) + left[1]
            if dx > 0:
                lower = max(lower, boundary)
            else:
                upper = min(upper, boundary)
        if not feasible:
            continue
        for other in range(_ceil_fraction(lower), upper.numerator // upper.denominator + 1):
            a1, a2 = (other, scan) if swap else (scan, other)
            profile = (atom_count(GROUP_BITS) - a1 - a2, a1, a2)
            if min(profile) < 0:
                raise ValueError("terminal triangle produced a profile outside the simplex")
            result.add(profile)
    return result


def _profiles_sha256(profiles: Iterable[tuple[int, int, int]]) -> str:
    payload = [list(profile) for profile in sorted(profiles)]
    encoded = json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def _canonical_triangle_vertices(row: dict[str, Any]) -> tuple[tuple[int, int, int], ...]:
    raw = row.get("vertices", row.get("vertex_profiles", row.get("profiles")))
    profiles = tuple(_profile_vertex(value) for value in raw)
    points = tuple(_point(profile) for profile in profiles)
    if _cross(points[0], points[1], points[2]) < 0:
        profiles = (profiles[0], profiles[2], profiles[1])
    # Match the producer's canonical triangle: CCW in (a1,a2), rotated to
    # begin at the lexicographically least (a1,a2) point.
    start = min(range(3), key=lambda index: _point(profiles[index]))
    return profiles[start:] + profiles[:start]


def verify_discrete_terminal_cover(
    ledger: dict[str, Any],
    assignments: list[tuple[tuple[int, int, int], tuple[tuple[str, Fraction], ...]]],
) -> dict[str, Any]:
    """Exhaustively match terminal-cell lattice points to fixed assignments."""

    terminals = _terminal_triangles(ledger)
    expected: set[tuple[int, int, int]] = set()
    occurrences = 0
    cell_counts = []
    cell_hashes = []
    tagged_cells = []
    for index, row in enumerate(terminals):
        profiles = _triangle_lattice_profiles(row)
        count = len(profiles)
        profile_hash = _profiles_sha256(profiles)
        reported_count = row.get("terminal_integer_lattice_count")
        if reported_count is not None and int(reported_count) != count:
            raise ValueError(
                f"terminal triangle {index} lattice count mismatch: {reported_count} != {count}"
            )
        reported_hash = row.get("terminal_integer_profiles_sha256")
        if reported_hash is not None and str(reported_hash) != profile_hash:
            raise ValueError(
                f"terminal triangle {index} profile digest mismatch: "
                f"{reported_hash} != {profile_hash}"
            )
        cell_counts.append(count)
        cell_hashes.append(profile_hash)
        tagged_cells.append(
            {
                "vertices": [list(profile) for profile in _canonical_triangle_vertices(row)],
                "profiles": [list(profile) for profile in sorted(profiles)],
            }
        )
        occurrences += count
        expected.update(profiles)

    tagged_cells.sort(key=lambda item: item["vertices"])
    cell_lattice_hash = hashlib.sha256(
        json.dumps(tagged_cells, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()
    unique_hash = _profiles_sha256(expected)
    reported_cell_hash = ledger.get("terminal_cell_lattice_sha256")
    if reported_cell_hash is not None and str(reported_cell_hash) != cell_lattice_hash:
        raise ValueError(
            "global terminal-cell lattice digest mismatch: "
            f"{reported_cell_hash} != {cell_lattice_hash}"
        )
    reported_unique_hash = ledger.get("terminal_unique_integer_profiles_sha256")
    if reported_unique_hash is not None and str(reported_unique_hash) != unique_hash:
        raise ValueError(
            "global terminal unique-profile digest mismatch: "
            f"{reported_unique_hash} != {unique_hash}"
        )

    assigned = {profile for profile, _mixture in assignments}
    outside = assigned - expected
    if outside:
        sample = min(outside)
        raise ValueError(f"discrete point assignment lies outside terminal cells: {sample}")
    missing = expected - assigned
    if missing:
        sample = min(missing)
        raise ValueError(
            f"terminal lattice cover misses {len(missing)} profiles; first is {sample}"
        )
    if assignments and not terminals:
        raise ValueError("discrete point assignments exist without terminal triangles")
    return {
        "terminal_cells": len(terminals),
        "terminal_lattice_occurrences": occurrences,
        "terminal_unique_lattice_profiles": len(expected),
        "discrete_assignment_rows_after_deduplication": len(assignments),
        "discrete_assigned_unique_profiles": len(assigned),
        "terminal_unique_profiles_sha256": unique_hash,
        "terminal_cell_lattice_sha256": cell_lattice_hash,
        "terminal_cell_lattice_counts": cell_counts,
        "terminal_cell_profile_sha256": cell_hashes,
    }


def verify_exact_geometry(ledger: dict[str, Any]) -> dict[str, Any]:
    """Verify the canonical hull partition, including hybrid terminal cells."""

    triangles = ledger.get("triangles") or ledger.get("triangle_to_witness_ledger") or ledger.get("cells", [])
    terminals = _terminal_triangles(ledger)
    partition = list(triangles) + terminals
    if not partition:
        raise ValueError("complete certificate needs explicit triangular geometry")
    for row in triangles:
        if any(
            key in row
            for key in (
                "vertex_mixture_weights",
                "vertex_mixture_witness_references",
                "vertex_mixtures",
            )
        ):
            raise ValueError("triangle geometry contains vertex-specific mixture weights")
        references = row.get("vertex_witnesses", row.get("witnesses"))
        if references is not None and len(set(map(str, references))) != 1:
            raise ValueError("triangle geometry assigns more than one fixed witness")
    total = atom_count(GROUP_BITS)
    pentagon = (
        (0, total, 0),
        (0, 0, total),
        (total - 11, 0, 11),
        (total - 11, 1, 10),
        (total - 21, 21, 0),
    )
    full_support = (
        (total - 20, 19, 1),
        (1, total - 2, 1),
        (1, 1, total - 2),
        (total - 11, 1, 10),
    )
    domain = str(ledger.get("geometry_domain", ledger.get("domain", ""))).lower()
    if not domain:
        floor = int(ledger.get("parameters", {}).get("interior_floor", 1))
        domain = "clipped_integer_hull" if floor == 0 else "full_support_plus_edges"
    if domain in {"clipped_integer_hull", "pentagon", "full_pentagon"}:
        result = _mesh_covers_polygon(partition, pentagon)
        result["domain"] = "clipped_integer_hull_pentagon"
    elif domain in {"full_support_hull", "full_support_plus_edges", ""}:
        result = _mesh_covers_polygon(partition, full_support)
        segments = ledger.get("segments", ledger.get("edge_segments", []))
        result.update(_edge_segments_cover(segments))
        result["domain"] = "full_support_hull_plus_three_edges"
    else:
        raise ValueError(f"unknown triangle geometry domain {domain!r}")
    result["covered_triangles"] = len(triangles)
    result["terminal_triangles"] = len(terminals)
    return result


def _cover_is_complete(ledger: dict[str, Any]) -> bool:
    if _terminal_triangles(ledger):
        # A hybrid artifact deliberately lacks a continuous cover.  Its
        # producer assertion is instead about all integer profiles, and is
        # accepted only after the verifier independently proves the partition
        # and exhausts every terminal-cell lattice point.
        return bool(
            ledger.get(
                "complete_global_integer_cover",
                ledger.get("complete_integer_cover", False),
            )
        )
    for key in ("complete_cover", "coverage_complete", "all_profiles_covered"):
        if key in ledger:
            return bool(ledger[key])
    status = str(ledger.get("status", "")).upper()
    return "COMPLETE" in status and "INCOMPLETE" not in status


def _normalization_worker(profile: tuple[int, int, int]):
    return profile, outward_normalization(GROUP_BITS, list(profile))


def precompute_outward_normalizations(
    profiles: list[tuple[int, int, int]], workers: int
) -> None:
    missing = sorted(
        set(profiles).difference(_OUTWARD_NORMALIZATION_CACHE)
    )
    if not missing:
        return
    if workers == 1:
        for index, profile in enumerate(missing, 1):
            _OUTWARD_NORMALIZATION_CACHE[profile] = outward_normalization(
                GROUP_BITS, list(profile)
            )
            if index % 1000 == 0 or index == len(missing):
                print(
                    f"outward normalizations={index}/{len(missing)}",
                    flush=True,
                )
        return
    with ProcessPoolExecutor(max_workers=min(workers, len(missing))) as executor:
        completed = 0
        for profile, value in executor.map(
            _normalization_worker,
            missing,
            chunksize=64,
        ):
            _OUTWARD_NORMALIZATION_CACHE[profile] = value
            completed += 1
            if completed % 1000 == 0 or completed == len(missing):
                print(
                    f"outward normalizations={completed}/{len(missing)}",
                    flush=True,
                )


def evaluate_vertex(profile: tuple[int, int, int], hardened: dict[str, Any]) -> Interval:
    value = hardened["constant"]
    for count, charge in zip(profile, hardened["charges"]):
        if not count:
            continue
        if charge is None:
            raise ValueError(
                f"witness {hardened['name']!r} has zero fugacity on active class"
            )
        value = value - charge.times_int(count)
    if hardened.get("subtract_normalization", True):
        normalization = _OUTWARD_NORMALIZATION_CACHE.get(profile)
        if normalization is None:
            normalization = outward_normalization(GROUP_BITS, list(profile))
            _OUTWARD_NORMALIZATION_CACHE[profile] = normalization
        value = value - normalization
    return value


def evaluate_mixture(
    profile: tuple[int, int, int],
    mixture: tuple[tuple[str, Fraction], ...],
    hardened: dict[str, dict[str, Any]],
) -> tuple[Interval, list[tuple[str, Fraction, Interval]]]:
    """Outward-evaluate every component, then form the fixed weighted sum."""

    total = Interval.exact(0)
    components = []
    for name, weight in mixture:
        # Evaluation is intentionally unconditional for every positive-weight
        # component.  A zero fugacity on any active class rejects the entire
        # mixture at this vertex instead of silently dropping that component.
        value = evaluate_vertex(profile, hardened[name])
        weight_interval = Interval.exact(weight.numerator) / Interval.exact(weight.denominator)
        total = total + value * weight_interval
        components.append((name, weight, value))
    return total, components


def harden_full_bijection() -> dict[str, Any]:
    constant = (
        Interval.exact(K)
        + log2_int(11).times_int(N)
        - log2_int(10).times_int(N - D)
    )
    zero = Interval.exact(0)
    return {
        "name": "full_bijection",
        "constant": constant,
        "charges": (zero, zero, zero),
        "subtract_normalization": True,
        "report": {
            "name": "full_bijection",
            "kind": "bijection",
            "constant_log2_interval": [str(constant.lo), str(constant.hi)],
        },
    }


def _source_paths(ledger_path: Path, ledger: dict[str, Any]) -> list[Path]:
    result = []
    for row in ledger.get("sources", []):
        raw = Path(str(row.get("path", "")))
        candidates = (raw, ledger_path.parent / raw, Path.cwd() / raw)
        resolved = next((path for path in candidates if path.exists()), None)
        if resolved is None:
            raise ValueError(f"cannot resolve ledger witness source {raw}")
        expected = row.get("sha256")
        if expected is not None:
            actual = hashlib.sha256(resolved.read_bytes()).hexdigest()
            if actual != str(expected):
                raise ValueError(
                    f"ledger witness source digest mismatch for {resolved}: "
                    f"{actual} != {expected}"
                )
        result.append(resolved)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--witness", type=Path, action="append", default=[])
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--checkpoint-dir", type=Path)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.iterations <= 0:
        raise SystemExit("triangle ledger certificate: iterations must be positive")
    if args.workers <= 0:
        raise SystemExit("triangle ledger certificate: workers must be positive")
    self_check()
    ledger = json.loads(args.ledger.read_text(encoding="utf-8"))
    if int(ledger.get("group_bits", GROUP_BITS)) != GROUP_BITS:
        raise SystemExit("triangle ledger certificate: only g=2 is supported")
    hybrid = bool(_terminal_triangles(ledger) or ledger.get("discrete_point_assignments"))
    producer_complete = _cover_is_complete(ledger)
    if not producer_complete and not args.allow_incomplete:
        raise SystemExit(
            "triangle ledger certificate: producer did not assert complete integer coverage; "
            "use --allow-incomplete only for smoke testing"
        )

    geometry_report = None
    if producer_complete or hybrid:
        geometry_report = verify_exact_geometry(ledger)
    continuous_assignments = []
    try:
        continuous_assignments = extract_assignments(ledger)
    except ValueError as error:
        if "no witness-at-vertex assignments" not in str(error) or not hybrid:
            if not args.allow_incomplete or "no witness-at-vertex assignments" not in str(error):
                raise
    continuous_cells = extract_triangle_cells(ledger)
    discrete_assignments = extract_discrete_assignments(ledger)
    discrete_report = None
    if hybrid:
        discrete_report = verify_discrete_terminal_cover(ledger, discrete_assignments)
    assignments = [
        ("affine_vertex", profile, mixture)
        for profile, mixture in continuous_assignments
    ] + [
        ("terminal_lattice_point", profile, mixture)
        for profile, mixture in discrete_assignments
    ]
    assignments = sorted(set(assignments))
    complete = bool(producer_complete and (not hybrid or discrete_report is not None))
    if not assignments:
        if not args.allow_incomplete:
            raise ValueError("triangle ledger contains no certifiable assignments")
        note = "triangle ledger contains no witness-at-vertex or discrete assignments"
        if hybrid and discrete_report is not None:
            note = "hybrid ledger has no certifiable assignments"
        report = {
            "status": "OUTWARD_G2_TRIANGULAR_LEDGER_INCOMPLETE_SMOKE",
            "group_bits": GROUP_BITS,
            "cover_complete": False,
            "integer_profile_cover_complete": False,
            "used_witnesses": 0,
            "used_witness_vertex_inequalities": 0,
            "used_mixture_vertex_inequalities": 0,
            "used_discrete_point_inequalities": 0,
            "passed": False,
            "note": note,
        }
        rendered = json.dumps(report, indent=2, sort_keys=True)
        print(rendered)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered + "\n", encoding="utf-8")
        return
    witness_paths = list(args.witness) or _source_paths(args.ledger, ledger)
    witnesses = load_witnesses(witness_paths)
    used_names = sorted(
        {name for _kind, _profile, mixture in assignments for name, _weight in mixture}
    )
    missing = [
        name for name in used_names if name != "full_bijection" and name not in witnesses
    ]
    if missing:
        raise SystemExit(
            "triangle ledger certificate: missing selected witnesses: " + ", ".join(missing)
        )
    hardened = harden_selected_witnesses(
        used_names,
        witnesses,
        args.iterations,
        args.workers,
        args.checkpoint_dir,
    )
    normalization_profiles = [
        profile
        for _kind, profile, mixture in assignments
        if any(hardened[name].get("subtract_normalization", True) for name, _weight in mixture)
    ]
    precompute_outward_normalizations(normalization_profiles, args.workers)

    digest = hashlib.sha256()
    maximum: Interval | None = None
    worst_pair = None
    worst_components = None
    usage = Counter()
    mixture_usage = Counter()
    component_evaluations = 0
    assignment_usage = Counter()
    evaluated: dict[
        tuple[tuple[int, int, int], tuple[tuple[str, Fraction], ...]], Interval
    ] = {}
    for assignment_kind, profile, mixture in assignments:
        value, components = evaluate_mixture(profile, mixture, hardened)
        evaluated[(profile, mixture)] = value
        assignment_usage[assignment_kind] += 1
        mixture_record = [
            {"witness": name, "weight": _weight_text(weight)}
            for name, weight in mixture
        ]
        component_record = []
        for name, weight, component_value in components:
            usage[name] += 1
            component_evaluations += 1
            component_record.append(
                [name, _weight_text(weight), str(component_value.lo), str(component_value.hi)]
            )
        mixture_key = json.dumps(mixture_record, separators=(",", ":"), sort_keys=True)
        mixture_usage[mixture_key] += 1
        digest.update(
            json.dumps(
                [
                    assignment_kind,
                    list(profile),
                    mixture_record,
                    component_record,
                    str(value.lo),
                    str(value.hi),
                ],
                separators=(",", ":"),
                sort_keys=True,
            ).encode()
        )
        if maximum is None or value.hi > maximum.hi:
            maximum = value
            worst_pair = (assignment_kind, profile, mixture)
            worst_components = components
    assert maximum is not None and worst_pair is not None and worst_components is not None

    profile_count_interval = log2_int(profile_count(GROUP_BITS, N))
    global_max_union = maximum + profile_count_interval

    # The legacy allocation charged the worst certified profile against every
    # profile in the simplex.  A much sharper bound follows directly from the
    # same convex-cell certificate: on each triangle the fixed mixture is
    # bounded by its maximum vertex value, so multiply that bound only by the
    # exact number of lattice profiles in the closed triangle.  Closed-cell
    # boundary duplicates are deliberately retained as a safe overcount.
    # Terminal cells are represented by their globally deduplicated singleton
    # assignments and are therefore summed one profile at a time.
    aggregation_digest = hashlib.sha256()
    cell_terms: list[Interval] = []
    cell_term_reports = []
    for cell_id, vertices, mixture, lattice_count in continuous_cells:
        vertex_values = []
        for vertex in vertices:
            value = evaluated.get((vertex, mixture))
            if value is None:
                value, _components = evaluate_mixture(vertex, mixture, hardened)
            vertex_values.append(value)
        cell_maximum = Interval(
            max(value.lo for value in vertex_values),
            max(value.hi for value in vertex_values),
        )
        contribution = cell_maximum + log2_int(lattice_count)
        mixture_record = [
            {"witness": name, "weight": _weight_text(weight)}
            for name, weight in mixture
        ]
        record = {
            "assignment_kind": "closed_convex_triangle",
            "cell_id": cell_id,
            "vertices": [list(vertex) for vertex in vertices],
            "mixture": mixture_record,
            "lattice_count": lattice_count,
            "maximum_log2_interval": [str(cell_maximum.lo), str(cell_maximum.hi)],
            "contribution_log2_interval": [str(contribution.lo), str(contribution.hi)],
        }
        aggregation_digest.update(
            json.dumps(record, separators=(",", ":"), sort_keys=True).encode()
        )
        cell_terms.append(contribution)
        cell_term_reports.append(record)

    for profile, mixture in discrete_assignments:
        value = evaluated[(profile, mixture)]
        record = {
            "assignment_kind": "terminal_lattice_point",
            "profile": list(profile),
            "mixture": [
                {"witness": name, "weight": _weight_text(weight)}
                for name, weight in mixture
            ],
            "lattice_count": 1,
            "contribution_log2_interval": [str(value.lo), str(value.hi)],
        }
        aggregation_digest.update(
            json.dumps(record, separators=(",", ":"), sort_keys=True).encode()
        )
        cell_terms.append(value)
        cell_term_reports.append(record)

    if continuous_cells:
        union = _log2_sum_exp(cell_terms)
        aggregation_method = "closed_cell_lattice_weighted_logsumexp"
    else:
        # Preserve support for the legacy flattened smoke-test schema, which
        # has no cell geometry from which to derive local multiplicities.
        union = global_max_union
        aggregation_method = "global_max_times_profile_count_fallback"
    security_margin = -union.hi
    slack_beyond_target = Decimal(-40) - union.hi
    passed = bool(complete and slack_beyond_target >= 0)
    cell_term_reports.sort(
        key=lambda row: Decimal(row["contribution_log2_interval"][1]), reverse=True
    )
    report = {
        "status": (
            (
                "OUTWARD_CERTIFIED_G2_HYBRID_INTEGER_2^-40_LEDGER"
                if hybrid
                else "OUTWARD_CERTIFIED_G2_TRIANGULAR_2^-40_LEDGER"
            )
            if passed
            else "OUTWARD_G2_TRIANGULAR_LEDGER_DID_NOT_CLOSE"
        ),
        "group_bits": GROUP_BITS,
        "cover_complete": complete,
        "integer_profile_cover_complete": complete,
        "producer_integer_cover_assertion": producer_complete if hybrid else None,
        "exact_geometry": geometry_report,
        "discrete_terminal_cover": discrete_report,
        "used_witnesses": len(used_names),
        "used_fixed_mixtures": len(mixture_usage),
        "used_mixture_vertex_inequalities": len(continuous_assignments),
        "used_discrete_point_inequalities": len(discrete_assignments),
        "used_witness_vertex_inequalities": component_evaluations,
        "assignment_usage": dict(sorted(assignment_usage.items())),
        "witness_usage": dict(sorted(usage.items())),
        "mixture_usage": dict(sorted(mixture_usage.items())),
        "inequality_sha256": digest.hexdigest(),
        "worst_assignment_kind": worst_pair[0],
        "worst_profile": list(worst_pair[1]),
        "worst_witness": worst_pair[2][0][0] if len(worst_pair[2]) == 1 else None,
        "worst_mixture": [
            {"witness": name, "weight": _weight_text(weight)}
            for name, weight in worst_pair[2]
        ],
        "worst_component_log2_intervals": [
            {
                "witness": name,
                "weight": _weight_text(weight),
                "interval": [str(component_value.lo), str(component_value.hi)],
            }
            for name, weight, component_value in worst_components
        ],
        "maximum_branch_log2_interval": [str(maximum.lo), str(maximum.hi)],
        "profile_count": profile_count(GROUP_BITS, N),
        "profile_count_log2_interval": [str(profile_count_interval.lo), str(profile_count_interval.hi)],
        "legacy_global_max_union_log2_interval": [
            str(global_max_union.lo),
            str(global_max_union.hi),
        ],
        "aggregation_method": aggregation_method,
        "aggregation_closed_cell_boundary_duplicates_retained": bool(continuous_cells),
        "aggregation_convex_cells": len(continuous_cells),
        "aggregation_terminal_unique_profiles": len(discrete_assignments),
        "aggregation_terms": len(cell_terms),
        "aggregation_sha256": aggregation_digest.hexdigest(),
        "dominant_aggregation_terms": cell_term_reports[:20],
        "union_log2_interval": [str(union.lo), str(union.hi)],
        "certified_margin_bits": str(security_margin),
        "certified_slack_beyond_40_bits": str(slack_beyond_target),
        "passed": passed,
        "arithmetic": (
            "exact Fraction local transportation and outer polynomials; guarded "
            "directed Decimal logarithm/exponential intervals"
        ),
        "witness_reports": [hardened[name]["report"] for name in used_names],
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    if complete and not passed:
        raise SystemExit("triangle ledger certificate: outward 40-bit union did not close")


if __name__ == "__main__":
    main()
