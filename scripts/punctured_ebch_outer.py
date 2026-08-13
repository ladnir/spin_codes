#!/usr/bin/env python3
"""Exact one-coordinate puncturing data for the committed extended BCH code.

The committed [128,64,22] generator is checked to be invariant under cyclic
scaling on the 127 nonzero field coordinates and under x -> x+1 on F_128.
These permutations generate a transitive affine action, so every coordinate
has the same weight-conditioned occupancy.  The punctured [127,64,>=21]
enumerator is therefore

    B_j = ((128-j) A_j + (j+1) A_{j+1}) / 128.

The module also provides truncated heterogeneous direct-sum coefficients and
floating Cauchy bounds for products of full and punctured local enumerators.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from analyze_nosinger_bch_kernel import generator_rows
from certificate_spectra import load_csv_spectrum
from certify_rm_outer_prefix_exact import (
    direct_sum_coefficients,
    truncated_convolution,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_SPECTRUM = ROOT / "ebch128_64_spectrum.csv"
FULL_LENGTH = 128
PUNCTURED_LENGTH = 127
DIMENSION = 64
MINIMUM_DISTANCE = 22


def _row_basis(rows: list[int]) -> dict[int, int]:
    basis: dict[int, int] = {}
    for row in rows:
        value = row
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return basis


def _in_span(value: int, basis: dict[int, int]) -> bool:
    while value:
        pivot = value.bit_length() - 1
        if pivot not in basis:
            return False
        value ^= basis[pivot]
    return True


def _permute_word(word: int, permutation: list[int]) -> int:
    result = 0
    for source, target in enumerate(permutation):
        if (word >> source) & 1:
            result |= 1 << target
    return result


def verify_coordinate_transitivity() -> None:
    """Machine-check the two generator permutations used for transitivity."""

    rows = list(generator_rows())
    basis = _row_basis(rows)
    if len(basis) != DIMENSION:
        raise SystemExit("punctured EBCH: committed generator does not have rank 64")

    # Coordinate i=0,...,126 is labelled alpha^i in GF(128), while coordinate
    # 127 is labelled zero.  The primitive polynomial is x^7+x+1 (0x83).
    powers = [1]
    for _ in range(1, 127):
        value = powers[-1] << 1
        if value & 0x80:
            value ^= 0x83
        powers.append(value)
    if len(set(powers)) != 127:
        raise SystemExit("punctured EBCH: x is not primitive modulo 0x83")
    logarithm = {value: exponent for exponent, value in enumerate(powers)}

    scaling = [(coordinate + 1) % 127 for coordinate in range(127)] + [127]
    translation: list[int] = []
    for coordinate in range(127):
        translated = powers[coordinate] ^ 1
        translation.append(127 if translated == 0 else logarithm[translated])
    translation.append(0)

    for name, permutation in (("scaling", scaling), ("translation", translation)):
        if sorted(permutation) != list(range(FULL_LENGTH)):
            raise SystemExit(f"punctured EBCH: {name} is not a permutation")
        if not all(_in_span(_permute_word(row, permutation), basis) for row in rows):
            raise SystemExit(f"punctured EBCH: generator is not invariant under {name}")


def full_spectrum(path: Path = DEFAULT_SPECTRUM) -> list[int]:
    loaded = load_csv_spectrum(
        path,
        name="extended BCH [128,64,22] spectrum",
        length=FULL_LENGTH,
        dimension=DIMENSION,
        minimum_distance=MINIMUM_DISTANCE,
    )
    result = [0] * (FULL_LENGTH + 1)
    for weight, count in loaded:
        result[weight] = count
    return result


def punctured_spectrum(path: Path = DEFAULT_SPECTRUM) -> list[int]:
    verify_coordinate_transitivity()
    full = full_spectrum(path)
    punctured = [0] * (PUNCTURED_LENGTH + 1)
    for weight in range(PUNCTURED_LENGTH + 1):
        numerator = (
            (FULL_LENGTH - weight) * full[weight]
            + (weight + 1) * full[weight + 1]
        )
        if numerator % FULL_LENGTH:
            raise SystemExit(
                f"punctured EBCH: nonintegral weight-{weight} coefficient"
            )
        punctured[weight] = numerator // FULL_LENGTH
    if sum(punctured) != 1 << DIMENSION:
        raise SystemExit("punctured EBCH: enumerator does not sum to 2^64")
    if any(
        punctured[weight] != punctured[PUNCTURED_LENGTH - weight]
        for weight in range(PUNCTURED_LENGTH + 1)
    ):
        raise SystemExit("punctured EBCH: enumerator is not symmetric")
    if next(weight for weight, count in enumerate(punctured) if weight and count) < 21:
        raise SystemExit("punctured EBCH: minimum-distance check failed")
    return punctured


def heterogeneous_coefficients(
    *,
    full_blocks: int,
    punctured_blocks: int,
    h_max: int,
    spectrum_path: Path = DEFAULT_SPECTRUM,
) -> list[int]:
    full = full_spectrum(spectrum_path)[: h_max + 1]
    punctured = punctured_spectrum(spectrum_path)[: h_max + 1]
    if len(full) < h_max + 1:
        full.extend([0] * (h_max + 1 - len(full)))
    if len(punctured) < h_max + 1:
        punctured.extend([0] * (h_max + 1 - len(punctured)))
    full_sum = direct_sum_coefficients(full, blocks=full_blocks, h_max=h_max)
    punctured_sum = direct_sum_coefficients(
        punctured, blocks=punctured_blocks, h_max=h_max
    )
    return truncated_convolution(full_sum, punctured_sum, h_max)


def optimize_heterogeneous_log2_many(
    message_weights: np.ndarray,
    *,
    full_blocks: int,
    punctured_blocks: int,
    spectrum_path: Path = DEFAULT_SPECTRUM,
) -> np.ndarray:
    """Floating Cauchy coefficient bounds for A(z)^F B(z)^P."""

    full = full_spectrum(spectrum_path)
    punctured = punctured_spectrum(spectrum_path)
    full_weights = np.flatnonzero(full).astype(float)
    punctured_weights = np.flatnonzero(punctured).astype(float)
    full_logs = np.log(np.array([full[int(weight)] for weight in full_weights]))
    punctured_logs = np.log(
        np.array([punctured[int(weight)] for weight in punctured_weights])
    )
    result = np.empty(len(message_weights), dtype=float)

    def moments(log_poles: np.ndarray, weights: np.ndarray, logs: np.ndarray):
        exponents = logs[:, None] + weights[:, None] * log_poles
        maxima = np.max(exponents, axis=0)
        terms = np.exp(exponents - maxima)
        normalizers = np.sum(terms, axis=0)
        means = np.sum(weights[:, None] * terms, axis=0) / normalizers
        log_values = maxima + np.log(normalizers)
        return means, log_values

    for start in range(0, len(message_weights), 65536):
        stop = min(start + 65536, len(message_weights))
        targets = np.asarray(message_weights[start:stop], dtype=float)
        low = np.full(len(targets), math.log(1e-12))
        high = np.full(len(targets), math.log(1.0 - 1e-12))
        for _ in range(70):
            log_poles = (low + high) / 2.0
            full_means, _ = moments(log_poles, full_weights, full_logs)
            punctured_means, _ = moments(
                log_poles, punctured_weights, punctured_logs
            )
            means = full_blocks * full_means + punctured_blocks * punctured_means
            go_right = means < targets
            low = np.where(go_right, log_poles, low)
            high = np.where(go_right, high, log_poles)
        log_poles = (low + high) / 2.0
        _, full_values = moments(log_poles, full_weights, full_logs)
        _, punctured_values = moments(log_poles, punctured_weights, punctured_logs)
        result[start:stop] = (
            full_blocks * full_values
            + punctured_blocks * punctured_values
            - targets * log_poles
        ) / math.log(2.0)
    return result


if __name__ == "__main__":
    punctured = punctured_spectrum()
    print("punctured_ebch_status,EXACT_FROM_TRANSITIVE_COMMITTED_GENERATOR")
    print(f"length,{PUNCTURED_LENGTH}")
    print(f"dimension,{DIMENSION}")
    print(f"minimum_distance,{next(i for i, value in enumerate(punctured) if i and value)}")
    print(f"cardinality,{sum(punctured)}")
