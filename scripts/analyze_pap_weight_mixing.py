#!/usr/bin/env python3
"""Global-consistency diagnostic for repeated PAP weight mixing.

The robust per-input-slice BCH caps can be mutually inconsistent.  This tool
instead combines exact boundary slices with the ordinary BCH column spectrum.
For unknown input weights s it bounds the k-round PAP density by

    P^k(t,s) <= R_t * C(64,s)/(2^64-1).

The resulting score bound converges to the Singer one-point law while retaining
the committed exact slices.  Long-double arithmetic makes this a design
diagnostic; the theorem-facing version still needs exact/outward arithmetic.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from analyze_nosinger_bch_kernel import (
    B,
    EXACT_SLICES,
    SPECTRUM,
    accumulator_weight_entries,
    exact_slice,
    generator_rows,
    load_exact_slices,
    split_rows,
)
from certificate_spectra import load_ebch128_spectrum


DEFAULT_POLES = (
    0.9800511634,
    0.9414884538,
    0.8178969571,
    0.5382199631,
    0.3204075791888,
)


def pap_weight_matrix() -> np.ndarray:
    matrix = np.zeros((B, B), dtype=np.longdouble)
    for input_weight in range(1, B + 1):
        denominator = math.comb(B, input_weight)
        for output_weight, count in accumulator_weight_entries(input_weight):
            matrix[input_weight - 1, output_weight - 1] = (
                np.longdouble(count) / denominator
            )
    return matrix


def global_consistency_bound(
    *,
    scores: dict[int, float],
    transition: np.ndarray,
    stationary: np.ndarray,
    spectrum: dict[int, int],
    exact: dict[int, dict[int, int]],
) -> tuple[float, int]:
    total = (1 << B) - 1
    known = tuple(sorted(exact))
    unknown = tuple(weight for weight in range(1, B + 1) if weight not in exact)
    stationary_average = sum(
        spectrum[weight] * scores[weight] for weight in spectrum
    ) / total
    exact_averages = {
        input_weight: sum(
            count * scores[output_weight]
            for output_weight, count in exact[input_weight].items()
        )
        / math.comb(B, input_weight)
        for input_weight in known
    }

    best = (-1.0, 0)
    for input_index in range(B):
        density = max(
            transition[input_index, weight - 1] / stationary[weight - 1]
            for weight in unknown
        )
        value = density * stationary_average
        value += sum(
            (
                transition[input_index, weight - 1]
                - density * stationary[weight - 1]
            )
            * exact_averages[weight]
            for weight in known
        )
        candidate = (float(value), input_index + 1)
        if candidate[0] > best[0]:
            best = candidate
    return best


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rounds", type=int, default=9)
    parser.add_argument("--exact-radius", type=int, default=5)
    parser.add_argument("--spectrum", type=Path, default=SPECTRUM)
    parser.add_argument("--exact-slices", type=Path, default=EXACT_SLICES)
    args = parser.parse_args()
    if not 1 <= args.rounds <= 32:
        raise SystemExit("PAP mixing diagnostic expects 1 <= rounds <= 32")
    if not 0 <= args.exact_radius <= B // 2:
        raise SystemExit("PAP mixing diagnostic: invalid exact radius")

    rows = generator_rows()
    exact_weights = tuple(range(1, args.exact_radius + 1)) + tuple(
        range(B - args.exact_radius, B + 1)
    )
    exact = {weight: exact_slice(rows, weight) for weight in exact_weights}
    for weight, histogram in load_exact_slices(args.exact_slices).items():
        exact[weight] = histogram

    loaded = load_ebch128_spectrum(args.spectrum)
    spectrum = {weight: count for weight, count in loaded if weight > 0 and count > 0}
    weights = tuple(sorted(spectrum))
    split = split_rows(weights)

    one_round = pap_weight_matrix()
    transition = np.linalg.matrix_power(one_round, args.rounds)
    stationary = np.array(
        [
            np.longdouble(math.comb(B, weight)) / ((1 << B) - 1)
            for weight in range(1, B + 1)
        ]
    )
    unknown = tuple(weight for weight in range(1, B + 1) if weight not in exact)
    density = max(
        transition[input_weight - 1, output_weight - 1]
        / stationary[output_weight - 1]
        for input_weight in range(1, B + 1)
        for output_weight in unknown
    )

    termination_best = (-1.0, 0, 0)
    for target_weight in range(B + 1):
        scores: dict[int, float] = {}
        for weight in weights:
            score = 0.0
            for _output_weight, state_weight, probability in split[weight]:
                if state_weight == 0:
                    score += float(probability)
                elif state_weight == target_weight:
                    score += float(probability) / math.comb(B, state_weight)
            scores[weight] = score
        value, input_weight = global_consistency_bound(
            scores=scores,
            transition=transition,
            stationary=stationary,
            spectrum=spectrum,
            exact=exact,
        )
        candidate = (value, target_weight, input_weight)
        if candidate[0] > termination_best[0]:
            termination_best = candidate

    print("PAP global-consistency weight-mixing diagnostic")
    print(f"rounds={args.rounds} exact_weights={tuple(sorted(exact))}")
    print(f"unknown_density_log2={math.log2(float(density)):.12f}")
    print(
        f"global_termination_log2={math.log2(termination_best[0]):.12f} "
        f"target_state_weight={termination_best[1]} "
        f"worst_input_weight={termination_best[2]}"
    )
    for pole in DEFAULT_POLES:
        scores = {
            weight: sum(
                float(probability) * pole**output_weight
                for output_weight, state_weight, probability in split[weight]
                if state_weight > 0
            )
            for weight in weights
        }
        value, input_weight = global_consistency_bound(
            scores=scores,
            transition=transition,
            stationary=stationary,
            spectrum=spectrum,
            exact=exact,
        )
        print(
            f"pole={pole:.12f} live_mgf_log2={math.log2(value):.12f} "
            f"worst_input_weight={input_weight}"
        )
    print("status=DIAGNOSTIC_ONLY_NEEDS_EXACT_OR_OUTWARD_ARITHMETIC")


if __name__ == "__main__":
    main()
