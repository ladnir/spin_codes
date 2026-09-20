#!/usr/bin/env python3
"""Probe genus-two MacWilliams LP bounds for a small BA constituent.

The unknown variables are the complete pair-type counts of a binary linear
code.  The LP fixes the ordinary spectrum and every rank-one pair count.  It
then imposes either nonnegativity of the dual genus-two transform or exact
self-duality.  The objective is the expected square of one output-shell count
after common-interleaver accumulator stages.

This is a floating small-length probe.  It tests whether genus-two constraints
can replace explicit pair enumeration; it is not a length-512 certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, linprog

from analyze_accumulator_pair_chain_small import compositions4, direct_sum_words
from probe_triangle_holder_pair_bound_small import (
    direct_sum_spectrum,
    one_word_transition,
)
from verify_accumulator_pair_type_kernel import kernel_counts, multinomial


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "genus2_macwilliams_lp_B8_probe.json"

HADAMARD = np.asarray(
    [
        [1, 1, 1, 1],
        [1, -1, 1, -1],
        [1, 1, -1, -1],
        [1, -1, -1, 1],
    ],
    dtype=np.int64,
)


def transform_column(kind: tuple[int, int, int, int]) -> dict[tuple[int, int, int, int], int]:
    polynomial: dict[tuple[int, int, int, int], int] = {(0, 0, 0, 0): 1}
    for symbol, multiplicity in enumerate(kind):
        for _ in range(multiplicity):
            following: dict[tuple[int, int, int, int], int] = {}
            for exponent, coefficient in polynomial.items():
                for output_symbol in range(4):
                    next_exponent = list(exponent)
                    next_exponent[output_symbol] += 1
                    key = tuple(next_exponent)
                    following[key] = following.get(key, 0) + coefficient * int(
                        HADAMARD[symbol, output_symbol]
                    )
            polynomial = following
    return polynomial


def macwilliams_matrix(types: list[tuple[int, int, int, int]]) -> np.ndarray:
    index = {kind: position for position, kind in enumerate(types)}
    matrix = np.zeros((len(types), len(types)), dtype=float)
    for column, kind in enumerate(types):
        for output, coefficient in transform_column(kind).items():
            matrix[index[output], column] = coefficient
    return matrix


def pair_transition(types: list[tuple[int, int, int, int]]) -> np.ndarray:
    index = {kind: position for position, kind in enumerate(types)}
    transition = np.zeros((len(types), len(types)), dtype=float)
    for row, kind in enumerate(types):
        denominator = multinomial(kind)
        for output, count in kernel_counts(kind).items():
            transition[row, index[output]] = count / denominator
    return transition


def pair_counts(words: list[int], length: int, types: list[tuple[int, int, int, int]]) -> np.ndarray:
    index = {kind: position for position, kind in enumerate(types)}
    counts = np.zeros(len(types), dtype=float)
    mask = (1 << length) - 1
    for first in words:
        for second in words:
            n11 = (first & second).bit_count()
            n10 = (first & (~second & mask)).bit_count()
            n01 = ((~first & mask) & second).bit_count()
            kind = (length - n01 - n10 - n11, n01, n10, n11)
            counts[index[kind]] += 1.0
    return counts


def marginal_rows(
    types: list[tuple[int, int, int, int]], spectrum: list[int], messages: int
) -> tuple[np.ndarray, np.ndarray]:
    rows = []
    rhs = []
    for character in range(3):
        for weight in range(len(spectrum)):
            row = np.zeros(len(types))
            for index, kind in enumerate(types):
                weights = (kind[2] + kind[3], kind[1] + kind[3], kind[1] + kind[2])
                if weights[character] == weight:
                    row[index] = 1.0
            rows.append(row / messages)
            rhs.append(float(spectrum[weight]))
    return np.asarray(rows), np.asarray(rhs)


def rank_one_bounds(
    types: list[tuple[int, int, int, int]], spectrum: list[int]
) -> Bounds:
    lower = np.zeros(len(types))
    upper = np.full(len(types), np.inf)
    for index, kind in enumerate(types):
        nonzero_symbols = sum(value > 0 for value in kind[1:])
        if nonzero_symbols <= 1:
            weight = sum(kind[1:])
            lower[index] = float(spectrum[weight])
            upper[index] = float(spectrum[weight])
    return Bounds(lower, upper)


def solve_objective(
    objective: np.ndarray,
    transform: np.ndarray,
    messages: int,
    marginal_matrix: np.ndarray,
    marginal_rhs: np.ndarray,
    bounds: Bounds,
    mode: str,
) -> tuple[float, str]:
    equalities = [marginal_matrix]
    equality_rhs = [marginal_rhs]
    inequalities = None
    inequality_rhs = None
    scaled_transform = transform / (messages * messages)
    if mode == "self-dual":
        equalities.append(scaled_transform - np.eye(transform.shape[0]))
        equality_rhs.append(np.zeros(transform.shape[0]))
    elif mode == "formal-dual":
        inequalities = -scaled_transform
        inequality_rhs = np.zeros(transform.shape[0])
    result = linprog(
        -objective,
        A_ub=inequalities,
        b_ub=inequality_rhs,
        A_eq=np.vstack(equalities),
        b_eq=np.concatenate(equality_rhs),
        bounds=list(zip(bounds.lb, bounds.ub)),
        method="highs",
        options={"presolve": True},
    )
    if not result.success:
        return math.nan, str(result.message)
    return -float(result.fun), "OPTIMAL"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blocks", type=int, choices=(1, 2), default=1)
    parser.add_argument("--maximum-stages", type=int, default=3)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    length = 8 * args.blocks
    dimension = 4 * args.blocks
    messages = 1 << dimension
    spectrum_array = direct_sum_spectrum(args.blocks)
    spectrum = [int(value) for value in spectrum_array]
    types = compositions4(length)
    words = direct_sum_words(args.blocks)
    actual_counts = pair_counts(words, length, types)

    print(f"building_macwilliams length={length} types={len(types)}", flush=True)
    transform = macwilliams_matrix(types)
    involution_error = float(
        np.max(np.abs(transform @ transform - (4**length) * np.eye(len(types))))
    )
    marginal_matrix, marginal_rhs = marginal_rows(types, spectrum, messages)
    bounds = rank_one_bounds(types, spectrum)
    actual_transform_error = float(
        np.max(np.abs(transform @ actual_counts / (messages * messages) - actual_counts))
    )

    print("building_pair_transition", flush=True)
    transition = pair_transition(types)
    one_word = one_word_transition(length)
    pair_power = np.eye(len(types))
    one_power = np.eye(length + 1)
    rows = []
    for stages in range(1, args.maximum_stages + 1):
        pair_power = pair_power @ transition
        one_power = one_power @ one_word
        expected_spectrum = np.asarray(spectrum, dtype=float) @ one_power
        for target_weight in range(1, length + 1):
            terminal = np.asarray(
                [
                    float(
                        kind[2] + kind[3] == target_weight
                        and kind[1] + kind[3] == target_weight
                    )
                    for kind in types
                ]
            )
            shell_objective = pair_power @ terminal
            mean = float(expected_spectrum[target_weight])
            actual_second = float(actual_counts @ shell_objective)
            actual_variance = max(0.0, actual_second - mean * mean)
            row: dict[str, object] = {
                "accumulator_stages": stages,
                "target_weight": target_weight,
                "mean": mean,
                "actual_variance_over_mean": actual_variance / mean if mean else None,
            }
            for mode in ("marginal", "formal-dual", "self-dual"):
                maximum_second, status = solve_objective(
                    shell_objective,
                    transform,
                    messages,
                    marginal_matrix,
                    marginal_rhs,
                    bounds,
                    mode,
                )
                maximum_variance = max(0.0, maximum_second - mean * mean)
                row[f"{mode}_status"] = status
                row[f"{mode}_variance_over_mean_upper"] = (
                    maximum_variance / mean if mean and math.isfinite(maximum_second) else None
                )
            rows.append(row)
            print(
                f"stage,{stages},weight,{target_weight},actual,{row['actual_variance_over_mean']},"
                f"marginal,{row['marginal_variance_over_mean_upper']},"
                f"formal,{row['formal-dual_variance_over_mean_upper']},"
                f"selfdual,{row['self-dual_variance_over_mean_upper']}",
                flush=True,
            )

    worst = {}
    for mode in ("actual", "marginal", "formal-dual", "self-dual"):
        field = (
            "actual_variance_over_mean"
            if mode == "actual"
            else f"{mode}_variance_over_mean_upper"
        )
        finite = [row for row in rows if row[field] is not None]
        worst[mode] = max(finite, key=lambda row: float(row[field])) if finite else None

    payload = {
        "schema": "genus2-macwilliams-lp-small-probe-v1",
        "status": "EXACT_COMBINATORIAL_INPUT_WITH_BINARY64_LP",
        "parameters": {
            "base_code": f"direct sum of {args.blocks} RM(1,3) [8,4,4] blocks",
            "length": length,
            "dimension": dimension,
            "pair_types": len(types),
        },
        "checks": {
            "macwilliams_involution_maximum_error": involution_error,
            "actual_self_dual_transform_maximum_error": actual_transform_error,
        },
        "worst": worst,
        "rows": rows,
        "limitations": [
            "The pair kernel and MacWilliams coefficients are generated from exact integers.",
            "HiGHS solves the LP in binary64; the reported optima are not outward certificates.",
            "The self-dual constraints apply to the RM(1,3) test code, not to the formally self-dual EBCH128 code.",
            "The small-length result does not prove the EBCH128 or length-512 target.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"checks": payload["checks"], "worst": worst}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
