#!/usr/bin/env python3
"""Finite one-active diagnostic for BA-3 plus bitwise RandomStepConv.

RandomStepConv has one input bit, one output bit, and ``memory_bits`` state
bits.  Setup samples one independent unrestricted binary linear map on the
input-plus-state space at every bit position.  This program evaluates the
exact two-state transfer for one active outer row after the existing
row-coordinate and transposed-region permutations.

The BA spectrum and numerical optimization use nearest binary64.  The output
is a diagnostic, not an outward certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
REPO_ROOT = WORKSTREAM.parents[1]
sys.path.insert(0, str(WORKSTREAM))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from analyze_golay_ba_rm2sub_joint import expected_ba_log_spectrum  # noqa: E402
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (  # noqa: E402
    weight_conditioned_log_moments,
)


B = 240
K = 120
L = 8832
N = B * L
D = (11 * N) // 100
LOG2 = math.log(2.0)
DEFAULT_OUTPUT = WORKSTREAM / "ba240_randomstepconv_g1_s19_q1_d11.json"


def matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return left @ right


def matrix_power(matrix: np.ndarray, exponent: int) -> np.ndarray:
    result = np.eye(2)
    factor = matrix.copy()
    while exponent:
        if exponent & 1:
            result = matmul(result, factor)
        exponent >>= 1
        if exponent:
            factor = matmul(factor, factor)
    return result


def marked_region(
    zero_step: np.ndarray, active_step: np.ndarray, positions: int
) -> tuple[np.ndarray, np.ndarray]:
    """Average one active input bit over all positions in one region."""
    powers = [np.eye(2)]
    for _ in range(positions):
        powers.append(powers[-1] @ zero_step)
    inactive = powers[positions]
    active = np.zeros((2, 2))
    for position in range(positions):
        active += powers[position] @ active_step @ powers[positions - position - 1]
    return inactive, active / positions


def step_matrices(z: float, memory_bits: int) -> tuple[np.ndarray, np.ndarray]:
    """Return exact zero-input and one-input tilted state-class kernels."""
    state_zero = math.ldexp(1.0, -memory_bits)
    output_moment = (1.0 + z) / 2.0
    terminate = state_zero * output_moment
    survive = (1.0 - state_zero) * output_moment
    zero = np.asarray(((1.0, 0.0), (terminate, survive)))
    active = np.asarray(((terminate, survive), (terminate, survive)))
    return zero, active


def evaluate(
    memory_bits: int,
    grid_min: float,
    grid_max: float,
    grid_step: float,
    lower_weight: int,
    upper_weight: int,
) -> dict[str, object]:
    spectrum = expected_ba_log_spectrum(B)
    best = np.full(B + 1, math.inf)
    best_log_surprisal = np.full(B + 1, math.nan)
    grid_count = int(math.floor((grid_max - grid_min) / grid_step + 0.5)) + 1

    for index in range(grid_count):
        log_surprisal = grid_min + index * grid_step
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero, active = step_matrices(z, memory_bits)
        inactive_region, active_region = marked_region(zero, active, L)
        moments = weight_conditioned_log_moments(
            outer_bits=B,
            zero_row=inactive_region,
            active_row=active_region,
        )
        values = moments + D * surprisal
        improved = values < best
        best[improved] = values[improved]
        best_log_surprisal[improved] = log_surprisal

    unconditional_terms = []
    rows = []
    for weight in range(1, B + 1):
        if not math.isfinite(float(spectrum[weight])):
            continue
        inner = min(0.0, float(best[weight]))
        term = math.log(L) + float(spectrum[weight]) + inner
        unconditional_terms.append(term)
        rows.append(
            {
                "outer_weight": weight,
                "log2_expected_ba_multiplicity": float(spectrum[weight]) / LOG2,
                "inner_log2_upper": inner / LOG2,
                "pointwise_log2_upper": term / LOG2,
                "best_log_surprisal": float(best_log_surprisal[weight]),
            }
        )

    unconditional = float(logsumexp(np.asarray(unconditional_terms)))
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    tail_terms = [
        float(spectrum[weight])
        for weight in list(range(1, lower_weight))
        + list(range(upper_weight + 1, B + 1))
        if math.isfinite(float(spectrum[weight]))
    ]
    tail_log = float(logsumexp(np.asarray(tail_terms)))
    tail_upper = min(1.0, math.exp(tail_log))
    if tail_upper >= 1.0:
        raise ArithmeticError("tail-free conditioning event has no positive lower bound")
    good_lower = 1.0 - tail_upper
    conditioned_terms = [
        float(row["pointwise_log2_upper"]) * LOG2 - math.log(good_lower)
        for row in rows
        if lower_weight <= int(row["outer_weight"]) <= upper_weight
    ]
    conditioned = float(logsumexp(np.asarray(conditioned_terms)))
    conditioned_dominant = max(
        (
            row
            for row in rows
            if lower_weight <= int(row["outer_weight"]) <= upper_weight
        ),
        key=lambda row: float(row["pointwise_log2_upper"]),
    )
    return {
        "schema": "ba240-randomstepconv-g1-one-active-v1",
        "status": "BINARY64_DIAGNOSTIC_EXACT_INNER_TRANSFER",
        "model": {
            "equation": "(y_t,s_{t+1}) = M_t (x_t,s_t)",
            "map_distribution": (
                "independent uniform unrestricted binary linear (M+1)-by-(M+1) "
                "matrix at every bit position"
            ),
            "shared_setup": "the sampled matrices are reused for every codeword",
            "initial_state": 0,
            "terminal_state": "discarded",
        },
        "parameters": {
            "outer_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "output_bits": N,
            "distance_cutoff": D,
            "relative_distance_cutoff": D / N,
            "memory_bits": memory_bits,
            "input_bits_per_step": 1,
            "output_bits_per_step": 1,
        },
        "probability_space": (
            "one BA-3 draw for the active row, its uniform coordinate permutation, "
            "the independent uniform row permutation in every transposed region, "
            "and the independent RandomStepConv matrices"
        ),
        "claim_scope": "occupation Q=1 only",
        "unconditional_expected_ba": {
            "log2_expected_bad_upper": unconditional / LOG2,
            "margin_bits": -unconditional / LOG2,
            "dominant": dominant,
        },
        "tail_free_conditioning": {
            "allowed_outer_weights": [lower_weight, upper_weight],
            "tail_word_expectation_log2": tail_log / LOG2,
            "tail_free_probability_lower": good_lower,
            "conditional_log2_expected_bad_upper": conditioned / LOG2,
            "conditional_margin_bits": -conditioned / LOG2,
            "dominant": conditioned_dominant,
        },
        "output_tilt_grid": {
            "log_surprisal_minimum": grid_min,
            "log_surprisal_maximum": grid_max,
            "log_surprisal_step": grid_step,
            "count": grid_count,
        },
        "weight_rows": rows,
        "limitations": [
            "Occupation Q>=2 is not included.",
            "The BA spectrum is averaged over its two interleavers.",
            "The calculation uses nearest binary64 rather than outward arithmetic.",
            "No implementation-cost claim is made for the dense random step maps.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=19)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.5)
    parser.add_argument("--grid-step", type=float, default=0.02)
    parser.add_argument("--lower-weight", type=int, default=23)
    parser.add_argument("--upper-weight", type=int, default=217)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = evaluate(
        args.memory_bits,
        args.grid_min,
        args.grid_max,
        args.grid_step,
        args.lower_weight,
        args.upper_weight,
    )
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["unconditional_expected_ba"], indent=2))
    print(json.dumps(payload["tail_free_conditioning"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
