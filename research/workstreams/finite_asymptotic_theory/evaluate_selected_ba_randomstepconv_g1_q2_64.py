#!/usr/bin/env python3
"""Finite Q=2..64 diagnostic for one selected BA code and RandomStepConv.

The selection argument fixes one repeated BA-3 constituent whose complete
allowed-shell spectrum is bounded by one common factor times the expected
BA spectrum. Independent row-coordinate and region permutations then make
the support sum factor through a Bernoulli envelope. RandomStepConv uses one
fresh random linear (M+1)-by-(M+1) map at each bit position.

All calculations in this program use nearest binary64. The output is a
diagnostic, not an outward certificate.
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
from analyze_riffle_striped_random_outer import (  # noqa: E402
    log_choose,
    log_matrix_power_moment,
)


B = 240
K = 120
L = 8832
N = B * L
D = (11 * N) // 100
LOG2 = math.log(2.0)
DEFAULT_OUTPUT = (
    WORKSTREAM / "selected_ba240_randomstepconv_g1_s30_q2_64_d11.json"
)


def step_matrices(z: float, memory_bits: int) -> tuple[np.ndarray, np.ndarray]:
    q = math.ldexp(1.0, -memory_bits)
    b = (1.0 + z) / 2.0
    terminate = q * b
    survive = (1.0 - q) * b
    zero = np.asarray((1.0, 0.0, terminate, survive))
    active = np.asarray((terminate, survive, terminate, survive))
    return zero, active


def multiply_batch(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Multiply a batch of flattened two-by-two matrices by one matrix."""
    a = left[:, 0]
    b = left[:, 1]
    c = left[:, 2]
    d = left[:, 3]
    e, f, g, h = right
    return np.column_stack(
        (a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h)
    )


def uniform_candidate_regions(
    *, zero: np.ndarray, candidate: np.ndarray, maximum_q: int
) -> np.ndarray:
    """Average Q candidate positions over a uniform Q-subset of one region."""
    current = np.zeros((maximum_q + 1, 4))
    current[0] = (1.0, 0.0, 0.0, 1.0)
    for completed in range(L):
        next_maximum = min(completed + 1, maximum_q)
        updated = np.zeros_like(current)

        old_zero = current[: next_maximum + 1]
        zero_products = multiply_batch(old_zero, zero)
        degrees = np.arange(next_maximum + 1, dtype=np.float64)
        updated[: next_maximum + 1] += zero_products * (
            ((completed + 1.0) - degrees) / (completed + 1.0)
        )[:, None]

        if next_maximum:
            old_candidate = current[:next_maximum]
            candidate_products = multiply_batch(old_candidate, candidate)
            selected_degrees = np.arange(1, next_maximum + 1, dtype=np.float64)
            updated[1 : next_maximum + 1] += candidate_products * (
                selected_degrees / (completed + 1.0)
            )[:, None]
        current = updated
    return current


def selected_spectrum(
    lower_weight: int, upper_weight: int
) -> tuple[np.ndarray, float, float, float]:
    expected = expected_ba_log_spectrum(B)
    tail = [
        float(expected[weight])
        for weight in list(range(1, lower_weight))
        + list(range(upper_weight + 1, B + 1))
        if math.isfinite(float(expected[weight]))
    ]
    tail_log = float(logsumexp(np.asarray(tail)))
    tail_upper = min(1.0, math.exp(tail_log))
    good_lower = 1.0 - tail_upper
    if good_lower <= 0.0:
        raise ArithmeticError("the selected tail-free event has no positive bound")

    allowed_count = upper_weight - lower_weight + 1
    selection_factor = allowed_count / good_lower
    selected = np.full(B + 1, -math.inf)
    selected[lower_weight : upper_weight + 1] = (
        expected[lower_weight : upper_weight + 1] + math.log(selection_factor)
    )
    return selected, tail_log, good_lower, selection_factor


def envelope(log_spectrum: np.ndarray, p: float) -> tuple[float, int]:
    best = -math.inf
    best_weight = -1
    for weight in range(1, B + 1):
        if not math.isfinite(float(log_spectrum[weight])):
            continue
        denominator = (
            log_choose(B, weight)
            + weight * math.log(p)
            + (B - weight) * math.log1p(-p)
        )
        value = float(log_spectrum[weight]) - denominator
        if value > best:
            best = value
            best_weight = weight
    return best, best_weight


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    spectrum, tail_log, good_lower, selection_factor = selected_spectrum(
        args.lower_weight, args.upper_weight
    )
    occupations = list(range(args.q_min, args.q_max + 1))
    best = {q: math.inf for q in occupations}
    witnesses: dict[int, dict[str, float | int]] = {}

    p_values = np.arange(
        args.p_min,
        args.p_max + 0.5 * args.p_step,
        args.p_step,
    )
    u_values = np.arange(
        args.grid_min,
        args.grid_max + 0.5 * args.grid_step,
        args.grid_step,
    )
    total_pairs = len(p_values) * len(u_values)
    pair_index = 0
    for p_raw in p_values:
        p = float(p_raw)
        log_envelope, envelope_weight = envelope(spectrum, p)
        for u_raw in u_values:
            pair_index += 1
            u = float(u_raw)
            surprisal = math.exp(u)
            z = math.exp(-surprisal)
            zero, active = step_matrices(z, args.memory_bits)
            candidate = (1.0 - p) * zero + p * active
            regions = uniform_candidate_regions(
                zero=zero,
                candidate=candidate,
                maximum_q=args.q_max,
            )
            correction = D * surprisal
            for q in occupations:
                inner = (
                    log_matrix_power_moment(tuple(regions[q]), B) + correction
                )
                inner = min(0.0, inner)
                value = log_choose(L, q) + q * log_envelope + inner
                if value < best[q]:
                    best[q] = value
                    witnesses[q] = {
                        "candidate_probability": p,
                        "log_surprisal": u,
                        "z": z,
                        "envelope_maximizing_weight": envelope_weight,
                        "log2_envelope": log_envelope / LOG2,
                        "inner_log2_upper": inner / LOG2,
                    }
            if args.progress_every and pair_index % args.progress_every == 0:
                print(f"pair,{pair_index},{total_pairs}", flush=True)

    rows = []
    for q in occupations:
        rows.append(
            {
                "active_outer_rows": q,
                "pointwise_log2_upper": best[q] / LOG2,
                **witnesses[q],
            }
        )
    total = float(logsumexp(np.asarray([best[q] for q in occupations])))
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    return {
        "schema": "selected-ba240-randomstepconv-g1-q2-64-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "claim_scope": {
            "outer": "one selected Golay--BA-3 constituent repeated in every row",
            "occupations": [args.q_min, args.q_max],
            "randomness": (
                "independent row-coordinate permutations, independent region "
                "permutations, and independent RandomStepConv maps"
            ),
            "log2_expected_bad_upper": total / LOG2,
            "margin_bits": -total / LOG2,
            "dominant": dominant,
        },
        "parameters": {
            "outer_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": args.memory_bits,
            "allowed_outer_weights": [args.lower_weight, args.upper_weight],
        },
        "constituent_selection": {
            "tail_word_expectation_log2": tail_log / LOG2,
            "tail_free_probability_lower": good_lower,
            "allowed_shell_count": args.upper_weight - args.lower_weight + 1,
            "common_spectrum_factor": selection_factor,
            "common_spectrum_factor_log2": math.log2(selection_factor),
            "argument": (
                "Condition on the tail-free event. The conditional mean of "
                "the sum over allowed weights of A_w/E[A_w] is at most the "
                "allowed-shell count divided by the event-probability lower "
                "bound. Hence one tail-free constituent obeys every listed "
                "pointwise spectrum bound simultaneously."
            ),
        },
        "grid": {
            "candidate_probability": [args.p_min, args.p_max, args.p_step],
            "log_surprisal": [args.grid_min, args.grid_max, args.grid_step],
            "evaluated_pairs": total_pairs,
        },
        "occupation_rows": rows,
        "limitations": [
            "The calculation uses nearest binary64 arithmetic.",
            "The witness grid is diagnostic and may be sharpened.",
            "Occupations outside the displayed range are not covered.",
            "The spectrum-selection argument proves existence but does not identify the constituent.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=30)
    parser.add_argument("--lower-weight", type=int, default=25)
    parser.add_argument("--upper-weight", type=int, default=215)
    parser.add_argument("--q-min", type=int, default=2)
    parser.add_argument("--q-max", type=int, default=64)
    parser.add_argument("--p-min", type=float, default=0.45)
    parser.add_argument("--p-max", type=float, default=0.65)
    parser.add_argument("--p-step", type=float, default=0.05)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-3.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--progress-every", type=int, default=5)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim_scope"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
