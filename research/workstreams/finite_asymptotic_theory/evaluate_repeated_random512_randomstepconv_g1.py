#!/usr/bin/env python3
"""All-occupation diagnostic for one repeated random [512,256] outer.

Select one full-rank random binary [512,256] constituent satisfying
A_w <= 4 E[A_w] at every nonzero weight. Repeat that one code in 4140 rows,
shorten 44 whole information rows to leave 4096 active message rows, sample
independent row-coordinate and region permutations, and apply the bitwise
RandomStepConv inner. A Bernoulli-half spectrum envelope and an explicit
global-zero exclusion cover all occupations in one finite sum.

Nearest binary64 log arithmetic is diagnostic, not an outward certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
B = 512
K = 256
L = 4140
ACTIVE_ROWS = 4096
N = B * L
D = (11 * N) // 100
MEMORY = 30
SPECTRUM_FACTOR = 4
LOG2 = math.log(2.0)
DEFAULT_OUTPUT = (
    WORKSTREAM / "repeated_random512_randomstepconv_g1_s30_allq_d11.json"
)


def log_two_power_minus_one(exponent: int) -> float:
    return exponent * LOG2 + math.log1p(-math.ldexp(1.0, -exponent))


def log_choose(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def log_matrix_entries(matrix: tuple[float, ...]) -> np.ndarray:
    return np.asarray(
        [math.log(value) if value > 0.0 else -math.inf for value in matrix],
        dtype=np.float64,
    ).reshape(2, 2)


def log_matmul_batch(coefficients: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.empty_like(coefficients)
    result[:, 0, 0] = np.logaddexp(
        coefficients[:, 0, 0] + right[0, 0],
        coefficients[:, 0, 1] + right[1, 0],
    )
    result[:, 0, 1] = np.logaddexp(
        coefficients[:, 0, 0] + right[0, 1],
        coefficients[:, 0, 1] + right[1, 1],
    )
    result[:, 1, 0] = np.logaddexp(
        coefficients[:, 1, 0] + right[0, 0],
        coefficients[:, 1, 1] + right[1, 0],
    )
    result[:, 1, 1] = np.logaddexp(
        coefficients[:, 1, 0] + right[0, 1],
        coefficients[:, 1, 1] + right[1, 1],
    )
    return result


def step_matrices(z: float, memory_bits: int) -> tuple[np.ndarray, np.ndarray]:
    q = math.ldexp(1.0, -memory_bits)
    b = (1.0 + z) / 2.0
    terminate = q * b
    survive = (1.0 - q) * b
    return (
        np.asarray((1.0, 0.0, terminate, survive)),
        np.asarray((terminate, survive, terminate, survive)),
    )


def log_uniform_candidate_regions(
    *, zero: np.ndarray, candidate: np.ndarray
) -> np.ndarray:
    log_zero = log_matrix_entries(tuple(zero))
    log_candidate = log_matrix_entries(tuple(candidate))
    current = np.full((ACTIVE_ROWS + 1, 2, 2), -math.inf)
    current[0, 0, 0] = 0.0
    current[0, 1, 1] = 0.0
    for completed in range(L):
        next_maximum = min(completed + 1, ACTIVE_ROWS)
        updated = np.full_like(current, -math.inf)

        zero_products = log_matmul_batch(
            current[: next_maximum + 1], log_zero
        )
        degrees = np.arange(next_maximum + 1, dtype=np.float64)
        weights = ((completed + 1.0) - degrees) / (completed + 1.0)
        positive = weights > 0.0
        updated[: next_maximum + 1][positive] = (
            zero_products[positive] + np.log(weights[positive])[:, None, None]
        )

        if next_maximum > 0:
            candidate_products = log_matmul_batch(
                current[:next_maximum], log_candidate
            )
            selected = np.arange(1, next_maximum + 1, dtype=np.float64)
            candidate_terms = candidate_products + np.log(
                selected / (completed + 1.0)
            )[:, None, None]
            updated[1 : next_maximum + 1] = np.logaddexp(
                updated[1 : next_maximum + 1], candidate_terms
            )
        current = updated
    return current


def logadd(*values: float) -> float:
    maximum = max(values)
    if maximum == -math.inf:
        return maximum
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def logdiffexp(high: float, low: float) -> float:
    if not high > low:
        raise ArithmeticError("positive activation mass was lost")
    return high + math.log(-math.expm1(low - high))


def nonzero_region_power(log_region: np.ndarray, occupation: int) -> float:
    log_zero_region = -occupation * LOG2
    activate_zero = logdiffexp(float(log_region[0, 0]), log_zero_region)
    activate_live = float(log_region[0, 1])
    inactive = 0.0
    active_zero = -math.inf
    active_live = -math.inf
    for _ in range(B):
        next_zero = logadd(
            inactive + activate_zero,
            active_zero + float(log_region[0, 0]),
            active_live + float(log_region[1, 0]),
        )
        next_live = logadd(
            inactive + activate_live,
            active_zero + float(log_region[0, 1]),
            active_live + float(log_region[1, 1]),
        )
        inactive += log_zero_region
        active_zero, active_live = next_zero, next_live
    return logadd(active_zero, active_live)


def spectrum_event_failure_upper() -> dict[str, float]:
    message_count = (1 << K) - 1
    tail = 0.0
    chebyshev = 0.0
    for weight in range(1, B + 1):
        mean = message_count * math.comb(B, weight) / (1 << B)
        if mean < 1.0 / SPECTRUM_FACTOR:
            tail += mean
        else:
            # Pairwise independence gives Var(A_w)<=E[A_w].
            deviation = (SPECTRUM_FACTOR - 1.0) * mean
            chebyshev += mean / (deviation * deviation)
    rank_failure = sum(math.ldexp(1.0, index - B) for index in range(K))
    total = tail + chebyshev + rank_failure
    return {
        "low_mean_markov_union_upper": tail,
        "central_chebyshev_union_upper": chebyshev,
        "rank_failure_union_upper": rank_failure,
        "total_failure_upper": total,
        "success_probability_lower": 1.0 - total,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    best = np.full(ACTIVE_ROWS + 1, math.inf)
    witnesses: list[dict[str, float] | None] = [None] * (ACTIVE_ROWS + 1)
    per_row_envelope = (
        log_two_power_minus_one(K) + math.log(SPECTRUM_FACTOR)
    )
    u_values = list(np.arange(args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step))
    u_values.extend(
        np.arange(
            args.refine_min,
            args.refine_max + args.refine_step / 2,
            args.refine_step,
        )
    )
    u_values = sorted(set(float(value) for value in u_values))

    for index, u in enumerate(u_values):
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = step_matrices(z, args.memory_bits)
        regions = log_uniform_candidate_regions(
            zero=zero, candidate=0.5 * zero + 0.5 * active
        )
        correction = D * surprisal
        for occupation in range(1, ACTIVE_ROWS + 1):
            inner = nonzero_region_power(regions[occupation], occupation)
            inner = min(0.0, inner + correction)
            value = (
                log_choose(ACTIVE_ROWS, occupation)
                + occupation * per_row_envelope
                + inner
            )
            if value < best[occupation]:
                best[occupation] = value
                witnesses[occupation] = {
                    "log_surprisal": u,
                    "surprisal": surprisal,
                    "z": z,
                    "inner_log2_upper": inner / LOG2,
                }
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    rows = [
        {
            "active_outer_rows": occupation,
            "pointwise_log2_upper": float(best[occupation]) / LOG2,
            **(witnesses[occupation] or {}),
        }
        for occupation in range(1, ACTIVE_ROWS + 1)
    ]
    aggregate = float(logsumexp(best[1:]))
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    setup = spectrum_event_failure_upper()
    return {
        "schema": "repeated-random512-randomstepconv-g1-allq-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "claim": {
            "log2_expected_bad_upper": aggregate / LOG2,
            "margin_bits": -aggregate / LOG2,
            "closes_40_bits": aggregate / LOG2 < -40.0,
            "occupations": [1, ACTIVE_ROWS],
            "dominant": dominant,
        },
        "parameters": {
            "outer_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "active_outer_rows": ACTIVE_ROWS,
            "shortened_outer_rows": L - ACTIVE_ROWS,
            "parent_message_bits": K * L,
            "message_bits": K * ACTIVE_ROWS,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": args.memory_bits,
            "spectrum_factor": SPECTRUM_FACTOR,
        },
        "constituent_event": {
            "event": "full rank and A_w <= 4 E[A_w] for every 1<=w<=512",
            **setup,
        },
        "probability_space": (
            "first select one random outer constituent satisfying the stated "
            "event; then sample independent row-coordinate permutations, "
            "independent region permutations, and all RandomStepConv maps"
        ),
        "bound": {
            "per_row_envelope_log2": per_row_envelope / LOG2,
            "global_zero_exclusion": (
                "the three-state region recurrence excludes the globally "
                "all-zero Bernoulli reference input"
            ),
            "zero_reference_rows": (
                "reference tuples with some zero rows remain included and "
                "therefore only enlarge the upper bound"
            ),
        },
        "tilts": u_values,
        "occupation_rows": rows,
        "limitations": [
            "The distance sum uses nearest binary64 log arithmetic.",
            "The constituent-event probability uses elementary union bounds.",
            "The result is not yet an outward-rounded certificate.",
            "Efficiently testing the constituent event is not supplied.",
            "The model is a proof baseline and has no implementation-cost claim.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=MEMORY)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--refine-min", type=float, default=0.65)
    parser.add_argument("--refine-max", type=float, default=0.80)
    parser.add_argument("--refine-step", type=float, default=0.025)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(json.dumps(payload["constituent_event"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
