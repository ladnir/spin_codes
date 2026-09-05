#!/usr/bin/env python3
"""All-message rank diagnostic for one repeated random outer and RandomStepConv.

Sample one uniform K-by-B binary generator and repeat it in L outer rows.
For a parent message matrix U, classify the message by r=rank(U). Choose r
rows that form an information set for the column space of U and delete every
one outside those rows. Input-one deletion can only increase the tilted
RandomStepConv moment. Each outer coordinate now contains r independent fair
bits, and the independent region permutation places them in a uniform
r-subset of the L bit positions.

The exact number of L-by-K matrices of rank r completes the all-message sum.
Nearest binary64 log arithmetic is diagnostic, not an outward certificate.
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

from analyze_riffle_striped_random_outer import (  # noqa: E402
    log_matrix_power_moment,
)
from evaluate_selected_ba_randomstepconv_g1_dense import (  # noqa: E402
    log_uniform_candidate_regions,
)
from evaluate_selected_ba_randomstepconv_g1_q2_64 import (  # noqa: E402
    B,
    D,
    K,
    L,
    LOG2,
    N,
    step_matrices,
)


DEFAULT_OUTPUT = (
    WORKSTREAM / "repeated_random240_randomstepconv_g1_s19_rank_d11.json"
)


def log_two_power_difference(high: int, low: int) -> float:
    if not 0 <= low < high:
        raise ValueError("expected 0 <= low < high")
    return high * LOG2 + math.log1p(-math.ldexp(1.0, low - high))


def log_rank_matrix_count(rows: int, columns: int, rank: int) -> float:
    if not 0 <= rank <= min(rows, columns):
        return -math.inf
    if rank == 0:
        return 0.0
    result = 0.0
    for index in range(rank):
        result += log_two_power_difference(rows, index)
        result += log_two_power_difference(columns, index)
        result -= log_two_power_difference(rank, index)
    return result


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    log_counts = np.asarray(
        [log_rank_matrix_count(L, K, rank) for rank in range(K + 1)]
    )
    if abs(log_counts[K] / LOG2 - L * K) > 1.0:
        raise ArithmeticError("rank-count scale check failed")

    best = np.full(K + 1, math.inf)
    witnesses: list[dict[str, float] | None] = [None] * (K + 1)
    grid_count = int(
        math.floor((args.grid_max - args.grid_min) / args.grid_step + 0.5)
    ) + 1
    for index in range(grid_count):
        u = args.grid_min + index * args.grid_step
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = step_matrices(z, args.memory_bits)
        fair_candidate = 0.5 * zero + 0.5 * active
        regions = log_uniform_candidate_regions(
            zero=zero,
            candidate=fair_candidate,
            maximum_q=K,
        )
        correction = D * surprisal
        for rank in range(1, K + 1):
            entries = regions[rank].reshape(4)
            scale = float(np.max(entries))
            normalized = tuple(
                math.exp(float(value - scale)) for value in entries
            )
            inner = (
                log_matrix_power_moment(normalized, B)
                + B * scale
                + correction
            )
            inner = min(0.0, inner)
            value = float(log_counts[rank]) + inner
            if value < best[rank]:
                best[rank] = value
                witnesses[rank] = {
                    "log_surprisal": u,
                    "surprisal": surprisal,
                    "z": z,
                    "message_count_log2": float(log_counts[rank]) / LOG2,
                    "inner_log2_upper": inner / LOG2,
                }
        print(f"tilt,{index + 1},{grid_count},u,{u:.6f}", flush=True)

    rows = [
        {
            "message_matrix_rank": rank,
            "pointwise_log2_upper": float(best[rank]) / LOG2,
            **(witnesses[rank] or {}),
        }
        for rank in range(1, K + 1)
    ]
    aggregate = float(logsumexp(best[1:]))
    dominant = max(rows, key=lambda row: float(row["pointwise_log2_upper"]))
    return {
        "schema": "repeated-random240-randomstepconv-g1-rank-v1",
        "status": "BINARY64_DIAGNOSTIC_EXACT_RANK_COUNT",
        "claim": {
            "log2_expected_bad_upper": aggregate / LOG2,
            "margin_bits": -aggregate / LOG2,
            "closes_40_bits": aggregate / LOG2 < -40.0,
            "all_nonzero_parent_messages_covered": True,
            "dominant": dominant,
        },
        "construction": {
            "outer": (
                "one uniform random binary K-by-B generator, sampled once "
                "and repeated in all L rows"
            ),
            "routing": (
                "transpose the B outer coordinates into B length-L regions "
                "and independently permute each region"
            ),
            "inner": (
                "independent uniform linear (M+1)-by-(M+1) maps at every "
                "bit position, sampled once and shared by all messages"
            ),
            "row_local_outer_coordinate_permutations": "omitted",
        },
        "parameters": {
            "outer_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "output_bits": N,
            "parent_dimension": K * L,
            "target_dimension_after_zero_shortening": 1 << 20,
            "distance_cutoff": D,
            "memory_bits": args.memory_bits,
        },
        "probability_space": (
            "the one shared random outer generator, the B independent region "
            "permutations, and the N independent RandomStepConv maps"
        ),
        "bound": {
            "rank_count": (
                "prod_{i=0}^{r-1} (2^L-2^i)(2^K-2^i)/(2^r-2^i)"
            ),
            "information_set": (
                "for each fixed rank-r parent message matrix, choose r row "
                "coordinates on which its column space projects bijectively"
            ),
            "monotonicity": (
                "delete all input ones outside the information-set rows"
            ),
            "resulting_region_law": (
                "r independent fair candidate bits placed in a uniform "
                "r-subset of the L region positions"
            ),
        },
        "grid": {
            "minimum": args.grid_min,
            "maximum": args.grid_max,
            "step": args.grid_step,
            "count": grid_count,
        },
        "rank_rows": rows,
        "limitations": [
            "The logarithmic calculation uses nearest binary64 arithmetic.",
            "The result is not yet an outward-rounded certificate.",
            "The outer is a random constituent, not Golay--BA-3.",
            "The result proves no implementation-cost bound.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=19)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = evaluate(args)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
