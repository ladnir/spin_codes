#!/usr/bin/env python3
"""Fast one-active-block contour for TransposePacketShuffle-RandomStepConv.

The one-active-block occupation is exact under the random-outer relaxation
and has dominated the certified g=4 point.  This scanner locates candidate
outer-size and inner-memory tradeoffs.  A positive sparse margin is necessary,
not sufficient; candidate boundary cells still require the full occupation
checker.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from analyze_riffle_striped_random_outer import (
    LOG2,
    log_matrix_power_moment,
    log_two_power_minus_one,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_transpose_packetshuffle_randomstepconv/"
    "receipts/g4_sparse_tradeoff.json"
)


def transition_matrices(packet_bits: int, sigma: int, z: float):
    state_zero = math.ldexp(1.0, -sigma)
    output_moment = math.ldexp((1.0 + z) ** packet_bits, -packet_bits)
    off = state_zero * output_moment
    live = (1.0 - state_zero) * output_moment
    zero = np.asarray(((1.0, 0.0), (off, live)))
    nonzero = np.asarray(((off, live), (off, live)))
    return zero, 0.5 * zero + 0.5 * nonzero


def zero_powers(zero: np.ndarray, count: int) -> np.ndarray:
    """Return M0^j for j=0,...,count-1 using the triangular form."""
    powers = np.zeros((count, 2, 2), dtype=np.float64)
    powers[:, 0, 0] = 1.0
    exponents = np.arange(count, dtype=np.float64)
    diagonal = np.power(zero[1, 1], exponents)
    powers[:, 1, 1] = diagonal
    denominator = 1.0 - zero[1, 1]
    powers[:, 1, 0] = zero[1, 0] * (1.0 - diagonal) / denominator
    return powers


def one_marked_row(zero: np.ndarray, marked: np.ndarray, positions: int):
    powers = zero_powers(zero, positions)
    products = np.matmul(np.matmul(powers, marked), powers[::-1])
    return np.mean(products, axis=0)


def sparse_point(
    *,
    message_bits: int,
    outer_bits: int,
    packet_bits: int,
    sigma: int,
    relative_distance: float,
) -> dict[str, float | int]:
    outer_dimension = outer_bits // 2
    blocks = message_bits // outer_dimension
    groups = blocks // packet_bits
    output_bits = outer_bits * groups * packet_bits
    distance = math.floor(relative_distance * output_bits)

    def objective(log_surprisal: float) -> float:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero, marked = transition_matrices(packet_bits, sigma, z)
        row = one_marked_row(zero, marked, groups)
        moment = log_matrix_power_moment(tuple(row.reshape(4)), outer_bits)
        return moment + distance * surprisal

    result = minimize_scalar(
        objective,
        method="bounded",
        bounds=(-16.0, 3.0),
        options={"xatol": 1e-10, "maxiter": 300},
    )
    conditioning_penalty = -math.log1p(-math.ldexp(1.0, -outer_bits))
    outer_log2 = (
        math.log(blocks)
        + log_two_power_minus_one(outer_dimension)
        + conditioning_penalty
    ) / LOG2
    inner_log2 = min(0.0, float(result.fun) / LOG2)
    return {
        "outer_bits": outer_bits,
        "outer_blocks": blocks,
        "fixed_block_groups": groups,
        "sigma": sigma,
        "distance": distance,
        "outer_log2": outer_log2,
        "inner_log2_upper": inner_log2,
        "pointwise_log2_upper": outer_log2 + inner_log2,
        "sparse_lambda_bits": -(outer_log2 + inner_log2),
        "optimal_log_surprisal": float(result.x),
        "optimizer_success": bool(result.success),
    }


def parse_ints(text: str) -> list[int]:
    return [int(item) for item in text.split(",") if item]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--g", type=int, default=4)
    parser.add_argument("--outer-bits", default="128,256,512,1024")
    parser.add_argument("--sigmas", default="8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32")
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = []
    for outer_bits in parse_ints(args.outer_bits):
        for sigma in parse_ints(args.sigmas):
            row = sparse_point(
                message_bits=args.message_bits,
                outer_bits=outer_bits,
                packet_bits=args.g,
                sigma=sigma,
                relative_distance=args.relative_distance,
            )
            rows.append(row)
            print(
                f"B,{outer_bits},sigma,{sigma},"
                f"sparse_lambda,{row['sparse_lambda_bits']:.6f}",
                flush=True,
            )
    payload = {
        "schema": "riffle-transpose-packetshuffle-sparse-curve-v1",
        "message_bits": args.message_bits,
        "packet_bits": args.g,
        "relative_distance": args.relative_distance,
        "scope": "Exact one-active-block diagnostic; full occupation certification is required at selected cells.",
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()

