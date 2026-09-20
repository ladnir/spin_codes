#!/usr/bin/env python3
"""Exact-shell Q=2 bound under a one-shot random-spectrum event.

The sampled constituent is fixed after setup.  On the spectrum event, its
weight-w multiplicity is at most the integer cap T_w.  The calculation uses
T_a T_b as an upper bound on the number of ordered pairs whose two images
have weights a and b.  Independent row-coordinate permutations make their
routed supports independent conditional on these weights.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients
from evaluate_single_random_constituent_highprob_renyi import (
    LOG2,
    MEMORY,
    spectrum_caps,
)
from evaluate_single_random_constituent_q2 import pair_support_coefficients


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_OUTPUT = WORKSTREAM / "single_random_constituent_B512_highprob_q2_s22.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=MEMORY)
    parser.add_argument("--block-bits", type=int, default=512)
    parser.add_argument("--dimension", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=1 << 21)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--per-shell-failure-exponent", type=int, default=51)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-7.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if args.output_bits % args.block_bits:
        parser.error("block length must divide output length")
    block_bits = args.block_bits
    outer_rows = args.output_bits // block_bits
    distance_cutoff = (
        args.distance_numerator * args.output_bits
        + args.distance_denominator
        - 1
    ) // args.distance_denominator
    caps, event = spectrum_caps(
        math.ldexp(1.0, -args.per_shell_failure_exponent),
        block_bits=block_bits,
        dimension=args.dimension,
    )
    log_caps = np.full(block_bits + 1, -math.inf)
    for weight in range(1, block_bits + 1):
        if int(caps[weight]):
            log_caps[weight] = math.log(int(caps[weight]))

    best = np.full((block_bits + 1, block_bits + 1), math.inf)
    witnesses = np.full_like(best, math.nan)
    u_values = np.arange(
        args.grid_min,
        args.grid_max + args.grid_step / 2.0,
        args.grid_step,
    )
    for index, u in enumerate(u_values):
        surprisal = math.exp(float(u))
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        regions = uniform_coefficients(
            transfer.log_entries(zero),
            transfer.log_entries(active),
            outer_rows,
            2,
        )
        coefficients = pair_support_coefficients(regions, block_bits)
        moments = np.logaddexp(coefficients[..., 0, 0], coefficients[..., 0, 1])
        values = np.minimum(0.0, moments + distance_cutoff * surprisal)
        improved = values < best
        best[improved] = values[improved]
        witnesses[improved] = u
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    row_pair_log = (
        math.lgamma(outer_rows + 1)
        - math.lgamma(3)
        - math.lgamma(outer_rows - 1)
    )
    shell_logs = row_pair_log + log_caps[:, None] + log_caps[None, :] + best
    finite = np.isfinite(shell_logs)
    aggregate = float(logsumexp(shell_logs[finite]))
    flat = shell_logs.ravel()
    top = np.argsort(flat)[-20:][::-1]
    width = block_bits + 1
    rows = [
        {
            "first_weight": int(index // width),
            "second_weight": int(index % width),
            "pointwise_log2_upper": float(flat[index] / LOG2),
            "log_surprisal": float(
                witnesses[index // width, index % width]
            ),
        }
        for index in top
        if math.isfinite(float(flat[index]))
    ]
    payload = {
        "schema": "single-random-constituent-highprob-q2-v1",
        "status": "BINARY64_DIAGNOSTIC_EXACT_SHELL_TRANSFER",
        "claim": {
            "q2_log2_upper": aggregate / LOG2,
            "q2_margin_bits": -aggregate / LOG2,
            "closes_40_bits": aggregate < -40.0 * LOG2,
            "dominant": rows[0],
        },
        "parameters": {
            "outer_code": (
                f"one uniform binary [{block_bits},{args.dimension}] generator"
            ),
            "outer_rows": outer_rows,
            "output_bits": args.output_bits,
            "distance_cutoff": distance_cutoff,
            "memory_bits": args.memory_bits,
            "occupation": 2,
        },
        "spectrum_event": event,
        "top_shell_pairs": rows,
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Only occupation two is covered.",
            "The cap product T_a T_b ignores linear-code pair correlations; it is a valid but potentially loose count bound.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
