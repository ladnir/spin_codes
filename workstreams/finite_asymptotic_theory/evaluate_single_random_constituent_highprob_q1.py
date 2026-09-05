#!/usr/bin/env python3
"""Exact-shell Q=1 bound for the high-probability [256,128] spectrum caps."""

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


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_OUTPUT = WORKSTREAM / "single_random_constituent_B256_highprob_q1.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-bits", type=int, default=MEMORY)
    parser.add_argument("--block-bits", type=int, default=256)
    parser.add_argument("--dimension", type=int, default=128)
    parser.add_argument("--output-bits", type=int, default=1 << 21)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--per-shell-failure-exponent", type=int, default=51)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--fine-min", type=float, default=-10.0)
    parser.add_argument("--fine-max", type=float, default=-7.0)
    parser.add_argument("--fine-step", type=float, default=0.05)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if args.output_bits % args.block_bits:
        parser.error("block length must divide the output length")
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
    weights = [
        weight for weight in range(1, block_bits + 1) if int(caps[weight])
    ]
    best = {weight: math.inf for weight in weights}
    witnesses = {weight: math.nan for weight in weights}
    coarse = np.arange(
        args.grid_min,
        args.grid_max + args.grid_step / 2.0,
        args.grid_step,
    )
    fine = np.arange(
        args.fine_min,
        args.fine_max + args.fine_step / 2.0,
        args.fine_step,
    )
    u_values = sorted(set(float(value) for value in coarse) | set(float(value) for value in fine))

    for index, u in enumerate(u_values):
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        region_coefficients = uniform_coefficients(
            transfer.log_entries(zero),
            transfer.log_entries(active),
            outer_rows,
            1,
        )
        coordinate_coefficients = uniform_coefficients(
            region_coefficients[0],
            region_coefficients[1],
            block_bits,
            block_bits,
        )
        moments = np.logaddexp(
            coordinate_coefficients[1:, 0, 0],
            coordinate_coefficients[1:, 0, 1],
        )
        for weight in weights:
            value = (
                math.log(outer_rows)
                + math.log(int(caps[weight]))
                + min(
                    0.0,
                    float(moments[weight - 1])
                    + distance_cutoff * surprisal,
                )
            )
            if value < best[weight]:
                best[weight] = value
                witnesses[weight] = u
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    aggregate = float(logsumexp(list(best.values())))
    rows = [
        {
            "weight": weight,
            "spectrum_cap": str(int(caps[weight])),
            "pointwise_log2_upper": best[weight] / LOG2,
            "log_surprisal": witnesses[weight],
        }
        for weight in weights
    ]
    payload = {
        "schema": "single-random-constituent-highprob-q1-v1",
        "status": "BINARY64_DIAGNOSTIC_EXACT_SHELL_TRANSFER",
        "claim": {
            "q1_log2_upper": aggregate / LOG2,
            "q1_margin_bits": -aggregate / LOG2,
            "closes_40_bits": aggregate < -40.0 * LOG2,
            "dominant": max(rows, key=lambda row: float(row["pointwise_log2_upper"])),
        },
        "parameters": {
            "outer_code": (
                f"one uniform binary [{block_bits},{args.dimension}] generator"
            ),
            "outer_rows": outer_rows,
            "output_bits": args.output_bits,
            "distance_cutoff": distance_cutoff,
            "memory_bits": args.memory_bits,
        },
        "spectrum_event": event,
        "weight_rows": rows,
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Only occupation one is covered.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
