#!/usr/bin/env python3
"""Occupation-one diagnostic for one repeated random constituent.

For each requested even block length B, sample one uniform linear injection
from F_2^(B/2) to F_2^B and reuse it in every outer row.  The calculation
averages over that one constituent, the row and region permutations, and the
RandomStepConv maps.  It covers occupation one exactly in the proof model.
Nearest binary64 arithmetic makes the result diagnostic.
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


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_OUTPUT = WORKSTREAM / "single_random_constituent_q1_d109.json"
LOG2 = math.log(2.0)


def log_two_power_minus_one(exponent: int) -> float:
    return exponent * LOG2 + math.log1p(-math.ldexp(1.0, -exponent))


def log_choose(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def evaluate_block(
    *,
    block_bits: int,
    output_bits: int,
    distance_cutoff: int,
    memory_bits: int,
    u_values: list[float],
) -> dict[str, object]:
    if block_bits <= 0 or block_bits % 2:
        raise ValueError("block length must be positive and even")
    if output_bits % block_bits:
        raise ValueError("block length must divide the output length")
    dimension = block_bits // 2
    outer_rows = output_bits // block_bits
    expected_shell_logs = np.asarray(
        [
            log_two_power_minus_one(dimension)
            + log_choose(block_bits, weight)
            - log_two_power_minus_one(block_bits)
            for weight in range(1, block_bits + 1)
        ]
    )
    best = np.full(block_bits, math.inf)
    witnesses = np.full(block_bits, math.nan)

    for index, u in enumerate(u_values):
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, memory_bits)
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
        values = (
            math.log(outer_rows)
            + expected_shell_logs
            + np.minimum(0.0, moments + distance_cutoff * surprisal)
        )
        improved = values < best
        best[improved] = values[improved]
        witnesses[improved] = u
        print(
            f"B,{block_bits},tilt,{index + 1},{len(u_values)},u,{u:.6f}",
            flush=True,
        )

    aggregate = float(logsumexp(best))
    dominant = int(np.argmax(best))
    return {
        "block_bits": block_bits,
        "dimension": dimension,
        "outer_rows": outer_rows,
        "message_bits": dimension * outer_rows,
        "output_bits": output_bits,
        "distance_cutoff": distance_cutoff,
        "memory_bits": memory_bits,
        "claim": {
            "occupation_one_log2_upper": aggregate / LOG2,
            "occupation_one_margin_bits": -aggregate / LOG2,
            "closes_40_bits": aggregate < -40.0 * LOG2,
            "dominant_weight": dominant + 1,
            "dominant_pointwise_log2_upper": float(best[dominant]) / LOG2,
            "dominant_log_surprisal": float(witnesses[dominant]),
        },
        "weight_rows": [
            {
                "weight": weight,
                "expected_multiplicity_log2": float(
                    expected_shell_logs[weight - 1] / LOG2
                ),
                "pointwise_log2_upper": float(best[weight - 1] / LOG2),
                "log_surprisal": float(witnesses[weight - 1]),
            }
            for weight in range(1, block_bits + 1)
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-bits", type=int, nargs="+", default=[128, 256, 512])
    parser.add_argument("--output-bits", type=int, default=1 << 21)
    parser.add_argument("--memory-bits", type=int, default=30)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--fine-min", type=float, default=-10.0)
    parser.add_argument("--fine-max", type=float, default=-7.0)
    parser.add_argument("--fine-step", type=float, default=0.05)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not 0 < args.distance_numerator < args.distance_denominator:
        raise ValueError("distance target must lie strictly between zero and one")
    distance_cutoff = (
        args.distance_numerator * args.output_bits
        + args.distance_denominator
        - 1
    ) // args.distance_denominator
    u_values = sorted(
        set(
            float(value)
            for value in np.arange(
                args.grid_min,
                args.grid_max + args.grid_step / 2.0,
                args.grid_step,
            )
        )
        | set(
            float(value)
            for value in np.arange(
                args.fine_min,
                args.fine_max + args.fine_step / 2.0,
                args.fine_step,
            )
        )
    )
    blocks = [
        evaluate_block(
            block_bits=block_bits,
            output_bits=args.output_bits,
            distance_cutoff=distance_cutoff,
            memory_bits=args.memory_bits,
            u_values=u_values,
        )
        for block_bits in args.block_bits
    ]
    payload = {
        "schema": "single-random-constituent-randomstepconv-q1-v1",
        "status": "BINARY64_DIAGNOSTIC_EXACT_OCCUPATION_ONE",
        "probability_space": (
            "for each block size, sample one uniform rate-half linear injection "
            "and reuse it in every outer row; independently sample every row "
            "coordinate permutation, region permutation, and RandomStepConv map"
        ),
        "outer_average": (
            "for occupation one, each fixed nonzero local message maps to a "
            "uniform nonzero block word; linearity of expectation is exact and "
            "does not resample the constituent across outer rows"
        ),
        "parameters": {
            "output_bits": args.output_bits,
            "message_bits": args.output_bits // 2,
            "distance_cutoff": distance_cutoff,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "memory_bits": args.memory_bits,
            "block_bits": args.block_bits,
        },
        "blocks": blocks,
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Only occupation one is covered.",
            "Occupations two and above require the joint law induced by reusing one constituent.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({block["block_bits"]: block["claim"] for block in blocks}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
