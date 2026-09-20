#!/usr/bin/env python3
"""Exact-region sparse three-band bound for the high-probability spectrum."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_single_random_constituent_highprob_bands import (
    BAND_LIMITS,
    band_majorant,
    log_multinomial,
)
from evaluate_single_random_constituent_highprob_renyi import (
    B,
    D,
    L,
    LOG2,
    MEMORY,
    spectrum_caps,
)


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    WORKSTREAM / "single_random_constituent_B256_highprob_sparse_bands.json"
)


def exact_region_matrices(
    zero: np.ndarray,
    candidates: list[np.ndarray],
    maximum_occupation: int,
    outer_rows: int,
) -> np.ndarray:
    """Average each ordered product with fixed counts of three row types."""

    side = maximum_occupation + 1
    current = np.full((side, side, side, 2, 2), -math.inf)
    current[0, 0, 0] = transfer.log_identity()
    low, central, high = np.indices((side, side, side))
    total = low + central + high
    log_zero = transfer.log_entries(zero)
    log_candidates = [transfer.log_entries(matrix) for matrix in candidates]

    for completed in range(outer_rows):
        denominator = float(completed + 1)
        updated = np.full_like(current, -math.inf)

        inactive_weight = (completed + 1 - total) / denominator
        valid = (
            (total <= completed)
            & (total <= maximum_occupation)
            & (inactive_weight > 0.0)
        )
        term = transfer.log_matmul(current, log_zero)
        updated[valid] = term[valid] + np.log(inactive_weight[valid])[:, None, None]

        for axis, log_candidate in enumerate(log_candidates):
            source = [slice(None), slice(None), slice(None)]
            target = [slice(None), slice(None), slice(None)]
            source[axis] = slice(0, maximum_occupation)
            target[axis] = slice(1, maximum_occupation + 1)
            source_tuple = tuple(source)
            target_tuple = tuple(target)
            products = transfer.log_matmul(current[source_tuple], log_candidate)
            selected_counts = (low, central, high)[axis][target_tuple]
            source_totals = total[source_tuple]
            valid = (
                (source_totals <= completed)
                & (source_totals < maximum_occupation)
            )
            terms = np.full_like(products, -math.inf)
            terms[valid] = products[valid] + np.log(
                selected_counts[valid] / denominator
            )[:, None, None]
            updated[target_tuple] = np.logaddexp(updated[target_tuple], terms)
        current = updated
    return current


def compositions(minimum: int, maximum: int) -> list[tuple[int, int, int]]:
    return [
        (low, central, occupation - low - central)
        for occupation in range(minimum, maximum + 1)
        for low in range(occupation + 1)
        for central in range(occupation - low + 1)
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-occupation", type=int, default=3)
    parser.add_argument("--maximum-occupation", type=int, default=8)
    parser.add_argument("--memory-bits", type=int, default=MEMORY)
    parser.add_argument("--block-bits", type=int, default=B)
    parser.add_argument("--dimension", type=int, default=B // 2)
    parser.add_argument("--output-bits", type=int, default=1 << 21)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument(
        "--band-limits",
        type=int,
        nargs=6,
        metavar=("LO1", "HI1", "LO2", "HI2", "LO3", "HI3"),
    )
    parser.add_argument("--per-shell-failure-exponent", type=int, default=51)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-5.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not 1 <= args.minimum_occupation <= args.maximum_occupation:
        parser.error("invalid occupation interval")
    if args.maximum_occupation > 16:
        parser.error("exact three-dimensional recurrence is limited to Q<=16")
    if args.output_bits % args.block_bits:
        parser.error("block length must divide the output length")

    block_bits = args.block_bits
    outer_rows = args.output_bits // block_bits
    distance_cutoff = (
        args.distance_numerator * args.output_bits
        + args.distance_denominator
        - 1
    ) // args.distance_denominator
    if args.band_limits is None:
        if block_bits != B:
            parser.error("nondefault block length requires --band-limits")
        band_limits = BAND_LIMITS
    else:
        band_limits = tuple(
            (args.band_limits[index], args.band_limits[index + 1])
            for index in range(0, 6, 2)
        )

    caps, event = spectrum_caps(
        math.ldexp(1.0, -args.per_shell_failure_exponent),
        block_bits=block_bits,
        dimension=args.dimension,
    )
    bands = [
        band_majorant(caps, *limits, block_bits=block_bits)
        for limits in band_limits
    ]
    candidates_to_check = compositions(
        args.minimum_occupation, args.maximum_occupation
    )
    best = {composition: math.inf for composition in candidates_to_check}
    witnesses = {composition: math.nan for composition in candidates_to_check}
    u_values = np.arange(
        args.grid_min,
        args.grid_max + args.grid_step / 2.0,
        args.grid_step,
    )

    for index, u in enumerate(u_values):
        surprisal = math.exp(float(u))
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        candidate_matrices = [
            (1.0 - float(band["value_probability"])) * zero
            + float(band["value_probability"]) * active
            for band in bands
        ]
        regions = exact_region_matrices(
            zero, candidate_matrices, args.maximum_occupation, outer_rows
        )
        for composition in candidates_to_check:
            powered = transfer.log_power(regions[composition], block_bits)
            moment = float(np.logaddexp(powered[0, 0], powered[0, 1]))
            counts = (outer_rows - sum(composition), *composition)
            outer = log_multinomial(counts) + sum(
                count * float(band["log_majorant"])
                for count, band in zip(composition, bands)
            )
            value = outer + min(0.0, moment + distance_cutoff * surprisal)
            if value < best[composition]:
                best[composition] = value
                witnesses[composition] = float(u)
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    rows = [
        {
            "occupation": sum(composition),
            "active_counts": list(composition),
            "pointwise_log2_upper": best[composition] / LOG2,
            "margin_bits": -best[composition] / LOG2,
            "log_surprisal": witnesses[composition],
        }
        for composition in candidates_to_check
    ]
    aggregate = float(logsumexp(list(best.values())))
    payload = {
        "schema": "single-random-constituent-highprob-sparse-bands-v1",
        "status": "BINARY64_EXACT_REGION_DIAGNOSTIC",
        "parameters": {
            "outer_code": (
                f"one uniform binary [{block_bits},{args.dimension}] generator"
            ),
            "outer_rows": outer_rows,
            "output_bits": args.output_bits,
            "minimum_occupation": args.minimum_occupation,
            "maximum_occupation": args.maximum_occupation,
            "distance_cutoff": distance_cutoff,
            "memory_bits": args.memory_bits,
        },
        "spectrum_event": event,
        "bands": bands,
        "claim": {
            "aggregate_log2_upper": aggregate / LOG2,
            "aggregate_margin_bits": -aggregate / LOG2,
            "worst_composition": max(
                rows, key=lambda row: float(row["pointwise_log2_upper"])
            ),
        },
        "rows": rows,
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Only the stated sparse occupations are covered.",
            "The three-band spectrum caps are pointwise upper bounds and may be loose.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
