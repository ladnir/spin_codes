#!/usr/bin/env python3
"""Sparse transfer with two defect-severity categories and one central band.

The input cap receipt is conditional on a shell-variance hypothesis.  This
program splits the symmetric defect tail into severe and moderate parts so
that rare low weights do not inherit the change-of-measure factor attained
near weight 79.  Row categories remain shared across all transposed regions.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_single_random_constituent_highprob_bands import (
    band_majorant,
    log_multinomial,
)
from evaluate_single_random_constituent_highprob_renyi import LOG2


WORKSTREAM = Path(__file__).resolve().parent
CAPS = WORKSTREAM / "ebch128x4_ba_variance_cap_requirements.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba_variance_proxy_split_sparse.json"
B = 512
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000


def exact_shared_regions(
    zero: np.ndarray,
    categories: list[np.ndarray],
    maximum_occupation: int,
) -> np.ndarray:
    """Average ordered products with one shared count per active category."""

    side = maximum_occupation + 1
    shape = (side,) * len(categories)
    current = np.full(shape + (2, 2), -math.inf)
    current[(0,) * len(categories)] = transfer.log_identity()
    indices = np.indices(shape)
    total = np.sum(indices, axis=0)
    log_zero = transfer.log_entries(zero)
    log_categories = [transfer.log_entries(category) for category in categories]

    for completed in range(L):
        denominator = float(completed + 1)
        updated = np.full_like(current, -math.inf)
        inactive_weight = (completed + 1 - total) / denominator
        valid = (total <= completed) & (total <= maximum_occupation) & (inactive_weight > 0)
        inactive_term = transfer.log_matmul(current, log_zero)
        updated[valid] = inactive_term[valid] + np.log(inactive_weight[valid])[:, None, None]

        for axis, candidate in enumerate(log_categories):
            source = [slice(None)] * len(categories)
            target = [slice(None)] * len(categories)
            source[axis] = slice(0, maximum_occupation)
            target[axis] = slice(1, maximum_occupation + 1)
            source_tuple = tuple(source)
            target_tuple = tuple(target)
            products = transfer.log_matmul(current[source_tuple], candidate)
            selected = indices[axis][target_tuple]
            source_totals = total[source_tuple]
            valid = (source_totals <= completed) & (source_totals < maximum_occupation)
            terms = np.full_like(products, -math.inf)
            terms[valid] = products[valid] + np.log(selected[valid] / denominator)[:, None, None]
            updated[target_tuple] = np.logaddexp(updated[target_tuple], terms)
        current = updated
    return current


def merged_band(
    caps: np.ndarray,
    low_limits: tuple[int, int],
    high_limits: tuple[int, int],
) -> dict[str, object]:
    low = band_majorant(caps, *low_limits, block_bits=B)
    high = band_majorant(caps, *high_limits, block_bits=B)
    probability = min(float(low["value_probability"]), 1.0 - float(high["value_probability"]))
    log_majorant = math.log(2.0) + max(float(low["log_majorant"]), float(high["log_majorant"]))
    return {
        "low": low,
        "high": high,
        "value_probability": probability,
        "log_majorant": log_majorant,
        "log2_majorant": log_majorant / LOG2,
    }


def compositions(categories: int, total: int):
    if categories == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for suffix in compositions(categories - 1, total - first):
            yield (first,) + suffix


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accumulator-stages", type=int, default=3)
    parser.add_argument("--variance-inflation-bits", type=float, default=9.0)
    parser.add_argument("--severe-upper", type=int, default=63)
    parser.add_argument("--minimum-occupation", type=int, default=3)
    parser.add_argument("--maximum-occupation", type=int, default=8)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-3.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if not 42 <= args.severe_upper < 79:
        parser.error("severe upper weight must lie from 42 through 78")

    receipt = json.loads(CAPS.read_text(encoding="utf-8"))
    source = next(
        row
        for row in receipt["rows"]
        if int(row["accumulator_stages"]) == args.accumulator_stages
        and float(row["variance_inflation_bits"]) == args.variance_inflation_bits
    )
    caps = np.asarray([int(value) for value in source["caps"]], dtype=object)
    severe = merged_band(
        caps,
        (42, args.severe_upper),
        (B - args.severe_upper, 470),
    )
    moderate = merged_band(
        caps,
        (args.severe_upper + 1, 79),
        (433, B - args.severe_upper - 1),
    )
    central = band_majorant(caps, 80, 432, block_bits=B)
    probabilities = [
        float(severe["value_probability"]),
        float(moderate["value_probability"]),
        float(central["value_probability"]),
    ]
    majorants = [
        float(severe["log_majorant"]),
        float(moderate["log_majorant"]),
        float(central["log_majorant"]),
    ]

    targets = [
        composition
        for occupation in range(args.minimum_occupation, args.maximum_occupation + 1)
        for composition in compositions(3, occupation)
    ]
    best = {target: math.inf for target in targets}
    witnesses = {target: math.nan for target in targets}
    u_values = np.arange(args.grid_min, args.grid_max + args.grid_step / 2, args.grid_step)
    for tilt_index, u_raw in enumerate(u_values):
        u = float(u_raw)
        surprisal = math.exp(u)
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        matrices = [(1.0 - probability) * zero + probability * active for probability in probabilities]
        regions = exact_shared_regions(zero, matrices, args.maximum_occupation)
        for target in targets:
            powered = transfer.log_power(regions[target], B)
            moment = float(np.logaddexp(powered[0, 0], powered[0, 1]))
            outer = log_multinomial((L - sum(target),) + target) + sum(
                count * majorant for count, majorant in zip(target, majorants)
            )
            value = outer + min(0.0, moment + D * surprisal)
            if value < best[target]:
                best[target] = value
                witnesses[target] = u
        print(f"tilt,{tilt_index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    rows = [
        {
            "occupation": sum(target),
            "severe_rows": target[0],
            "moderate_rows": target[1],
            "central_rows": target[2],
            "pointwise_log2_upper": best[target] / LOG2,
            "margin_bits": -best[target] / LOG2,
            "log_surprisal": witnesses[target],
        }
        for target in targets
    ]
    aggregate = float(logsumexp(list(best.values())))
    result = {
        "schema": "ebch128x4-ba-variance-proxy-split-sparse-v1",
        "status": "CONDITIONAL_BINARY64_EXACT_SHARED_CATEGORY_DIAGNOSTIC",
        "parameters": {
            "accumulator_stages": args.accumulator_stages,
            "variance_inflation_bits": args.variance_inflation_bits,
            "variance_hypothesis": source["variance_hypothesis"],
            "memory_bits": args.memory_bits,
            "occupations": [args.minimum_occupation, args.maximum_occupation],
            "distance_cutoff": D,
            "severe_upper": args.severe_upper,
        },
        "cap_event_failure_log2_upper": source["event_failure_log2_upper"],
        "bands": {"severe": severe, "moderate": moderate, "central": central},
        "claim": {
            "aggregate_log2_upper": aggregate / LOG2,
            "aggregate_margin_bits": -aggregate / LOG2,
            "worst": max(rows, key=lambda row: float(row["pointwise_log2_upper"])),
        },
        "rows": rows,
        "limitations": [
            "The shell-variance hypothesis is not proved.",
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Only the stated sparse occupation interval is covered.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
