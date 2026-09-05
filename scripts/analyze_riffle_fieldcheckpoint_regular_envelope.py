#!/usr/bin/env python3
"""Random-support envelope for regular FieldCheckpointAccumulate outer words."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_fieldcheckpoint_accumulate_oneblock import (
    occupancy_epoch_matrices,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import LOG2, log_choose, log_two_power_minus_one


DEFAULT_OUTPUT = Path(
    "constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/"
    "receipts/goal06_regular_envelope_a64.json"
)


def spectrum_density_envelope_log(
    outer_bits: int, dimension: int, spectrum: np.ndarray
) -> tuple[float, int]:
    """Return log eta for mu(S) <= eta (2^k-1)/2^B."""
    random_log_density = log_two_power_minus_one(dimension) - outer_bits * LOG2
    best = -math.inf
    best_weight = -1
    for weight in range(1, outer_bits):
        if not math.isfinite(float(spectrum[weight])):
            continue
        log_density = float(spectrum[weight]) - log_choose(outer_bits, weight)
        ratio = log_density - random_log_density
        if ratio > best:
            best = ratio
            best_weight = weight
    return best, best_weight


def candidate_epoch_matrices(
    state_bits: int, epoch_bits: int, maximum_candidates: int, z: float
) -> list[np.ndarray]:
    """Average an epoch with c uniform candidate positions and fair bits."""
    impulses = occupancy_epoch_matrices(state_bits, epoch_bits, z)
    candidates = []
    for count in range(maximum_candidates + 1):
        matrix = np.zeros((2, 2), dtype=np.float64)
        scale = math.ldexp(1.0, -count)
        for active in range(count + 1):
            matrix += math.comb(count, active) * scale * impulses[active]
        candidates.append(matrix)
    return candidates


def regular_region_matrices(
    state_bits: int,
    epoch_bits: int,
    epochs_per_region: int,
    maximum_candidates: int,
    z: float,
) -> list[np.ndarray]:
    """Average one region for each regular active-block count."""
    candidate_epoch = candidate_epoch_matrices(
        state_bits, epoch_bits, maximum_candidates, z
    )
    weighted_epoch = [
        math.comb(epoch_bits, count) * candidate_epoch[count]
        for count in range(maximum_candidates + 1)
    ]
    coefficients = [np.zeros((2, 2)) for _ in range(maximum_candidates + 1)]
    coefficients[0] = np.eye(2)
    for _ in range(epochs_per_region):
        updated = [np.zeros((2, 2)) for _ in range(maximum_candidates + 1)]
        for total in range(maximum_candidates + 1):
            for current in range(total + 1):
                updated[total] += coefficients[total - current] @ weighted_epoch[current]
        coefficients = updated
    region_bits = epoch_bits * epochs_per_region
    return [
        coefficients[count] / math.comb(region_bits, count)
        for count in range(maximum_candidates + 1)
    ]


def self_test() -> dict[str, float]:
    regions = regular_region_matrices(3, 6, 2, 6, 1.0)
    stochastic_error = max(
        float(np.max(np.abs(matrix.sum(axis=1) - 1.0)))
        for matrix in regions
    )
    if stochastic_error > 4e-13:
        raise AssertionError("regular-envelope region transfer is not stochastic")
    return {"maximum_small_region_stochastic_error": stochastic_error}


def log_row_power_moment(matrix: np.ndarray, power: int) -> float:
    """Return log(e_0^T matrix^power 1) with per-factor rescaling."""
    row = np.asarray((1.0, 0.0), dtype=np.float64)
    accumulated = 0.0
    for _ in range(power):
        row = row @ matrix
        scale = float(np.max(row))
        if scale <= 0.0:
            return -math.inf
        row /= scale
        accumulated += math.log(scale)
    return accumulated + math.log(float(row.sum()))


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    per_block_log_mass = log_two_power_minus_one(dimension) + log_eta

    best = np.full(args.maximum_active_blocks + 1, math.inf)
    best_tilt = np.full(args.maximum_active_blocks + 1, math.nan)
    grid_count = int(
        math.floor((args.grid_max - args.grid_min) / args.grid_step + 0.5)
    ) + 1
    for grid_index in range(grid_count):
        log_surprisal = args.grid_min + grid_index * args.grid_step
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        regions = regular_region_matrices(
            args.state_bits,
            args.epoch_bits,
            args.epochs_per_region,
            args.maximum_active_blocks,
            z,
        )
        for active_blocks in range(1, args.maximum_active_blocks + 1):
            log_moment = log_row_power_moment(
                regions[active_blocks], args.outer_bits
            )
            candidate = log_moment + distance * surprisal
            if candidate < best[active_blocks]:
                best[active_blocks] = candidate
                best_tilt[active_blocks] = log_surprisal
        print(
            f"grid,{grid_index + 1},{grid_count},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    contributions = []
    for active_blocks in range(1, args.maximum_active_blocks + 1):
        outer_log = (
            log_choose(outer_blocks, active_blocks)
            + active_blocks * per_block_log_mass
        )
        inner_log = min(0.0, float(best[active_blocks]))
        contribution = outer_log + inner_log
        contributions.append(contribution)
        rows.append(
            {
                "active_regular_outer_blocks": active_blocks,
                "best_log_surprisal": float(best_tilt[active_blocks]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
            }
        )
    partial_log = float(logsumexp(np.asarray(contributions)))
    dominant = sorted(
        rows, key=lambda row: float(row["pointwise_log2_upper"]), reverse=True
    )[:20]
    return {
        "schema": "riffle-fieldcheckpoint-regular-envelope-v1",
        "candidate": "Riffle FieldCheckpointAccumulate t=32 s=64 K=32",
        "envelope": {
            "definition": (
                "For every regular support S, expected permuted support "
                "multiplicity mu(S) is at most eta*(2^128-1)/2^256."
            ),
            "eta": math.exp(log_eta),
            "log2_eta": log_eta / LOG2,
            "maximizing_weight": envelope_weight,
            "special_support_excluded": "the unique all-one outer word",
        },
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "state_bits": args.state_bits,
            "epoch_bits": args.epoch_bits,
            "epochs_per_region": args.epochs_per_region,
            "maximum_active_blocks": args.maximum_active_blocks,
        },
        "self_test": self_test(),
        "partial_regular_log2_upper": partial_log / LOG2,
        "partial_regular_lambda_bits_lower_float": -partial_log / LOG2,
        "dominant_rows": dominant,
        "occupation_rows": rows,
        "scope": (
            "Floating-point density-envelope diagnostic for configurations "
            f"with one through {args.maximum_active_blocks} regular active outer "
            "blocks. Configurations containing the all-one word and larger "
            "occupations are excluded. Arithmetic is not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--epoch-bits", type=int, default=1024)
    parser.add_argument("--epochs-per-region", type=int, default=8)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--maximum-active-blocks", type=int, default=64)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-5.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "partial_regular_lambda_bits,"
        f"{payload['partial_regular_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
