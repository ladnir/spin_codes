#!/usr/bin/env python3
"""All-trajectory regular-word bound using a worst one-region allocation.

For every active-block count a, this checker upper-bounds a complete region by
maximizing over all 32-epoch allocations A_0+...+A_31=a with 0<=A_i<=256.
It deliberately forgets the favorable packet-shuffle probability.  A
max-convolution dynamic program retains the exact two-state epoch transfer and
therefore charges both output weight before a reset and every later
reactivation.  Entrywise maxima are propagated independently, which is an
additional valid relaxation.  The resulting 2x2 region matrix is raised to
the 256th power and combined with the complete regular outer multiplicity.

Arithmetic is nearest binary64.  The unique all-one outer word and outward
rounding are outside this receipt.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp

from analyze_riffle_fieldcheckpoint_regular_envelope import (
    spectrum_density_envelope_log,
)
from analyze_riffle_packet4_fieldcheckpoint_s256_holder import log_matmul, log_matpow
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import LOG2, log_two_power_minus_one


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_packetshuffle_fieldcheckpoint_g4_s256/"
    "receipts/all_profiles_regular_all_trajectories_regionworst_delta09.json"
)


def log_difference(log_larger: float, log_smaller: float) -> float:
    if log_smaller == -math.inf:
        return log_larger
    if not log_smaller < log_larger:
        return -math.inf
    return log_larger + math.log1p(-math.exp(log_smaller - log_larger))


def epoch_transfer_logs(z: float, state_bits: int, maximum_a: int) -> np.ndarray:
    transfers = np.full((maximum_a + 1, 2, 2), -math.inf, dtype=np.float64)
    log_d = math.log(math.ldexp(1.0, state_bits) - 1.0)
    log_s_minus_one = math.log(math.expm1(state_bits * math.log1p(z)))
    log_q_base = math.log1p(z) - LOG2
    for active in range(maximum_a + 1):
        log_p = -active * LOG2
        log_q = active * log_q_base
        transfers[active, 0, 0] = log_p
        if active:
            transfers[active, 0, 1] = log_difference(log_q, log_p)
            transfers[active, 1, 0] = log_difference(0.0, log_p) - log_d
        log_q_minus_p = (
            log_difference(log_q, log_p) if active else -math.inf
        )
        transfers[active, 1, 1] = (
            log_difference(log_s_minus_one, log_q_minus_p) - log_d
        )
    return transfers


def worst_region_logs(
    z: float, *, state_bits: int, epoch_bits: int, epochs_per_region: int
) -> np.ndarray:
    """Entrywise upper bound over every epoch allocation at every total a."""
    maximum_total = epoch_bits * epochs_per_region
    transfers = epoch_transfer_logs(z, state_bits, epoch_bits)
    previous = np.full((maximum_total + 1, 2, 2), -math.inf, dtype=np.float64)
    previous[0, 0, 0] = 0.0
    previous[0, 1, 1] = 0.0
    previous_limit = 0
    for _ in range(epochs_per_region):
        current_limit = previous_limit + epoch_bits
        current = np.full_like(previous, -math.inf)
        for active in range(epoch_bits + 1):
            source = previous[: previous_limit + 1]
            factor = np.broadcast_to(transfers[active], source.shape)
            candidate = log_matmul(source, factor)
            target = current[active : active + previous_limit + 1]
            np.maximum(target, candidate, out=target)
        previous = current
        previous_limit = current_limit
    return previous


def self_test() -> dict[str, float]:
    """The relaxed DP must dominate exact enumeration in a toy region."""
    z = 0.57
    epoch_bits = 2
    epochs = 3
    upper = worst_region_logs(
        z, state_bits=5, epoch_bits=epoch_bits, epochs_per_region=epochs
    )
    transfers = np.exp(epoch_transfer_logs(z, 5, epoch_bits))
    maximum_violation = 0.0
    for a0 in range(epoch_bits + 1):
        for a1 in range(epoch_bits + 1):
            for a2 in range(epoch_bits + 1):
                total = a0 + a1 + a2
                exact = transfers[a0] @ transfers[a1] @ transfers[a2]
                finite = exact > 0.0
                maximum_violation = max(
                    maximum_violation,
                    float(np.max(np.log(exact[finite]) - upper[total][finite])),
                )
    if maximum_violation > 2e-12:
        raise AssertionError("region maximum DP fell below an exact allocation")
    return {"maximum_log_upper_violation": maximum_violation}


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    if args.epoch_bits * args.epochs_per_region != outer_blocks:
        raise ValueError("one region must contain exactly one slot per outer block")

    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(dimension) + log_eta
    active_counts = np.arange(1, outer_blocks + 1, dtype=np.float64)
    log_block_choices = (
        gammaln(outer_blocks + 1.0)
        - gammaln(active_counts + 1.0)
        - gammaln(outer_blocks - active_counts + 1.0)
    )

    best = np.full(outer_blocks, math.inf)
    best_tilt = np.full(outer_blocks, math.nan)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        region = worst_region_logs(
            z,
            state_bits=args.state_bits,
            epoch_bits=args.epoch_bits,
            epochs_per_region=args.epochs_per_region,
        )[1:]
        all_regions = log_matpow(region, args.outer_bits)
        inner = np.logaddexp(all_regions[:, 0, 0], all_regions[:, 0, 1])
        candidates = (
            log_block_choices
            + active_counts * regular_log_mass
            + inner
            + distance * surprisal
        )
        improved = candidates < best
        best[improved] = candidates[improved]
        best_tilt[improved] = log_surprisal
        print(
            f"tilt,{tilt_index + 1},{len(args.log_surprisals)},"
            f"log_surprisal,{log_surprisal:.6f},"
            f"current_lambda,{-float(logsumexp(best)) / LOG2:.6f}",
            flush=True,
        )

    total_log = float(logsumexp(best))
    rows = [
        {
            "active_regular_outer_blocks": active_blocks,
            "best_log_surprisal": float(best_tilt[active_blocks - 1]),
            "pointwise_log2_upper": float(best[active_blocks - 1]) / LOG2,
        }
        for active_blocks in range(1, outer_blocks + 1)
    ]
    return {
        "schema": "riffle-packet4-fieldcheckpoint-s256-regular-regionworst-v1",
        "candidate": "BCHPerm-TransposePacketShuffle-FieldCheckpoint-s256-g4",
        "parameters": vars(args) | {
            "output_bits": output_bits,
            "distance": distance,
            "outer_blocks": outer_blocks,
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
            "special_support_excluded": "the unique all-one outer word",
        },
        "method": {
            "trajectory_class": "all zero/live trajectories",
            "packet_allocation": "entrywise worst case over every 32-epoch allocation",
            "arithmetic": "nearest binary64 log semiring",
            "outward_rounded": False,
        },
        "self_test": self_test(),
        "regular_all_trajectories_log2_upper": total_log / LOG2,
        "regular_all_trajectories_lambda_bits_lower_float": -total_log / LOG2,
        "dominant_rows": sorted(
            rows,
            key=lambda row: float(row["pointwise_log2_upper"]),
            reverse=True,
        )[:30],
        "occupation_rows": rows,
        "scope": (
            "Every regular active-block count, packet-group histogram, and state "
            "trajectory.  The unique all-one outer word and outward rounding "
            "remain open."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--state-bits", type=int, default=256)
    parser.add_argument("--epoch-bits", type=int, default=256)
    parser.add_argument("--epochs-per-region", type=int, default=32)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
                 -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0,
                 -1.5, -1.0, -0.5, 0.0, 0.5, 0.75, 1.0),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    output = args.output
    payload["parameters"]["output"] = str(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    print(
        "regular_all_trajectories_lambda_bits,"
        f"{payload['regular_all_trajectories_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{output}")


if __name__ == "__main__":
    main()
