#!/usr/bin/env python3
"""Floating-point all-trajectory regular-word bound for the s=256 variant.

The 256 outer coordinates define 256 statistically identical regions.  For a
fixed active-block subset, condition on the zero/live state at every region
boundary and apply Holder's inequality across the regions.  The required
256th one-region moments are bounded by their first moments because every
entry is in [0,1].  Summing the resulting two-state boundary paths is a 2x2
matrix power.

For one region, the grand-canonical active-subset sum factors over its 32
epochs.  Each epoch contains 256 block slots, so its matrix is

  M(x,z) = sum_A binom(256,A) x^A F_A(z).

Thus M(x,z)^32 generates the one-region first moments.  A positive-radius
Cauchy bound encloses each active-count coefficient.  The radius grid is only
used to tighten the bound; every grid point is independently valid.

This checker covers every regular active-block count, every packet-group rank
histogram, and every zero/live state trajectory.  The unique all-one outer
word is deliberately excluded.
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
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import LOG2, log_two_power_minus_one


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_packetshuffle_fieldcheckpoint_g4_s256/"
    "receipts/all_profiles_regular_all_trajectories_holder_delta09.json"
)


def log_difference(log_larger: float, log_smaller: float) -> float:
    if log_smaller == -math.inf:
        return log_larger
    if not log_smaller < log_larger:
        return -math.inf
    return log_larger + math.log1p(-math.exp(log_smaller - log_larger))


def log_matmul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Log-semiring product, with arbitrary leading batch dimensions."""
    result = np.empty_like(left)
    result[..., 0, 0] = np.logaddexp(
        left[..., 0, 0] + right[..., 0, 0],
        left[..., 0, 1] + right[..., 1, 0],
    )
    result[..., 0, 1] = np.logaddexp(
        left[..., 0, 0] + right[..., 0, 1],
        left[..., 0, 1] + right[..., 1, 1],
    )
    result[..., 1, 0] = np.logaddexp(
        left[..., 1, 0] + right[..., 0, 0],
        left[..., 1, 1] + right[..., 1, 0],
    )
    result[..., 1, 1] = np.logaddexp(
        left[..., 1, 0] + right[..., 0, 1],
        left[..., 1, 1] + right[..., 1, 1],
    )
    return result


def log_matpow(matrix: np.ndarray, exponent: int) -> np.ndarray:
    shape = matrix.shape
    result = np.full(shape, -math.inf, dtype=np.float64)
    result[..., 0, 0] = 0.0
    result[..., 1, 1] = 0.0
    base = matrix.copy()
    power = exponent
    while power:
        if power & 1:
            result = log_matmul(result, base)
        power >>= 1
        if power:
            base = log_matmul(base, base)
    return result


def epoch_log_matrix(
    *, log_radius: float, z: float, epoch_bits: int, state_bits: int
) -> np.ndarray:
    """Return log M(x,z) by direct positive summation over A."""
    counts = np.arange(epoch_bits + 1, dtype=np.float64)
    log_choose = (
        gammaln(epoch_bits + 1.0)
        - gammaln(counts + 1.0)
        - gammaln(epoch_bits - counts + 1.0)
    )
    log_weights = log_choose + counts * log_radius
    log_d = math.log(math.ldexp(1.0, state_bits) - 1.0)
    log_s = state_bits * math.log1p(z)
    log_s_minus_one = math.log(math.expm1(log_s))
    log_q_base = math.log1p(z) - LOG2

    entries = np.full((4, epoch_bits + 1), -math.inf, dtype=np.float64)
    for index, count_float in enumerate(counts):
        count = int(count_float)
        log_p = -count * LOG2
        log_q = count * log_q_base
        entries[0, index] = log_p
        if count:
            entries[1, index] = log_difference(log_q, log_p)
            entries[2, index] = log_difference(0.0, log_p) - log_d

        # S - q - 1 + p = (S-1) - (q-p).
        log_q_minus_p = (
            log_difference(log_q, log_p) if count else -math.inf
        )
        entries[3, index] = (
            log_difference(log_s_minus_one, log_q_minus_p) - log_d
        )

    sums = logsumexp(log_weights[None, :] + entries, axis=1)
    return np.asarray([[sums[0], sums[1]], [sums[2], sums[3]]])


def one_region_log_generator(
    *, log_radius: float, z: float, epoch_bits: int, state_bits: int,
    epochs_per_region: int,
) -> np.ndarray:
    epoch = epoch_log_matrix(
        log_radius=log_radius,
        z=z,
        epoch_bits=epoch_bits,
        state_bits=state_bits,
    )
    return log_matpow(epoch, epochs_per_region)


def self_test() -> dict[str, float]:
    """Compare the generator with exhaustive active subsets in a toy region."""
    epoch_bits = 2
    epochs = 2
    z = 0.61
    x = 0.37
    generated = np.exp(
        one_region_log_generator(
            log_radius=math.log(x),
            z=z,
            epoch_bits=epoch_bits,
            state_bits=5,
            epochs_per_region=epochs,
        )
    )

    # Directly enumerate the active count in each of the two disjoint epochs.
    direct = np.zeros((2, 2), dtype=np.float64)
    for a0 in range(epoch_bits + 1):
        for a1 in range(epoch_bits + 1):
            matrices = []
            for count in (a0, a1):
                p = 2.0 ** (-count)
                q = ((1.0 + z) * 0.5) ** count
                d = 2.0**5 - 1.0
                s = (1.0 + z) ** 5
                matrices.append(
                    np.asarray(
                        [[p, q - p], [(1.0 - p) / d, (s - q - 1.0 + p) / d]]
                    )
                )
            multiplicity = (
                math.comb(epoch_bits, a0)
                * math.comb(epoch_bits, a1)
                * x ** (a0 + a1)
            )
            direct += multiplicity * (matrices[0] @ matrices[1])
    relative_error = float(
        np.max(np.abs(generated - direct) / np.maximum(direct, 1e-300))
    )
    if relative_error > 2e-12:
        raise AssertionError("one-region generator failed exhaustive toy check")
    return {"maximum_relative_generator_error": relative_error}


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    regions = args.outer_bits
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
    best = np.full(outer_blocks, math.inf)
    best_tilt = np.full(outer_blocks, math.nan)
    best_radius = np.full((4, outer_blocks), math.nan)

    radius_grid = np.linspace(
        args.log_radius_min, args.log_radius_max, args.radius_grid_points
    )
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        coefficient_upper = np.full((4, outer_blocks), math.inf)
        coefficient_radius = np.full((4, outer_blocks), math.nan)

        for log_radius in radius_grid:
            region = one_region_log_generator(
                log_radius=float(log_radius),
                z=z,
                epoch_bits=args.epoch_bits,
                state_bits=args.state_bits,
                epochs_per_region=args.epochs_per_region,
            ).reshape(4)
            candidates = region[:, None] - active_counts[None, :] * log_radius
            improved = candidates < coefficient_upper
            coefficient_upper[improved] = candidates[improved]
            coefficient_radius[improved] = log_radius

        # Holder across regions.  Each entry is the Rth root of its regular
        # outer mass times the one-region coefficient upper bound.
        boundary = (
            coefficient_upper.T.reshape(outer_blocks, 2, 2)
            + active_counts[:, None, None] * regular_log_mass
        ) / regions
        all_paths = log_matpow(boundary, regions)
        inner_outer = np.logaddexp(all_paths[:, 0, 0], all_paths[:, 0, 1])
        candidates = inner_outer + distance * surprisal
        improved = candidates < best
        best[improved] = candidates[improved]
        best_tilt[improved] = log_surprisal
        best_radius[:, improved] = coefficient_radius[:, improved]
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
            "best_log_radii_by_boundary_entry": [
                float(value) for value in best_radius[:, active_blocks - 1]
            ],
            "pointwise_log2_upper": float(best[active_blocks - 1]) / LOG2,
        }
        for active_blocks in range(1, outer_blocks + 1)
    ]
    return {
        "schema": "riffle-packet4-fieldcheckpoint-s256-regular-holder-v1",
        "candidate": "BCHPerm-TransposePacketShuffle-FieldCheckpoint-s256-g4",
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
            "regions": regions,
            "log_surprisals": args.log_surprisals,
            "log_radius_min": args.log_radius_min,
            "log_radius_max": args.log_radius_max,
            "radius_grid_points": args.radius_grid_points,
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
            "special_support_excluded": "the unique all-one outer word",
        },
        "method": {
            "trajectory_class": "all zero/live trajectories",
            "group_histograms": "all histograms summed by a Cauchy coefficient bound",
            "region_inequality": "Holder across 256 regions, then E[X]^256 <= E[X]",
            "coefficient_bound": "positive-radius Cauchy grid",
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
            "Every regular active-block count, group-rank histogram, and state "
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
    parser.add_argument("--log-radius-min", type=float, default=-14.0)
    parser.add_argument("--log-radius-max", type=float, default=14.0)
    parser.add_argument("--radius-grid-points", type=int, default=1401)
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
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "regular_all_trajectories_lambda_bits,"
        f"{payload['regular_all_trajectories_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
