#!/usr/bin/env python3
"""All-profile regular-word bound for reset-free s=256 trajectories.

The checker sums every regular active-block count and every group-rank
histogram.  A trajectory starts in the zero checkpoint state, activates once,
and never returns to zero.  Before activation, epoch factors are 2^-A.  The
activation factor is q^A-2^-A, where q=(1+z)/2.  After activation, every live
epoch is bounded by ((1+z)^256-1)/(2^256-1).

For a fixed active-block count, the outer group sum reduces to coefficients of
three binomial factors.  A positive-radius Cauchy bound encloses each
coefficient.  The radius is optimized in floating point, but validity does not
depend on optimality.  Return-to-zero trajectories and all-one outer words are
outside this receipt.
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
    "receipts/all_profiles_regular_noreset_delta09.json"
)


def log_choose(total: int, count: np.ndarray | float | int) -> np.ndarray:
    values = np.asarray(count, dtype=np.float64)
    return (
        gammaln(total + 1.0)
        - gammaln(values + 1.0)
        - gammaln(total - values + 1.0)
    )


def log_difference(log_larger: np.ndarray, log_smaller: np.ndarray) -> np.ndarray:
    """Return log(exp(log_larger)-exp(log_smaller)) stably."""
    gap = np.minimum(log_smaller - log_larger, -1e-15)
    return log_larger + np.log(-np.expm1(gap))


def activation_coefficient_cauchy_logs(
    *,
    total_blocks: int,
    prefix_blocks: int,
    activation_blocks: int,
    p: float,
    q: float,
    active_counts: np.ndarray,
) -> np.ndarray:
    """Upper-bound all coefficients of the activation polynomial.

    The polynomial is

      (1+p x)^prefix * ((1+q x)^activation-(1+p x)^activation)
      * (1+x)^later.

    Every positive radius gives a valid Cauchy coefficient bound.  Newton
    iteration selects a tight radius separately for each requested degree.
    """
    later_blocks = total_blocks - prefix_blocks - activation_blocks
    degrees = np.asarray(active_counts, dtype=np.float64)
    log_radius = np.log(degrees / (total_blocks - degrees + 1e-12))
    log_radius = np.clip(log_radius, -40.0, 40.0)

    for _ in range(12):
        radius = np.exp(log_radius)
        xp = p * radius / (1.0 + p * radius)
        x1 = radius / (1.0 + radius)
        xq = q * radius / (1.0 + q * radius)
        log_q_term = activation_blocks * np.log1p(q * radius)
        log_p_term = activation_blocks * np.log1p(p * radius)
        ratio = np.exp(np.minimum(log_p_term - log_q_term, -1e-15))
        denominator = 1.0 - ratio

        mean_q = activation_blocks * xq
        mean_p = activation_blocks * xp
        variance_q = activation_blocks * xq * (1.0 - xq)
        variance_p = activation_blocks * xp * (1.0 - xp)
        mean_difference = (mean_q - ratio * mean_p) / denominator
        second_difference = (
            variance_q
            + mean_q * mean_q
            - ratio * (variance_p + mean_p * mean_p)
        ) / denominator
        variance_difference = np.maximum(
            second_difference - mean_difference * mean_difference,
            1e-8,
        )

        mean = (
            prefix_blocks * xp
            + later_blocks * x1
            + mean_difference
        )
        variance = (
            prefix_blocks * xp * (1.0 - xp)
            + later_blocks * x1 * (1.0 - x1)
            + variance_difference
        )
        log_radius = np.clip(
            log_radius + (degrees - mean) / np.maximum(variance, 1e-8),
            -40.0,
            40.0,
        )

    radius = np.exp(log_radius)
    log_q_term = activation_blocks * np.log1p(q * radius)
    log_p_term = activation_blocks * np.log1p(p * radius)
    log_activation_difference = log_difference(log_q_term, log_p_term)
    return (
        prefix_blocks * np.log1p(p * radius)
        + later_blocks * np.log1p(radius)
        + log_activation_difference
        - degrees * log_radius
    )


def exact_activation_log_coefficient(
    *,
    total_blocks: int,
    prefix_blocks: int,
    activation_blocks: int,
    p: float,
    q: float,
    degree: int,
) -> float:
    """Compute a small-degree reference coefficient in the log domain."""
    later_blocks = total_blocks - prefix_blocks - activation_blocks
    terms = []
    for activation_degree in range(1, min(activation_blocks, degree) + 1):
        log_activation = (
            float(log_choose(activation_blocks, activation_degree))
            + float(
                log_difference(
                    np.asarray(activation_degree * math.log(q)),
                    np.asarray(activation_degree * math.log(p)),
                )
            )
        )
        base_degree = degree - activation_degree
        base_terms = []
        lower = max(0, base_degree - later_blocks)
        upper = min(prefix_blocks, base_degree)
        for prefix_degree in range(lower, upper + 1):
            base_terms.append(
                float(log_choose(prefix_blocks, prefix_degree))
                + prefix_degree * math.log(p)
                + float(log_choose(later_blocks, base_degree - prefix_degree))
            )
        terms.append(log_activation + float(logsumexp(np.asarray(base_terms))))
    return float(logsumexp(np.asarray(terms)))


def self_test() -> dict[str, float]:
    maximum_violation = 0.0
    p = 0.5
    q = 0.73
    degrees = np.arange(1, 17, dtype=np.float64)
    for prefix_blocks in (0, 256, 1280, 7936):
        upper = activation_coefficient_cauchy_logs(
            total_blocks=8192,
            prefix_blocks=prefix_blocks,
            activation_blocks=256,
            p=p,
            q=q,
            active_counts=degrees,
        )
        for index, degree in enumerate(range(1, 17)):
            exact = exact_activation_log_coefficient(
                total_blocks=8192,
                prefix_blocks=prefix_blocks,
                activation_blocks=256,
                p=p,
                q=q,
                degree=degree,
            )
            maximum_violation = max(maximum_violation, exact - upper[index])
    if maximum_violation > 2e-10:
        raise AssertionError("Cauchy coefficient bound fell below exact coefficient")
    return {"maximum_log_coefficient_upper_violation": maximum_violation}


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    groups = outer_blocks // args.packet_bits
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    if args.state_bits != args.epoch_bits:
        raise ValueError("the reset-free reduction requires one visit per state lane")
    if groups != args.epoch_packet_slots * args.epochs_per_region:
        raise ValueError("epoch packet geometry does not fill one region")
    total_epochs = args.outer_bits * args.epochs_per_region

    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(dimension) + log_eta

    active_counts = np.arange(1, outer_blocks + 1, dtype=np.float64)
    log_block_choices = log_choose(outer_blocks, active_counts)
    best = np.full(outer_blocks, math.inf)
    best_tilt = np.full(outer_blocks, math.nan)

    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        p = 0.5
        q = (1.0 + z) * 0.5
        live_log_moment = (
            math.log(math.expm1(args.state_bits * math.log1p(z)))
            - math.log(math.ldexp(1.0, args.state_bits) - 1.0)
        )

        # The proxy words never activate during their 256 coordinate bits.
        terms = [
            log_block_choices - args.outer_bits * active_counts * LOG2
        ]

        # If activation occurs in epoch k of region u, the u-dependence is a
        # geometric series.  Only 32 activation-epoch coefficient bounds are
        # required for the complete 256-region sum.
        geometric_ratio_log = (
            -active_counts * LOG2
            - args.epochs_per_region * live_log_moment
        )
        geometric_logs = logsumexp(
            np.arange(args.outer_bits, dtype=np.float64)[:, None]
            * geometric_ratio_log[None, :],
            axis=0,
        )
        blocks_per_epoch = args.epoch_packet_slots * args.packet_bits
        for epoch in range(args.epochs_per_region):
            prefix_blocks = blocks_per_epoch * epoch
            coefficient_logs = activation_coefficient_cauchy_logs(
                total_blocks=outer_blocks,
                prefix_blocks=prefix_blocks,
                activation_blocks=blocks_per_epoch,
                p=p,
                q=q,
                active_counts=active_counts,
            )
            live_tail_epochs = total_epochs - epoch - 1
            terms.append(
                coefficient_logs
                + live_tail_epochs * live_log_moment
                + geometric_logs
            )

        inner_and_placements = np.logaddexp.reduce(np.stack(terms))
        candidates = (
            active_counts * regular_log_mass
            + inner_and_placements
            + distance * surprisal
        )
        improved = candidates < best
        best[improved] = candidates[improved]
        best_tilt[improved] = log_surprisal
        print(
            f"tilt,{tilt_index + 1},{len(args.log_surprisals)},"
            f"log_surprisal,{log_surprisal:.6f}",
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
        "schema": "riffle-packet4-fieldcheckpoint-s256-regular-noreset-v1",
        "candidate": "BCHPerm-TransposePacketShuffle-FieldCheckpoint-s256-g4",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "packet_bits": args.packet_bits,
            "state_bits": args.state_bits,
            "epoch_bits": args.epoch_bits,
            "epoch_packet_slots": args.epoch_packet_slots,
            "epochs_per_region": args.epochs_per_region,
            "total_epochs": total_epochs,
            "log_surprisals": args.log_surprisals,
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
            "special_support_excluded": "the unique all-one outer word",
        },
        "method": {
            "trajectory_class": "one activation and no return to zero",
            "group_histograms": "all histograms summed through block coefficients",
            "coefficient_bound": "positive-radius Cauchy bound",
            "arithmetic": "nearest binary64 log semiring",
            "outward_rounded": False,
        },
        "self_test": self_test(),
        "regular_noreset_log2_upper": total_log / LOG2,
        "regular_noreset_lambda_bits_lower_float": -total_log / LOG2,
        "dominant_rows": sorted(
            rows,
            key=lambda row: float(row["pointwise_log2_upper"]),
            reverse=True,
        )[:30],
        "occupation_rows": rows,
        "scope": (
            "All regular active-block counts and group-rank histograms for "
            "trajectories that never return to zero after activation. "
            "Return-to-zero trajectories, all-one words, and outward rounding "
            "remain open."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--packet-bits", type=int, default=4)
    parser.add_argument("--state-bits", type=int, default=256)
    parser.add_argument("--epoch-bits", type=int, default=256)
    parser.add_argument("--epoch-packet-slots", type=int, default=64)
    parser.add_argument("--epochs-per-region", type=int, default=32)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0,
                 -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0,
                 0.5, 0.75, 1.0),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "regular_noreset_lambda_bits,"
        f"{payload['regular_noreset_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
