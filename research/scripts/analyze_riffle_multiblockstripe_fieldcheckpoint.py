#!/usr/bin/env python3
"""Focused 9% diagnostic for MultiBlockStripe-q FieldCheckpoint.

The regular-spectrum envelope replaces every active outer word by 256 fair
candidate bits.  MultiBlockStripe-q divides each 8192-position transposed
region into q consecutive lane slices.  Inside a slice, the candidates are a
uniform subset; their distribution among slices may retain q-block structure.

For one active outer block we evaluate the proposed cyclic lane schedule
exactly.  For general occupation a we use a proof-safe relaxation: in every
region an adversary may choose any composition of a candidates among the q
slices.  Entrywise maxima during the slice product give a nonnegative matrix
that dominates every such region transfer, hence its 256th power dominates
arbitrary correlations between regions.

This is a floating-point proof diagnostic, not an outward-rounded certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from analyze_riffle_fieldcheckpoint_accumulate_oneblock import (
    epoch_matrices,
)
from analyze_riffle_fieldcheckpoint_regular_envelope import (
    regular_region_matrices,
    spectrum_density_envelope_log,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import (
    LOG2,
    log_choose,
    log_two_power_minus_one,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_multiblockstripe_fieldcheckpoint/"
    "receipts/focused_delta09.json"
)


def scaled_row_product(matrices: list[np.ndarray]) -> float:
    """Return log(e_0^T product(matrices) 1) with rescaling."""
    row = np.asarray((1.0, 0.0), dtype=np.float64)
    accumulated = 0.0
    for matrix in matrices:
        row = row @ matrix
        scale = float(np.max(row))
        if scale <= 0.0:
            return -math.inf
        row /= scale
        accumulated += math.log(scale)
    return accumulated + math.log(float(row.sum()))


def matrix_power_log_moment(matrix: np.ndarray, power: int) -> float:
    return scaled_row_product([matrix] * power)


def exact_one_active_log_moment(
    *,
    q: int,
    outer_bits: int,
    epoch_bits: int,
    epochs_per_region: int,
    state_bits: int,
    z: float,
) -> float:
    """Exact moment for one fair candidate following the cyclic lane rule."""
    zero_epoch, one_epoch = epoch_matrices(state_bits, epoch_bits, z)
    epochs_per_lane = epochs_per_region // q
    zero_powers = [np.eye(2)]
    for _ in range(epochs_per_region):
        zero_powers.append(zero_powers[-1] @ zero_epoch)

    inactive = zero_powers[epochs_per_region]
    lane_matrices: list[np.ndarray] = []
    for lane in range(q):
        first = lane * epochs_per_lane
        active = np.zeros((2, 2), dtype=np.float64)
        for epoch in range(first, first + epochs_per_lane):
            active += (
                zero_powers[epoch]
                @ one_epoch
                @ zero_powers[epochs_per_region - epoch - 1]
            )
        active /= epochs_per_lane
        # The candidate outer bit is fair under the spectrum envelope.
        lane_matrices.append(0.5 * (inactive + active))

    # The proposed route advances the lane by one for every outer coordinate.
    # Average the setup shift outside the complete-word moment.
    shifted_logs = []
    for shift in range(q):
        matrices = [lane_matrices[(coordinate + shift) % q] for coordinate in range(outer_bits)]
        shifted_logs.append(scaled_row_product(matrices))
    return float(logsumexp(np.asarray(shifted_logs))) - math.log(q)


def dominating_region_matrices(
    *, q: int, lane_matrices: list[np.ndarray], maximum_occupation: int
) -> list[np.ndarray]:
    """Entrywise upper bounds for all q-lane compositions of each occupation."""
    negative = np.full((2, 2), -math.inf, dtype=np.float64)
    current = [negative.copy() for _ in range(maximum_occupation + 1)]
    current[0] = np.eye(2)
    for _ in range(q):
        updated = [negative.copy() for _ in range(maximum_occupation + 1)]
        for total in range(maximum_occupation + 1):
            upper = np.zeros((2, 2), dtype=np.float64)
            have_value = False
            for here in range(total + 1):
                before = current[total - here]
                if not np.isfinite(before).any():
                    continue
                candidate = before @ lane_matrices[here]
                if not have_value:
                    upper = candidate
                    have_value = True
                else:
                    upper = np.maximum(upper, candidate)
            if have_value:
                updated[total] = upper
        current = updated
    return current


def independent_lane_region_matrices(
    *, q: int, lane_matrices: list[np.ndarray], maximum_occupation: int
) -> list[np.ndarray]:
    """Average iid lane choices, retaining without-replacement positions per lane.

    This is exact when the active blocks lie in distinct q-block groups and
    every group samples an independent lane permutation (or shift) in every
    outer-coordinate region.  Coefficients use the exponential generating
    function: a! q^-a [u^a] (sum_k L_k u^k/k!)^q.
    """
    coefficients = [np.zeros((2, 2), dtype=np.float64) for _ in range(maximum_occupation + 1)]
    coefficients[0] = np.eye(2)
    inverse_factorials = [1.0 / math.factorial(k) for k in range(maximum_occupation + 1)]
    for _ in range(q):
        updated = [np.zeros((2, 2), dtype=np.float64) for _ in range(maximum_occupation + 1)]
        for total in range(maximum_occupation + 1):
            for here in range(total + 1):
                updated[total] += (
                    coefficients[total - here]
                    @ lane_matrices[here]
                    * inverse_factorials[here]
                )
        coefficients = updated
    return [
        coefficients[a] * math.factorial(a) / (q**a)
        for a in range(maximum_occupation + 1)
    ]


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    if args.outer_bits != 256:
        raise ValueError("the focused schedule currently assumes 256 outer coordinates")
    if args.region_bits % args.epoch_bits:
        raise ValueError("epoch size must divide a transposed region")
    epochs_per_region = args.region_bits // args.epoch_bits
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

    q_values = [int(value) for value in args.q_values.split(",")]
    results = []
    for q in q_values:
        if epochs_per_region % q:
            raise ValueError(f"q={q} must divide {epochs_per_region} epochs per region")
        lane_epochs = epochs_per_region // q
        best_exact_one = math.inf
        best_exact_one_tilt = math.nan
        best_relaxed = np.full(args.maximum_active_blocks + 1, math.inf)
        best_relaxed_tilt = np.full(args.maximum_active_blocks + 1, math.nan)
        best_independent = np.full(args.maximum_active_blocks + 1, math.inf)
        best_independent_tilt = np.full(args.maximum_active_blocks + 1, math.nan)

        grid_count = int(
            math.floor((args.grid_max - args.grid_min) / args.grid_step + 0.5)
        ) + 1
        for grid_index in range(grid_count):
            log_surprisal = args.grid_min + grid_index * args.grid_step
            surprisal = math.exp(log_surprisal)
            z = math.exp(-surprisal)

            exact_one = exact_one_active_log_moment(
                q=q,
                outer_bits=args.outer_bits,
                epoch_bits=args.epoch_bits,
                epochs_per_region=epochs_per_region,
                state_bits=args.state_bits,
                z=z,
            ) + distance * surprisal
            if exact_one < best_exact_one:
                best_exact_one = exact_one
                best_exact_one_tilt = log_surprisal

            lane = regular_region_matrices(
                args.state_bits,
                args.epoch_bits,
                lane_epochs,
                args.maximum_active_blocks,
                z,
            )
            dominating = dominating_region_matrices(
                q=q,
                lane_matrices=lane,
                maximum_occupation=args.maximum_active_blocks,
            )
            independent = independent_lane_region_matrices(
                q=q,
                lane_matrices=lane,
                maximum_occupation=args.maximum_active_blocks,
            )
            for occupation in range(1, args.maximum_active_blocks + 1):
                candidate = (
                    matrix_power_log_moment(dominating[occupation], args.outer_bits)
                    + distance * surprisal
                )
                if candidate < best_relaxed[occupation]:
                    best_relaxed[occupation] = candidate
                    best_relaxed_tilt[occupation] = log_surprisal
                independent_candidate = (
                    matrix_power_log_moment(independent[occupation], args.outer_bits)
                    + distance * surprisal
                )
                if independent_candidate < best_independent[occupation]:
                    best_independent[occupation] = independent_candidate
                    best_independent_tilt[occupation] = log_surprisal

        rows = []
        contributions = []
        for occupation in range(1, args.maximum_active_blocks + 1):
            outer_log = (
                log_choose(outer_blocks, occupation)
                + occupation * per_block_log_mass
            )
            inner_log = min(0.0, float(best_relaxed[occupation]))
            contribution = outer_log + inner_log
            contributions.append(contribution)
            rows.append(
                {
                    "active_regular_outer_blocks": occupation,
                    "best_log_surprisal": float(best_relaxed_tilt[occupation]),
                    "outer_log2_envelope": outer_log / LOG2,
                    "inner_log2_upper": inner_log / LOG2,
                    "pointwise_log2_upper": contribution / LOG2,
                }
            )
        exact_one_contribution = (
            math.log(outer_blocks) + per_block_log_mass + min(0.0, best_exact_one)
        )
        partial_log = float(logsumexp(np.asarray(contributions)))
        independent_rows = []
        independent_contributions = []
        for occupation in range(1, args.maximum_active_blocks + 1):
            outer_log = (
                log_choose(outer_blocks, occupation)
                + occupation * per_block_log_mass
            )
            inner_log = min(0.0, float(best_independent[occupation]))
            contribution = outer_log + inner_log
            independent_contributions.append(contribution)
            independent_rows.append(
                {
                    "active_regular_outer_blocks": occupation,
                    "best_log_surprisal": float(best_independent_tilt[occupation]),
                    "outer_log2_envelope": outer_log / LOG2,
                    "inner_log2_upper": inner_log / LOG2,
                    "pointwise_log2_upper": contribution / LOG2,
                }
            )
        independent_partial_log = float(logsumexp(np.asarray(independent_contributions)))
        results.append(
            {
                "q": q,
                "lane_bits": args.region_bits // q,
                "epochs_per_lane": lane_epochs,
                "exact_one_active": {
                    "best_log_surprisal": best_exact_one_tilt,
                    "pointwise_log2_upper": exact_one_contribution / LOG2,
                    "lambda_bits_lower_float": -exact_one_contribution / LOG2,
                },
                "adversarial_per_region_relaxation": {
                    "maximum_active_blocks": args.maximum_active_blocks,
                    "partial_log2_upper": partial_log / LOG2,
                    "partial_lambda_bits_lower_float": -partial_log / LOG2,
                    "dominant_rows": sorted(
                        rows,
                        key=lambda row: float(row["pointwise_log2_upper"]),
                        reverse=True,
                    )[:10],
                    "occupation_rows": rows,
                },
                "independent_per_coordinate_lane_model": {
                    "status": (
                        "exact for active blocks in distinct q-block groups; "
                        "shared-group domination remains a proof obligation"
                    ),
                    "maximum_active_blocks": args.maximum_active_blocks,
                    "partial_log2_upper": independent_partial_log / LOG2,
                    "partial_lambda_bits_lower_float": -independent_partial_log / LOG2,
                    "dominant_rows": sorted(
                        independent_rows,
                        key=lambda row: float(row["pointwise_log2_upper"]),
                        reverse=True,
                    )[:10],
                    "occupation_rows": independent_rows,
                },
            }
        )
        print(
            f"q,{q},exact_one_lambda,{-exact_one_contribution / LOG2:.9f},"
            f"relaxed_partial_lambda,{-partial_log / LOG2:.9f},"
            f"independent_lane_partial_lambda,{-independent_partial_log / LOG2:.9f}",
            flush=True,
        )

    return {
        "schema": "riffle-multiblockstripe-fieldcheckpoint-focused-v1",
        "candidate": "Riffle MultiBlockStripe-q FieldCheckpoint",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "region_bits": args.region_bits,
            "state_bits": args.state_bits,
            "epoch_bits": args.epoch_bits,
            "epochs_per_region": epochs_per_region,
            "maximum_active_blocks": args.maximum_active_blocks,
            "q_values": q_values,
        },
        "spectrum_envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_weight": envelope_weight,
            "special_support_excluded": "the unique all-one outer word",
        },
        "results": results,
        "scope": (
            "Exact cyclic-lane calculation for one regular active block and a "
            "proof-safe adversarial per-region lane-composition relaxation for "
            "the requested low occupations. Floating-point arithmetic is not "
            "outward rounded; all-one words and higher occupations are excluded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--region-bits", type=int, default=8192)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--epoch-bits", type=int, default=256)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--q-values", default="1,2,4,8,16,32")
    parser.add_argument("--maximum-active-blocks", type=int, default=16)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-2.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
