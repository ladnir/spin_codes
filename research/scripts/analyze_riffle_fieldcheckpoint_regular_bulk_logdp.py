#!/usr/bin/env python3
"""Exponent-safe full regular-occupation diagnostic in the log semiring."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp

from analyze_riffle_fieldcheckpoint_accumulate_oneblock import (
    occupancy_epoch_matrices,
)
from analyze_riffle_fieldcheckpoint_regular_bulk_tiled import binomial_transform
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
    log_matmul_batch,
    log_matrix_entries,
    log_two_power_minus_one,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/"
    "receipts/goal09_regular_bulk_logdp.json"
)


def log_choose_vector(total: int) -> np.ndarray:
    counts = np.arange(total + 1, dtype=np.float64)
    return (
        gammaln(total + 1.0)
        - gammaln(counts + 1.0)
        - gammaln(total - counts + 1.0)
    )


def candidate_epoch_logs(
    state_bits: int,
    epoch_bits: int,
    z: float,
    transform: np.ndarray,
) -> np.ndarray:
    impulses = np.stack(occupancy_epoch_matrices(state_bits, epoch_bits, z))
    candidates = (transform @ impulses.reshape(epoch_bits + 1, 4)).reshape(
        epoch_bits + 1, 2, 2
    )
    candidates[candidates < 0.0] = 0.0
    with np.errstate(divide="ignore"):
        return np.log(candidates)


def regular_region_log_matrices(
    candidate_logs: np.ndarray,
    epoch_bits: int,
    epochs_per_region: int,
    maximum_occupation: int,
    progress: bool = False,
) -> np.ndarray:
    """Return log averaged region matrices for all candidate occupations."""
    current = candidate_logs.copy()
    current_maximum = epoch_bits
    epoch_choose = log_choose_vector(epoch_bits)
    for completed_epochs in range(1, epochs_per_region):
        next_maximum = min(
            (completed_epochs + 1) * epoch_bits, maximum_occupation
        )
        updated = np.full((next_maximum + 1, 2, 2), -math.inf)
        current_choose = log_choose_vector(completed_epochs * epoch_bits)
        next_choose = log_choose_vector((completed_epochs + 1) * epoch_bits)
        for next_epoch_count in range(epoch_bits + 1):
            source_count = min(
                current_maximum,
                next_maximum - next_epoch_count,
            ) + 1
            if source_count <= 0:
                break
            products = log_matmul_batch(
                current[:source_count], candidate_logs[next_epoch_count]
            )
            destinations = slice(
                next_epoch_count, next_epoch_count + source_count
            )
            weights = (
                epoch_choose[next_epoch_count]
                + current_choose[:source_count]
                - next_choose[next_epoch_count : next_epoch_count + source_count]
            )
            terms = products + weights[:, None, None]
            updated[destinations] = np.logaddexp(updated[destinations], terms)
        current = updated
        current_maximum = next_maximum
        if progress:
            print(
                f"region_epoch,{completed_epochs + 1},{epochs_per_region},"
                f"occupations,{current_maximum + 1}",
                flush=True,
            )
    return current


def log_matrix_power_moments_batch(
    matrices: np.ndarray, power: int
) -> np.ndarray:
    rows = np.full((matrices.shape[0], 2), -math.inf)
    rows[:, 0] = 0.0
    for _ in range(power):
        next_zero = np.logaddexp(
            rows[:, 0] + matrices[:, 0, 0],
            rows[:, 1] + matrices[:, 1, 0],
        )
        next_live = np.logaddexp(
            rows[:, 0] + matrices[:, 0, 1],
            rows[:, 1] + matrices[:, 1, 1],
        )
        rows[:, 0] = next_zero
        rows[:, 1] = next_live
    return np.logaddexp(rows[:, 0], rows[:, 1])


def one_forced_from_regular_logs(regular_logs: np.ndarray) -> np.ndarray:
    """Return F_{a,1}=2 F_{a+1,0}-F_{a,0} entrywise in log form."""
    minuend = LOG2 + regular_logs[1:]
    subtrahend = regular_logs[:-1]
    both_zero = np.isneginf(subtrahend) & np.isneginf(minuend)
    with np.errstate(invalid="ignore"):
        gap = subtrahend - minuend
    gap[both_zero] = -math.inf
    if np.any(gap > 2e-12):
        raise ArithmeticError("one-forced finite difference became negative")
    gap = np.minimum(gap, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        correction = np.log1p(-np.exp(gap))
    return minuend + correction


def self_test(transform: np.ndarray) -> dict[str, float]:
    from analyze_riffle_fieldcheckpoint_mixed_envelope import (
        active_region_matrices,
        mixed_region_matrix,
    )

    z = 0.73
    # Keep the total occupation cap no larger than one toy epoch.  The
    # ordinary-arithmetic reference routine uses that cap both for an epoch
    # and for the whole region, whereas the log recurrence separates them.
    candidate_logs = candidate_epoch_logs(3, 12, z, transform[:13, :13])
    recovered = regular_region_log_matrices(candidate_logs, 12, 2, 12)
    reference = regular_region_matrices(3, 12, 2, 12, z)
    maximum_error = 0.0
    for occupation in range(13):
        recovered_matrix = np.exp(recovered[occupation])
        maximum_error = max(
            maximum_error,
            float(np.max(np.abs(recovered_matrix - reference[occupation]))),
        )
    forced_recovered = one_forced_from_regular_logs(recovered)
    active_reference = active_region_matrices(3, 12, 2, 12, z)
    forced_error = 0.0
    for regular in range(12):
        forced_reference = mixed_region_matrix(active_reference, regular, 1)
        forced_error = max(
            forced_error,
            float(
                np.max(
                    np.abs(np.exp(forced_recovered[regular]) - forced_reference)
                )
            ),
        )
    stochastic_logs = regular_region_log_matrices(
        candidate_epoch_logs(3, 12, 1.0, transform[:13, :13]), 12, 2, 12
    )
    stochastic_error = float(
        np.max(
            np.abs(
                np.exp(logsumexp(stochastic_logs, axis=2)) - 1.0
            )
        )
    )
    if max(maximum_error, forced_error, stochastic_error) > 8e-13:
        raise AssertionError("log-domain region recurrence failed")
    return {
        "maximum_small_reference_error": maximum_error,
        "maximum_small_one_forced_error": forced_error,
        "maximum_small_stochastic_error": stochastic_error,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    region_bits = args.epoch_bits * args.epochs_per_region
    if region_bits != outer_blocks:
        raise ValueError(
            "epoch_bits * epochs_per_region must equal one transposed region"
        )
    if args.epoch_bits % args.state_bits != 0:
        raise ValueError("epoch_bits must be divisible by state_bits")
    checkpoint_steps = args.epoch_bits // args.step_bits
    candidate_name = (
        "Riffle FieldCheckpointAccumulate "
        f"t={args.step_bits} s={args.state_bits} K={checkpoint_steps}"
    )
    maximum_occupation = min(args.maximum_active_blocks, region_bits)
    computed_maximum_occupation = min(
        maximum_occupation + int(args.include_one_all_one), region_bits
    )
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(dimension) + log_eta
    transform = binomial_transform(args.epoch_bits)

    best = np.full(maximum_occupation + 1, math.inf)
    best_tilt = np.full(maximum_occupation + 1, math.nan)
    one_all_one_best = np.full(maximum_occupation + 1, math.inf)
    one_all_one_best_tilt = np.full(maximum_occupation + 1, math.nan)
    moment_grid = []
    surprisals = []
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        candidate_logs = candidate_epoch_logs(
            args.state_bits, args.epoch_bits, z, transform
        )
        regions = regular_region_log_matrices(
            candidate_logs,
            args.epoch_bits,
            args.epochs_per_region,
            computed_maximum_occupation,
            progress=True,
        )
        moments = log_matrix_power_moments_batch(
            regions[: maximum_occupation + 1], args.outer_bits
        )
        moment_grid.append(moments)
        surprisals.append(surprisal)
        candidates = moments + distance * surprisal
        improved = candidates < best
        best[improved] = candidates[improved]
        best_tilt[improved] = log_surprisal
        if args.include_one_all_one:
            one_all_one_regions = one_forced_from_regular_logs(regions)
            one_all_one_moments = log_matrix_power_moments_batch(
                one_all_one_regions[: maximum_occupation + 1], args.outer_bits
            )
            one_all_one_candidates = (
                one_all_one_moments + distance * surprisal
            )
            one_all_one_improved = one_all_one_candidates < one_all_one_best
            one_all_one_best[one_all_one_improved] = one_all_one_candidates[
                one_all_one_improved
            ]
            one_all_one_best_tilt[one_all_one_improved] = log_surprisal
        print(
            f"tilt,{tilt_index + 1},{len(args.log_surprisals)},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    contributions = []
    for occupation in range(args.minimum_active_blocks, maximum_occupation + 1):
        outer_log = (
            log_choose(outer_blocks, occupation)
            + occupation * regular_log_mass
        )
        inner_log = min(0.0, float(best[occupation]))
        contribution = outer_log + inner_log
        contributions.append(contribution)
        rows.append(
            {
                "active_regular_outer_blocks": occupation,
                "best_log_surprisal": float(best_tilt[occupation]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
            }
        )
    partial_log = float(logsumexp(np.asarray(contributions)))
    dominant = sorted(
        rows, key=lambda row: float(row["pointwise_log2_upper"]), reverse=True
    )[:30]
    one_all_one_rows = []
    one_all_one_log = -math.inf
    if args.include_one_all_one:
        one_all_one_contributions = []
        for occupation in range(
            args.minimum_active_blocks, maximum_occupation + 1
        ):
            if occupation >= outer_blocks:
                continue
            outer_log = (
                log_choose(outer_blocks, occupation)
                + math.log(outer_blocks - occupation)
                + occupation * regular_log_mass
            )
            inner_log = min(0.0, float(one_all_one_best[occupation]))
            contribution = outer_log + inner_log
            one_all_one_contributions.append(contribution)
            one_all_one_rows.append(
                {
                    "regular_active_blocks": occupation,
                    "all_one_active_blocks": 1,
                    "best_log_surprisal": float(
                        one_all_one_best_tilt[occupation]
                    ),
                    "outer_log2_envelope": outer_log / LOG2,
                    "inner_log2_upper": inner_log / LOG2,
                    "pointwise_log2_upper": contribution / LOG2,
                }
            )
        one_all_one_log = float(logsumexp(np.asarray(one_all_one_contributions)))
    scan_profiles = []
    if args.scan_relative_distances:
        moment_array = np.stack(moment_grid)
        surprisal_array = np.asarray(surprisals)[:, None]
        for relative_distance in args.scan_relative_distances:
            scan_distance = math.floor(relative_distance * output_bits)
            scan_candidates = moment_array + scan_distance * surprisal_array
            scan_best = np.min(scan_candidates, axis=0)
            scan_best_tilts = np.argmin(scan_candidates, axis=0)
            scan_contributions = []
            scan_rows = []
            for occupation in range(
                args.minimum_active_blocks, maximum_occupation + 1
            ):
                outer_log = (
                    log_choose(outer_blocks, occupation)
                    + occupation * regular_log_mass
                )
                inner_log = min(0.0, float(scan_best[occupation]))
                contribution = outer_log + inner_log
                scan_contributions.append(contribution)
                scan_rows.append(
                    {
                        "active_regular_outer_blocks": occupation,
                        "best_log_surprisal": args.log_surprisals[
                            int(scan_best_tilts[occupation])
                        ],
                        "pointwise_log2_upper": contribution / LOG2,
                    }
                )
            scan_partial_log = float(logsumexp(np.asarray(scan_contributions)))
            scan_dominant = max(
                scan_rows, key=lambda row: float(row["pointwise_log2_upper"])
            )
            scan_profiles.append(
                {
                    "relative_distance": relative_distance,
                    "distance": scan_distance,
                    "partial_log2_upper": scan_partial_log / LOG2,
                    "partial_lambda_bits_lower_float": -scan_partial_log / LOG2,
                    "dominant_row": scan_dominant,
                }
            )

    return {
        "schema": "riffle-fieldcheckpoint-regular-bulk-logdp-v1",
        "candidate": candidate_name,
        "method": {
            "region_recurrence": "exact hypergeometric conditioning after each epoch",
            "arithmetic": "nearest binary64 log semiring",
            "outward_rounded": False,
            "fft_used": False,
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
        },
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "state_bits": args.state_bits,
            "step_bits": args.step_bits,
            "checkpoint_steps": checkpoint_steps,
            "epoch_bits": args.epoch_bits,
            "epochs_per_region": args.epochs_per_region,
            "minimum_active_blocks": args.minimum_active_blocks,
            "maximum_active_blocks": maximum_occupation,
            "log_surprisals": args.log_surprisals,
        },
        "self_test": self_test(transform),
        "partial_log2_upper": partial_log / LOG2,
        "partial_lambda_bits_lower_float": -partial_log / LOG2,
        "dominant_rows": dominant,
        "occupation_rows": rows,
        "one_all_one": {
            "included": args.include_one_all_one,
            "partial_log2_upper": (
                one_all_one_log / LOG2 if args.include_one_all_one else None
            ),
            "partial_lambda_bits_lower_float": (
                -one_all_one_log / LOG2 if args.include_one_all_one else None
            ),
            "dominant_rows": sorted(
                one_all_one_rows,
                key=lambda row: float(row["pointwise_log2_upper"]),
                reverse=True,
            )[:30],
            "occupation_rows": one_all_one_rows,
        },
        "distance_scan_profiles": scan_profiles,
        "scope": (
            "Complete exponent-safe floating-point diagnostic over the reported "
            "regular occupation range. The log semiring is not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--step-bits", type=int, default=32)
    parser.add_argument("--epoch-bits", type=int, default=1024)
    parser.add_argument("--epochs-per-region", type=int, default=8)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--minimum-active-blocks", type=int, default=129)
    parser.add_argument("--maximum-active-blocks", type=int, default=1024)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-5.0, -4.5, -4.0, -3.5, -3.0),
    )
    parser.add_argument(
        "--scan-relative-distances",
        type=float,
        nargs="+",
        default=(),
        help="reuse every computed moment to report additional distance profiles",
    )
    parser.add_argument(
        "--include-one-all-one",
        action="store_true",
        help="also evaluate configurations with one exact all-one outer word",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "regular_bulk_lambda_bits,"
        f"{payload['partial_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
