#!/usr/bin/env python3
"""Exponent-safe tiled bulk diagnostic for regular outer occupations.

Each tile exponentially tilts the candidate-count polynomial so coefficients
near one occupation remain macroscopic.  Three normalized matrix-polynomial
squares combine the eight epochs in a region.  Overlapping tiles provide a
numerical consistency check.  This is an exploratory checker, not yet an
outward-rounded certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve
from scipy.special import gammaln

from analyze_riffle_fieldcheckpoint_accumulate_oneblock import (
    occupancy_epoch_matrices,
)
from analyze_riffle_fieldcheckpoint_regular_envelope import (
    log_row_power_moment,
    spectrum_density_envelope_log,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import LOG2, log_choose, log_two_power_minus_one


DEFAULT_OUTPUT = Path(
    "constructions/riffle_fieldcheckpoint_accumulate_t32_s64_k32/"
    "receipts/goal08_regular_bulk_tiled.json"
)


def log_binomial_pmf(trials: int, probability: float) -> np.ndarray:
    counts = np.arange(trials + 1, dtype=np.float64)
    return (
        gammaln(trials + 1.0)
        - gammaln(counts + 1.0)
        - gammaln(trials - counts + 1.0)
        + counts * math.log(probability)
        + (trials - counts) * math.log1p(-probability)
    )


def binomial_transform(maximum: int) -> np.ndarray:
    """Return T[c,r] = C(c,r)/2^c."""
    transform = np.zeros((maximum + 1, maximum + 1), dtype=np.float64)
    for candidates in range(maximum + 1):
        active = np.arange(candidates + 1, dtype=np.float64)
        logs = (
            gammaln(candidates + 1.0)
            - gammaln(active + 1.0)
            - gammaln(candidates - active + 1.0)
            - candidates * LOG2
        )
        transform[candidates, : candidates + 1] = np.exp(logs)
    return transform


def candidate_epoch_matrices_all(
    state_bits: int,
    epoch_bits: int,
    z: float,
    transform: np.ndarray,
) -> np.ndarray:
    impulses = np.stack(occupancy_epoch_matrices(state_bits, epoch_bits, z))
    flattened = impulses.reshape(epoch_bits + 1, 4)
    return (transform @ flattened).reshape(epoch_bits + 1, 2, 2)


def matrix_polynomial_square(
    polynomial: np.ndarray, log_scale: float
) -> tuple[np.ndarray, float, float]:
    degree = 2 * polynomial.shape[2] - 1
    result = np.zeros((2, 2, degree), dtype=np.float64)
    for row in range(2):
        for column in range(2):
            for inner in range(2):
                result[row, column] += fftconvolve(
                    polynomial[row, inner], polynomial[inner, column]
                )
    negative_floor = float(np.min(result))
    result[result < 0.0] = 0.0
    maximum = float(np.max(result))
    if maximum <= 0.0:
        raise ArithmeticError("matrix polynomial vanished")
    result /= maximum
    return result, 2.0 * log_scale + math.log(maximum), negative_floor


def region_polynomial_tile(
    candidate_epoch: np.ndarray,
    epoch_bits: int,
    center: int,
    region_bits: int,
) -> tuple[np.ndarray, float, np.ndarray, float]:
    probability = min(1.0 - 1.0 / region_bits, max(1.0 / region_bits, center / region_bits))
    epoch_log_pmf = log_binomial_pmf(epoch_bits, probability)
    epoch_pmf = np.exp(epoch_log_pmf)
    polynomial = np.moveaxis(candidate_epoch * epoch_pmf[:, None, None], 0, 2)
    maximum = float(np.max(polynomial))
    polynomial /= maximum
    log_scale = math.log(maximum)
    negative_floor = 0.0
    for _ in range(3):
        polynomial, log_scale, stage_floor = matrix_polynomial_square(
            polynomial, log_scale
        )
        negative_floor = min(negative_floor, stage_floor)
    return polynomial, log_scale, log_binomial_pmf(region_bits, probability), negative_floor


def adaptive_centers(region_bits: int, minimum: int) -> list[int]:
    centers = []
    center = minimum
    while center < region_bits - 1:
        centers.append(center)
        probability = center / region_bits
        deviation = math.sqrt(region_bits * probability * (1.0 - probability))
        center += max(16, int(round(4.0 * deviation)))
    if centers[-1] != region_bits - 1:
        centers.append(region_bits - 1)
    return centers


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    region_bits = args.epoch_bits * args.epochs_per_region
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(dimension) + log_eta
    transform = binomial_transform(args.epoch_bits)
    centers = adaptive_centers(region_bits, args.minimum_active_blocks)

    best = np.full(region_bits + 1, math.inf)
    best_tilt = np.full(region_bits + 1, math.nan)
    best_center = np.full(region_bits + 1, -1, dtype=np.int64)
    maximum_overlap_gap_bits = 0.0
    minimum_fft_coefficient = 0.0
    evaluation_counts = np.zeros(region_bits + 1, dtype=np.int64)

    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        candidate_epoch = candidate_epoch_matrices_all(
            args.state_bits, args.epoch_bits, z, transform
        )
        estimates_by_occupation: dict[int, list[float]] = {}
        denominators_by_occupation: dict[int, list[float]] = {}
        for center_index, center in enumerate(centers):
            polynomial, tile_log_scale, denominator_logs, negative_floor = (
                region_polynomial_tile(
                    candidate_epoch, args.epoch_bits, center, region_bits
                )
            )
            minimum_fft_coefficient = min(minimum_fft_coefficient, negative_floor)
            covered = np.flatnonzero(denominator_logs >= args.minimum_log_denominator)
            for occupation in covered:
                if occupation < args.minimum_active_blocks:
                    continue
                matrix = polynomial[:, :, occupation]
                if float(np.max(matrix)) <= 0.0:
                    continue
                region_factor_log = tile_log_scale - float(denominator_logs[occupation])
                inner_log = (
                    args.outer_bits * region_factor_log
                    + log_row_power_moment(matrix, args.outer_bits)
                )
                candidate = inner_log + distance * surprisal
                estimates_by_occupation.setdefault(int(occupation), []).append(candidate)
                denominators_by_occupation.setdefault(int(occupation), []).append(
                    float(denominator_logs[occupation])
                )
            print(
                f"tilt,{tilt_index + 1},{len(args.log_surprisals)},"
                f"tile,{center_index + 1},{len(centers)},center,{center}",
                flush=True,
            )

        for occupation, estimates in estimates_by_occupation.items():
            denominator_logs = denominators_by_occupation[occupation]
            preferred = int(np.argmax(np.asarray(denominator_logs)))
            candidate = estimates[preferred]
            evaluation_counts[occupation] += 1
            finite = [value for value in estimates if math.isfinite(value)]
            if len(finite) > 1:
                maximum_overlap_gap_bits = max(
                    maximum_overlap_gap_bits,
                    (max(finite) - min(finite)) / LOG2,
                )
            if candidate < best[occupation]:
                best[occupation] = candidate
                best_tilt[occupation] = log_surprisal
                best_center[occupation] = centers[preferred]

    rows = []
    for occupation in range(args.minimum_active_blocks, region_bits + 1):
        if not math.isfinite(float(best[occupation])):
            continue
        outer_log = (
            log_choose(outer_blocks, occupation)
            + occupation * regular_log_mass
        )
        inner_log = min(0.0, float(best[occupation]))
        rows.append(
            {
                "active_regular_outer_blocks": occupation,
                "best_log_surprisal": float(best_tilt[occupation]),
                "preferred_tile_center": int(best_center[occupation]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": (outer_log + inner_log) / LOG2,
                "tilt_evaluation_count": int(evaluation_counts[occupation]),
            }
        )
    dominant = sorted(
        rows, key=lambda row: float(row["pointwise_log2_upper"]), reverse=True
    )[:30]
    uncovered = [
        occupation
        for occupation in range(args.minimum_active_blocks, region_bits + 1)
        if not math.isfinite(float(best[occupation]))
    ]
    return {
        "schema": "riffle-fieldcheckpoint-regular-bulk-tiled-v1",
        "candidate": "Riffle FieldCheckpointAccumulate t=32 s=64 K=32",
        "method": {
            "candidate_count_tiling": "binomial exponential tilts",
            "region_composition": "three normalized FFT matrix-polynomial squares",
            "arithmetic": "nearest binary64 exploratory arithmetic",
            "outward_rounded": False,
            "maximum_overlap_gap_bits": maximum_overlap_gap_bits,
            "minimum_fft_coefficient_before_clamp": minimum_fft_coefficient,
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
            "epoch_bits": args.epoch_bits,
            "epochs_per_region": args.epochs_per_region,
            "minimum_active_blocks": args.minimum_active_blocks,
            "log_surprisals": args.log_surprisals,
            "minimum_log_denominator": args.minimum_log_denominator,
            "tile_count": len(centers),
        },
        "covered_occupation_count": len(rows),
        "uncovered_occupations": uncovered,
        "dominant_rows": dominant,
        "occupation_rows": rows,
        "scope": (
            "Exploratory regular-occupation bulk diagnostic. FFT rounding and "
            "clamping are not one-sided, so these values cannot support a "
            "certificate without independent outward-rounded replacement."
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
    parser.add_argument("--minimum-active-blocks", type=int, default=129)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-5.0, -4.0, -3.0, -2.0, -1.0, 0.0),
    )
    parser.add_argument("--minimum-log-denominator", type=float, default=-18.0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"covered_occupations,{payload['covered_occupation_count']}")
    if payload["dominant_rows"]:
        print(
            "bulk_dominant_pointwise_log2,"
            f"{payload['dominant_rows'][0]['pointwise_log2_upper']:.12f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
