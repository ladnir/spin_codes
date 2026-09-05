#!/usr/bin/env python3
"""Bound one outer occupation shell with three BCH-weight buckets.

The buckets are light (22--50), normal (52--76), and heavy (78--128).
For one bucket pattern, choose any two occupied outer coordinates as dependent
coordinates.  The MDS property makes their two columns invertible.  Summing
over the other coordinates and maximizing over the two determined values gives
a rigorous weighted count for the pattern.

Each pattern uses one optimized coefficient tilt.  The packet calculation
retains only packet support and the 16-state accumulator state.  The script
prints one JSON receipt and does not modify the workspace.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

import analyze_riffle_shiftalpha64_compressed_return as compressed
import analyze_riffle_shiftalpha64_outer_occupation as outer


BANDS = (
    ("light", 22, 50),
    ("normal", 52, 76),
    ("heavy", 78, 128),
)


def compositions(total: int) -> list[tuple[int, int, int]]:
    return [
        (light, normal, total - light - normal)
        for light in range(total + 1)
        for normal in range(total - light + 1)
    ]


def log_multinomial(counts: tuple[int, int, int]) -> float:
    total = sum(counts)
    return math.lgamma(total + 1) - sum(math.lgamma(count + 1) for count in counts)


def band_statistics(
    spectrum: dict[int, int], weight_tilt: float
) -> list[dict[str, float | int | str]]:
    log_x = math.log(weight_tilt)
    result = []
    for name, lower, upper in BANDS:
        values = [
            (
                weight,
                count,
                -compressed.kernel.log_binom(compressed.PART_BITS, weight)
                - weight * log_x,
            )
            for weight, count in spectrum.items()
            if lower <= weight <= upper and count > 0
        ]
        if not values:
            raise RuntimeError(f"empty BCH band: {name}")
        log_mass = float(
            logsumexp([math.log(count) + log_value for _, count, log_value in values])
        )
        max_weight, _, log_maximum = max(values, key=lambda item: item[2])
        result.append(
            {
                "name": name,
                "lower_weight": lower,
                "upper_weight": upper,
                "field_value_count": sum(count for _, count, _ in values),
                "log_mass": log_mass,
                "log_maximum": log_maximum,
                "maximum_weight": max_weight,
            }
        )
    return result


def best_elimination(
    counts: tuple[int, int, int], statistics: list[dict[str, float | int | str]]
) -> tuple[float, tuple[int, int, int]]:
    """Eliminate the best two coordinates using the MDS equations."""
    best_log = math.inf
    best_dependent = None
    for dependent in itertools.product(range(3), repeat=2):
        dependent_counts = tuple(dependent.count(index) for index in range(3))
        if any(dependent_counts[index] > counts[index] for index in range(3)):
            continue
        candidate = 0.0
        for index in range(3):
            free_count = counts[index] - dependent_counts[index]
            candidate += free_count * float(statistics[index]["log_mass"])
            candidate += dependent_counts[index] * float(
                statistics[index]["log_maximum"]
            )
        if candidate < best_log:
            best_log = candidate
            best_dependent = dependent_counts
    if best_dependent is None:
        raise RuntimeError("no valid pair of dependent coordinates")
    return best_log, best_dependent


def evaluate_pattern(
    *,
    counts: tuple[int, int, int],
    distance: int,
    packet_positions: int,
    spectrum: dict[int, int],
    scaled_cost: float,
    zero_scale: float,
    weight_tilt: float,
) -> dict[str, object]:
    occupation = sum(counts)
    if (
        scaled_cost <= 0.0
        or not 0.0 < zero_scale < packet_positions
        or weight_tilt <= 0.0
    ):
        return {"raw_log_bound": math.inf}

    compressed.configure_kernel(occupation, 22)
    u = scaled_cost / distance
    theta = zero_scale / packet_positions
    state_factor = np.empty(compressed.STATE_COUNT, dtype=np.float64)
    state_factor[0] = 1.0 / theta
    for state in range(1, compressed.STATE_COUNT):
        exponent = u * state.bit_count()
        state_factor[state] = math.exp(-exponent) / -math.expm1(-exponent)

    log_path_counts = compressed.weighted_path_log_counts_with_weight_tilt(
        state_factor, weight_tilt
    )
    statistics = band_statistics(spectrum, weight_tilt)
    log_value_sum, dependent_counts = best_elimination(counts, statistics)
    log_outer_factor = (
        compressed.kernel.log_binom(outer.OUTER_LENGTH, occupation)
        + log_multinomial(counts)
        + log_value_sum
    )
    zero_gap_log_factor = -math.log1p(-theta)
    terms = []
    for support, log_count in enumerate(log_path_counts):
        if not math.isfinite(log_count):
            continue
        term = (
            log_count
            + log_outer_factor
            + (packet_positions - support + 1) * zero_gap_log_factor
            - compressed.kernel.log_binom(packet_positions, support)
        )
        terms.append((support, term))

    raw_log_bound = scaled_cost + float(logsumexp([term for _, term in terms]))
    return {
        "raw_log_bound": raw_log_bound,
        "log2_bound": raw_log_bound / math.log(2.0),
        "band_counts": {
            BANDS[index][0]: counts[index] for index in range(3)
        },
        "dependent_band_counts": {
            BANDS[index][0]: dependent_counts[index] for index in range(3)
        },
        "scaled_positive_cost": scaled_cost,
        "scaled_zero_return_parameter": zero_scale,
        "binary_weight_coefficient_tilt": weight_tilt,
        "band_statistics": [
            {
                "name": item["name"],
                "field_value_count": item["field_value_count"],
                "maximum_weight_for_dependent_coordinate": item["maximum_weight"],
                "log2_free_coordinate_mass": float(item["log_mass"]) / math.log(2.0),
                "log2_dependent_coordinate_maximum": float(item["log_maximum"])
                / math.log(2.0),
            }
            for item in statistics
        ],
        "top_packet_support": max(terms, key=lambda item: item[1])[0],
    }


def optimize_pattern(
    *,
    counts: tuple[int, int, int],
    distance: int,
    packet_positions: int,
    spectrum: dict[int, int],
    maxiter: int,
) -> dict[str, object]:
    occupation = sum(counts)
    representatives = (40.0, 64.0, 90.0)
    average_weight = sum(
        counts[index] * representatives[index] for index in range(3)
    ) / occupation
    start_weight_tilt = min(4.0, max(0.25, average_weight / (128.0 - average_weight)))
    start = np.log(
        np.asarray([22.0 * occupation, 8.0 * occupation, start_weight_tilt])
    )

    def objective(log_point: np.ndarray) -> float:
        scaled_cost, zero_scale, weight_tilt = np.exp(log_point)
        return float(
            evaluate_pattern(
                counts=counts,
                distance=distance,
                packet_positions=packet_positions,
                spectrum=spectrum,
                scaled_cost=float(scaled_cost),
                zero_scale=float(zero_scale),
                weight_tilt=float(weight_tilt),
            )["raw_log_bound"]
        )

    result = minimize(
        objective,
        start,
        method="Nelder-Mead",
        options={"xatol": 3e-5, "fatol": 3e-7, "maxiter": maxiter},
    )
    scaled_cost, zero_scale, weight_tilt = np.exp(result.x)
    selected = evaluate_pattern(
        counts=counts,
        distance=distance,
        packet_positions=packet_positions,
        spectrum=spectrum,
        scaled_cost=float(scaled_cost),
        zero_scale=float(zero_scale),
        weight_tilt=float(weight_tilt),
    )
    selected.update(
        {
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "optimizer_evaluations": int(result.nfev),
        }
    )
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--occupation", type=int, choices=(3, 4), required=True)
    parser.add_argument("--distance", type=int, default=188_766)
    parser.add_argument(
        "--packet-positions",
        type=int,
        default=compressed.DEFAULT_PACKET_POSITIONS,
    )
    parser.add_argument("--spectrum", type=Path, default=outer.DEFAULT_SPECTRUM)
    parser.add_argument(
        "--exclude-all-one",
        action="store_true",
        help="exclude the unique BCH word of weight 128",
    )
    parser.add_argument("--optimizer-maxiter", type=int, default=60)
    args = parser.parse_args()

    spectrum = outer.load_spectrum(args.spectrum)
    if args.exclude_all_one:
        spectrum = dict(spectrum)
        if spectrum.get(128) != 1:
            raise RuntimeError("all-one BCH spectrum entry failed validation")
        spectrum[128] = 0
    pattern_results = []
    for counts in compositions(args.occupation):
        pattern_results.append(
            optimize_pattern(
                counts=counts,
                distance=args.distance,
                packet_positions=args.packet_positions,
                spectrum=spectrum,
                maxiter=args.optimizer_maxiter,
            )
        )

    total_raw_log = float(
        logsumexp([float(item["raw_log_bound"]) for item in pattern_results])
    )
    exact_mds_count = outer.exact_mds_occupation_count(args.occupation)
    payload = {
        "schema": "riffle-shiftalpha64-threeband-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "parameters": {
            "occupied_field_symbols": args.occupation,
            "bad_output_weight_inclusive": args.distance,
            "bands": [
                {"name": name, "minimum": lower, "maximum": upper}
                for name, lower, upper in BANDS
            ],
            "excluded_bch_weights": [128] if args.exclude_all_one else [],
        },
        "complete_shell": {
            "log2_bound": total_raw_log / math.log(2.0),
            "exact_mds_log2_multiplicity": math.log2(exact_mds_count),
            "log2_bound_after_trivial_message_count_cap": min(
                total_raw_log / math.log(2.0), math.log2(exact_mds_count)
            ),
            "beats_trivial_message_count": (
                total_raw_log / math.log(2.0) < math.log2(exact_mds_count)
            ),
            "pattern_count": len(pattern_results),
        },
        "patterns": sorted(
            pattern_results, key=lambda item: float(item["raw_log_bound"]), reverse=True
        ),
        "method": {
            "free_coordinates": "sum the exact BCH-spectrum mass inside each band",
            "dependent_coordinates": "maximize their band-restricted coefficient weight",
            "outer_checks": "choose two coordinates and solve them using the invertible MDS column pair",
            "inner_state": "packet support and current GF(2)^4 accumulator state",
        },
        "validation": {
            "exact_bch_spectrum": "PASS",
            "patterns_partition_every_included_nonzero_bch_weight": True,
            "optimizer_warning_does_not_invalidate_selected_parameters": True,
        },
        "scope": (
            "Rigorous upper bound for the included part of one outer occupation shell. "
            "The two solved values are bounded separately by their declared bands."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
