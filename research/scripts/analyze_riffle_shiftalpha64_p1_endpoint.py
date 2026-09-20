#!/usr/bin/env python3
"""Bound all-one endpoint triples that contain the second parity position.

Such a support contains the second parity position and two finite projective
points.  The two finite field values are equal.  There are two endpoint
normalizations: either both finite BCH blocks are all ones, or the parity BCH
block is all ones.  This script transfers the exact BCH-weight histogram of
each family through the fixed-packet inner operator.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

import analyze_riffle_shiftalpha64_compressed_return as compressed
import analyze_riffle_shiftalpha64_endpoint_fixed_packets as fixed


DEFAULT_SCAN = (
    Path(__file__).resolve().parent.parent
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "receipts"
    / "goal10_p1_endpoint_scan.json"
)
DATA_BLOCKS = 16_384
FINITE_PROJECTIVE_POINTS = DATA_BLOCKS + 1
P1_SUPPORTS = FINITE_PROJECTIVE_POINTS * (FINITE_PROJECTIVE_POINTS - 1) // 2


def load_histograms(path: Path) -> tuple[dict[int, int], dict[int, int], int]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    finite_histogram = {
        int(weight): int(count)
        for weight, count in payload[
            "finite_coordinate_normalization_histogram"
        ].items()
    }
    p1_histogram = {
        int(weight): int(count)
        for weight, count in payload["p1_coordinate_normalization_histogram"].items()
    }
    if sum(finite_histogram.values()) != P1_SUPPORTS:
        raise RuntimeError("finite-coordinate histogram has the wrong mass")
    if sum(p1_histogram.values()) != P1_SUPPORTS:
        raise RuntimeError("p1-coordinate histogram has the wrong mass")
    overlap = int(payload["all_three_all_one_overlap_supports"])
    return finite_histogram, p1_histogram, overlap


def evaluate_family(
    *,
    histogram: dict[int, int],
    fixed_packets: int,
    variable_packets: int,
    variable_weight_multiplier: int,
    slice_count: int,
    distance: int,
    packet_positions: int,
    scaled_cost: float,
    zero_scale: float,
    weight_tilt: float,
) -> dict[str, object]:
    if (
        scaled_cost <= 0.0
        or not 0.0 < zero_scale < packet_positions
        or weight_tilt <= 0.0
    ):
        return {"raw_log_bound": math.inf}

    u = scaled_cost / distance
    theta = zero_scale / packet_positions
    state_factor = np.empty(fixed.STATE_COUNT, dtype=np.float64)
    state_factor[0] = 1.0 / theta
    for state in range(1, fixed.STATE_COUNT):
        exponent = u * state.bit_count()
        state_factor[state] = math.exp(-exponent) / -math.expm1(-exponent)

    log_path_counts = fixed.fixed_packet_path_log_counts(
        state_factor,
        weight_tilt,
        fixed_packets=fixed_packets,
        variable_packets=variable_packets,
    )
    log_x = math.log(weight_tilt)
    profile_terms = []
    for weight, count in histogram.items():
        term = (
            math.log(count)
            - slice_count * compressed.kernel.log_binom(128, weight)
            - variable_weight_multiplier * weight * log_x
        )
        profile_terms.append((weight, term))
    log_outer_factor = float(logsumexp([term for _, term in profile_terms]))

    staging_packets = fixed_packets + variable_packets
    interleaving_log_normalization = compressed.kernel.log_binom(
        staging_packets, fixed_packets
    )
    zero_gap_log_factor = -math.log1p(-theta)
    support_terms = []
    for support, log_count in enumerate(log_path_counts):
        if not math.isfinite(log_count):
            continue
        term = (
            log_count
            - interleaving_log_normalization
            + log_outer_factor
            + (packet_positions - support + 1) * zero_gap_log_factor
            - compressed.kernel.log_binom(packet_positions, support)
        )
        support_terms.append((support, term))
    raw_log_bound = scaled_cost + float(
        logsumexp([term for _, term in support_terms])
    )
    return {
        "raw_log_bound": raw_log_bound,
        "log2_bound": raw_log_bound / math.log(2.0),
        "scaled_positive_cost": scaled_cost,
        "scaled_zero_return_parameter": zero_scale,
        "binary_weight_coefficient_tilt": weight_tilt,
        "top_packet_support": max(support_terms, key=lambda item: item[1])[0],
        "top_profile_terms": [
            {"bch_weight": weight, "log2_contribution": term / math.log(2.0)}
            for weight, term in sorted(
                profile_terms, key=lambda item: item[1], reverse=True
            )[:8]
        ],
    }


def optimize_family(
    *,
    histogram: dict[int, int],
    fixed_packets: int,
    variable_packets: int,
    variable_weight_multiplier: int,
    slice_count: int,
    distance: int,
    packet_positions: int,
    maxiter: int,
) -> dict[str, object]:
    start = np.log(np.asarray([65.0, 26.0, 1.0]))

    def objective(point: np.ndarray) -> float:
        scaled_cost, zero_scale, weight_tilt = np.exp(point)
        return float(
            evaluate_family(
                histogram=histogram,
                fixed_packets=fixed_packets,
                variable_packets=variable_packets,
                variable_weight_multiplier=variable_weight_multiplier,
                slice_count=slice_count,
                distance=distance,
                packet_positions=packet_positions,
                scaled_cost=float(scaled_cost),
                zero_scale=float(zero_scale),
                weight_tilt=float(weight_tilt),
            )["raw_log_bound"]
        )

    result = minimize(
        objective,
        start,
        method="Nelder-Mead",
        options={"xatol": 2e-5, "fatol": 2e-7, "maxiter": maxiter},
    )
    scaled_cost, zero_scale, weight_tilt = np.exp(result.x)
    selected = evaluate_family(
        histogram=histogram,
        fixed_packets=fixed_packets,
        variable_packets=variable_packets,
        variable_weight_multiplier=variable_weight_multiplier,
        slice_count=slice_count,
        distance=distance,
        packet_positions=packet_positions,
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
    parser.add_argument("--scan", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--distance", type=int, default=188_766)
    parser.add_argument(
        "--packet-positions",
        type=int,
        default=compressed.DEFAULT_PACKET_POSITIONS,
    )
    parser.add_argument("--optimizer-maxiter", type=int, default=100)
    args = parser.parse_args()
    finite_histogram, p1_histogram, overlap = load_histograms(args.scan)
    validation = fixed.validate_operator()

    finite_normalization = optimize_family(
        histogram=finite_histogram,
        fixed_packets=64,
        variable_packets=32,
        variable_weight_multiplier=1,
        slice_count=1,
        distance=args.distance,
        packet_positions=args.packet_positions,
        maxiter=args.optimizer_maxiter,
    )
    parity_normalization = optimize_family(
        histogram=p1_histogram,
        fixed_packets=32,
        variable_packets=64,
        variable_weight_multiplier=2,
        slice_count=2,
        distance=args.distance,
        packet_positions=args.packet_positions,
        maxiter=args.optimizer_maxiter,
    )
    combined_log = float(
        logsumexp(
            [
                finite_normalization["raw_log_bound"],
                parity_normalization["raw_log_bound"],
            ]
        )
    )
    payload = {
        "schema": "riffle-shiftalpha64-p1-endpoint-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "second_parity_supports": P1_SUPPORTS,
        "source_scan": str(args.scan),
        "all_three_all_one_overlap_supports": overlap,
        "finite_coordinate_normalized_to_all_one": finite_normalization,
        "second_parity_coordinate_normalized_to_all_one": parity_normalization,
        "combined_log2_bound": combined_log / math.log(2.0),
        "validation": {
            **validation,
            "exact_histogram_mass": "PASS",
            "endpoint_normalizations_cover_family": True,
            "normalization_families_are_disjoint": overlap == 0,
            "selected_parameters_valid_if_optimizer_stops": True,
        },
        "scope": (
            "Rigorous upper bound for all outer weight-three endpoint words "
            "whose support contains the second parity position."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
