#!/usr/bin/env python3
"""Transfer the exact finite-support all-one triple histogram through the inner kernel."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

import analyze_riffle_shiftalpha64_compressed_return as compressed


DEFAULT_SCAN = (
    Path(__file__).resolve().parent.parent
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "receipts"
    / "goal09_endpoint_triple_scan.json"
)


def load_histogram(path: Path) -> dict[int, int]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    histogram = {
        int(weight): int(count)
        for weight, count in payload["candidate_weight_histogram"].items()
    }
    expected = 3 * (
        int(payload["data_only_supports"])
        + int(payload["p0_and_two_data_supports"])
    )
    if sum(histogram.values()) != expected:
        raise RuntimeError("endpoint histogram failed total-count validation")
    if min(histogram) != 32 or max(histogram) != 96:
        raise RuntimeError("endpoint histogram failed range validation")
    return histogram


def evaluate(
    *,
    histogram: dict[int, int],
    distance: int,
    packet_positions: int,
    scaled_cost: float,
    zero_scale: float,
    weight_tilt: float,
) -> dict[str, object]:
    compressed.configure_kernel(3, 22)
    if (
        scaled_cost <= 0.0
        or not 0.0 < zero_scale < packet_positions
        or weight_tilt <= 0.0
    ):
        return {"raw_log_bound": math.inf}

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

    log_x = math.log(weight_tilt)
    profile_terms = []
    for weight, count in histogram.items():
        complement = 128 - weight
        term = (
            math.log(count)
            - compressed.kernel.log_binom(128, weight)
            - compressed.kernel.log_binom(128, complement)
            - 256 * log_x
        )
        profile_terms.append((weight, term))
    log_outer_factor = float(logsumexp([term for _, term in profile_terms]))

    zero_gap_log_factor = -math.log1p(-theta)
    support_terms = []
    for support, log_count in enumerate(log_path_counts):
        if not math.isfinite(log_count):
            continue
        term = (
            log_count
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
        "finite_endpoint_shell_log2_bound": raw_log_bound / math.log(2.0),
        "scaled_positive_cost": scaled_cost,
        "scaled_zero_return_parameter": zero_scale,
        "binary_weight_coefficient_tilt": weight_tilt,
        "top_profile_terms": [
            {"candidate_weight": weight, "log2_contribution": term / math.log(2.0)}
            for weight, term in sorted(
                profile_terms, key=lambda item: item[1], reverse=True
            )[:8]
        ],
        "top_packet_support": max(support_terms, key=lambda item: item[1])[0],
    }


def optimize(
    histogram: dict[int, int], distance: int, packet_positions: int, maxiter: int
) -> dict[str, object]:
    start = np.log(np.asarray([65.0, 24.0, 1.3]))

    def objective(point: np.ndarray) -> float:
        scaled_cost, zero_scale, weight_tilt = np.exp(point)
        return float(
            evaluate(
                histogram=histogram,
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
    selected = evaluate(
        histogram=histogram,
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
    parser.add_argument("--optimizer-maxiter", type=int, default=80)
    args = parser.parse_args()
    histogram = load_histogram(args.scan)
    result = optimize(
        histogram, args.distance, args.packet_positions, args.optimizer_maxiter
    )
    payload = {
        "schema": "riffle-shiftalpha64-endpoint-triple-transfer-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "source_scan": str(args.scan),
        "endpoint_words_in_finite_support_scan": sum(histogram.values()),
        "bound": result,
        "validation": {
            "exact_histogram_total": "PASS",
            "exact_candidate_weight_range": [min(histogram), max(histogram)],
            "selected_parameters_valid_if_optimizer_stops": True,
        },
        "scope": (
            "Rigorous transfer for outer weight-three words on finite supports "
            "that contain the all-one BCH word. Supports containing p1 are excluded."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
