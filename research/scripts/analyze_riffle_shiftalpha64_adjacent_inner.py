#!/usr/bin/env python3
"""Compute rigorous one-word inner bounds for two adjacent BCH profiles."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_shiftalpha64_compressed_return as compressed  # noqa: E402
from analyze_riffle_shiftalpha64_occ3_no_endpoint import (  # noqa: E402
    Profile,
    optimize_family,
)


PROFILES = ((22, 106, 106), (22, 22, 24))
PACKET_BITS = 4
PACKETS_PER_BLOCK = 32
STATE_COUNT = 16


def log_multinomial_packet_labels(block_count: int) -> float:
    return float(
        gammaln(block_count * PACKETS_PER_BLOCK + 1)
        - block_count * gammaln(PACKETS_PER_BLOCK + 1)
    )


def grouped_profile(weights: tuple[int, int, int]) -> list[tuple[int, int]]:
    groups: list[tuple[int, int]] = []
    for weight in weights:
        for index, (multiplicity, existing) in enumerate(groups):
            if existing == weight:
                groups[index] = (multiplicity + 1, existing)
                break
        else:
            groups.append((1, weight))
    return groups


def source_preserving_path_counts(
    state_factor: np.ndarray,
    groups: list[tuple[int, int]],
    source_tilts: np.ndarray,
    weight_tilts: np.ndarray,
) -> np.ndarray:
    """Tilt source packet counts and source-local binary weights separately."""
    packet_count = PACKETS_PER_BLOCK * sum(count for count, _ in groups)
    transition = np.zeros(STATE_COUNT, dtype=np.float64)
    for value in range(STATE_COUNT):
        value_weight = value.bit_count()
        transition[value] = sum(
            count * source_tilt * weight_tilt**value_weight
            for (count, _), source_tilt, weight_tilt in zip(
                groups, source_tilts, weight_tilts
            )
        )

    current = np.zeros((packet_count + 1, STATE_COUNT), dtype=np.float64)
    current[0, 0] = 1.0
    states = np.arange(STATE_COUNT)
    log_scale = 0.0
    for slot in range(packet_count):
        following = np.zeros_like(current)
        support_stop = slot + 1
        for value in range(STATE_COUNT):
            increment = int(value != 0)
            source = current[:support_stop, states ^ value]
            if value:
                source = source * state_factor[states]
            following[
                increment : increment + support_stop, states
            ] += source * transition[value]
        scale = float(np.max(following))
        if not math.isfinite(scale) or scale <= 0.0:
            return np.full(packet_count + 1, math.inf)
        following /= scale
        log_scale += math.log(scale)
        current = following

    totals = np.sum(current, axis=1)
    result = np.full(packet_count + 1, -math.inf)
    positive = totals > 0.0
    result[positive] = np.log(totals[positive]) + log_scale
    return result


def evaluate_source_preserving(
    weights: tuple[int, int, int],
    *,
    distance: int,
    packet_positions: int,
    point: np.ndarray,
) -> dict[str, object]:
    groups = grouped_profile(weights)
    scaled_cost = math.exp(float(point[0]))
    zero_scale = math.exp(float(point[1]))
    if not 0.0 < zero_scale < packet_positions:
        return {"raw_log_bound": math.inf}
    source_tilts = np.ones(len(groups), dtype=np.float64)
    if len(groups) > 1:
        source_tilts[:-1] = np.exp(point[2 : 1 + len(groups)])
    weight_offset = 1 + len(groups)
    weight_tilts = np.exp(point[weight_offset : weight_offset + len(groups)])

    u = scaled_cost / distance
    theta = zero_scale / packet_positions
    state_factor = np.empty(STATE_COUNT, dtype=np.float64)
    state_factor[0] = 1.0 / theta
    for state in range(1, STATE_COUNT):
        exponent = u * state.bit_count()
        state_factor[state] = math.exp(-exponent) / -math.expm1(-exponent)
    paths = source_preserving_path_counts(
        state_factor, groups, source_tilts, weight_tilts
    )

    coefficient_cost = 0.0
    normalization = log_multinomial_packet_labels(len(weights))
    for (multiplicity, weight), source_tilt, weight_tilt in zip(
        groups, source_tilts, weight_tilts
    ):
        coefficient_cost -= (
            multiplicity * PACKETS_PER_BLOCK * math.log(source_tilt)
            + multiplicity * weight * math.log(weight_tilt)
        )
        normalization += multiplicity * compressed.kernel.log_binom(128, weight)

    zero_gap = -math.log1p(-theta)
    terms = []
    for support, log_count in enumerate(paths):
        if math.isfinite(log_count):
            terms.append(
                (
                    support,
                    log_count
                    + coefficient_cost
                    - normalization
                    + (packet_positions - support + 1) * zero_gap
                    - compressed.kernel.log_binom(packet_positions, support),
                )
            )
    raw = scaled_cost + float(logsumexp([term for _, term in terms]))
    return {
        "raw_log_bound": raw,
        "log2_bound": min(0.0, raw) / math.log(2.0),
        "scaled_positive_cost": scaled_cost,
        "scaled_zero_return_parameter": zero_scale,
        "source_packet_count_tilts": source_tilts.tolist(),
        "source_binary_weight_tilts": weight_tilts.tolist(),
        "grouped_sources": [
            {"multiplicity": count, "binary_weight": weight}
            for count, weight in groups
        ],
        "top_packet_support": max(terms, key=lambda item: item[1])[0],
    }


def optimize_source_preserving(
    weights: tuple[int, int, int],
    *,
    distance: int,
    packet_positions: int,
    maxiter: int,
) -> dict[str, object]:
    groups = grouped_profile(weights)
    initial_source_logs = np.zeros(max(0, len(groups) - 1))
    initial_weight_logs = np.asarray(
        [math.log(weight / (128 - weight)) for _, weight in groups]
    )
    start = np.concatenate(
        [np.log(np.asarray([40.0, 18.0])), initial_source_logs, initial_weight_logs]
    )

    def objective(point: np.ndarray) -> float:
        return float(
            evaluate_source_preserving(
                weights,
                distance=distance,
                packet_positions=packet_positions,
                point=point,
            )["raw_log_bound"]
        )

    result = minimize(
        objective,
        start,
        method="Nelder-Mead",
        options={"xatol": 2e-5, "fatol": 2e-7, "maxiter": maxiter},
    )
    selected = evaluate_source_preserving(
        weights,
        distance=distance,
        packet_positions=packet_positions,
        point=result.x,
    )
    source_tilts = np.asarray(
        selected["source_packet_count_tilts"], dtype=np.float64
    )
    weight_tilts = np.asarray(
        selected["source_binary_weight_tilts"], dtype=np.float64
    )
    neutral_paths = source_preserving_path_counts(
        np.ones(STATE_COUNT, dtype=np.float64),
        groups,
        source_tilts,
        weight_tilts,
    )
    packet_count = PACKETS_PER_BLOCK * sum(count for count, _ in groups)
    expected_neutral_log = packet_count * math.log(
        sum(
            count * source_tilt * (1.0 + weight_tilt) ** PACKET_BITS
            for (count, _), source_tilt, weight_tilt in zip(
                groups, source_tilts, weight_tilts
            )
        )
    )
    neutral_error = abs(float(logsumexp(neutral_paths)) - expected_neutral_log)
    if neutral_error > 2e-10:
        raise RuntimeError("source-preserving neutral-mass identity failed")
    selected.update(
        {
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "optimizer_evaluations": int(result.nfev),
            "certificate": (
                "multivariate positive-coefficient envelope retaining source "
                "packet counts and each source-group binary weight"
            ),
            "neutral_mass_log_error": neutral_error,
            "neutral_mass_identity": "PASS",
        }
    )
    return selected


def analyze_profile(
    weights: tuple[int, int, int],
    *,
    distance: int,
    packet_positions: int,
    maxiter: int,
) -> dict[str, object]:
    profile = Profile(
        weights=weights,
        total_weight=sum(weights),
        log_count_bound=0.0,
        log_slice_normalization=sum(
            compressed.kernel.log_binom(128, weight) for weight in weights
        ),
    )
    result = optimize_family(
        profiles=[profile],
        support_count=1,
        distance=distance,
        packet_positions=packet_positions,
        maxiter=maxiter,
    )
    return {
        "weights": list(weights),
        "total_binary_input_weight": sum(weights),
        "pooled_one_word_log2_bound": min(0.0, result["log2_bound"]),
        "pooled_optimizer": result,
        "source_preserving": optimize_source_preserving(
            weights,
            distance=distance,
            packet_positions=packet_positions,
            maxiter=maxiter,
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance", type=int, default=188_766)
    parser.add_argument(
        "--packet-positions",
        type=int,
        default=compressed.DEFAULT_PACKET_POSITIONS,
    )
    parser.add_argument("--optimizer-maxiter", type=int, default=150)
    args = parser.parse_args()

    results = [
        analyze_profile(
            weights,
            distance=args.distance,
            packet_positions=args.packet_positions,
            maxiter=args.optimizer_maxiter,
        )
        for weights in PROFILES
    ]
    payload = {
        "schema": "riffle-shiftalpha64-adjacent-inner-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "parameters": {
            "bad_output_weight_inclusive": args.distance,
            "global_packet_positions": args.packet_positions,
        },
        "profiles": results,
        "validation": {
            "uniform_slice_normalization_gate": "PASS",
            "packet_support_marginal_gate": "PASS",
            "each_outer_word_counted_once": True,
            "source_packet_counts_retained_by_coefficient_tilts": True,
            "source_binary_weights_retained_by_coefficient_tilts": True,
        },
        "scope": (
            "Rigorous inner bad-output probability bounds for one fixed outer "
            "word of each declared BCH-weight profile. Outer multiplicity is "
            "not included."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
