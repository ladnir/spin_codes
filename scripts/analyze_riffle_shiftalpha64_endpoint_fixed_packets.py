#!/usr/bin/env python3
"""Transfer the exact endpoint histogram while preserving the all-one block.

An endpoint triple contains one all-one BCH block and two complementary BCH
blocks.  The all-one block contributes 32 packets of value 15.  This script
retains those packets exactly.  It pools only the remaining 256 bits and
applies one coefficient tilt to their fixed total weight 128.
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


DEFAULT_SCAN = (
    Path(__file__).resolve().parent.parent
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "receipts"
    / "goal09_endpoint_triple_scan.json"
)
STATE_COUNT = 16
FIXED_VALUE = 15
FIXED_PACKETS = 32
VARIABLE_PACKETS = 64
VARIABLE_WEIGHT = 128
STAGING_PACKETS = FIXED_PACKETS + VARIABLE_PACKETS


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


def fixed_packet_path_log_counts(
    state_factor: np.ndarray,
    weight_tilt: float,
    *,
    fixed_packets: int = FIXED_PACKETS,
    variable_packets: int = VARIABLE_PACKETS,
) -> np.ndarray:
    """Sum paths by support with an exact number of fixed value-15 packets."""
    if weight_tilt <= 0.0:
        return np.full(fixed_packets + variable_packets + 1, math.inf)

    total_packets = fixed_packets + variable_packets
    shape = (fixed_packets + 1, total_packets + 1, STATE_COUNT)
    current = np.zeros(shape, dtype=np.float64)
    current[0, 0, 0] = 1.0
    states = np.arange(STATE_COUNT)
    log_scale = 0.0

    active_transition = np.zeros((STATE_COUNT, STATE_COUNT), dtype=np.float64)
    for previous in range(STATE_COUNT):
        for following in range(STATE_COUNT):
            value = previous ^ following
            if value:
                active_transition[previous, following] = (
                    weight_tilt ** value.bit_count() * state_factor[following]
                )

    for slot in range(total_packets):
        support_stop = slot + 1
        following = np.zeros_like(current)

        # A zero-valued variable packet is absent from the active walk.
        following[:, :support_stop, :] += current[:, :support_stop, :]

        # A nonzero variable packet advances the walk and the active support.
        following[:, 1 : support_stop + 1, :] += np.matmul(
            current[:, :support_stop, :], active_transition
        )

        # A fixed packet has value 15 and is always active.
        following[1:, 1 : support_stop + 1, :] += (
            current[:-1, :support_stop, states ^ FIXED_VALUE]
            * state_factor[states]
        )

        processed = slot + 1
        minimum_fixed = max(0, processed - variable_packets)
        maximum_fixed = min(fixed_packets, processed)
        if minimum_fixed:
            following[:minimum_fixed, :, :] = 0.0
        if maximum_fixed < fixed_packets:
            following[maximum_fixed + 1 :, :, :] = 0.0

        scale = float(np.max(following))
        if not math.isfinite(scale) or scale <= 0.0:
            return np.full(total_packets + 1, math.inf)
        following /= scale
        log_scale += math.log(scale)
        current = following

    totals = np.sum(current[fixed_packets, :, :], axis=1)
    result = np.full(total_packets + 1, -math.inf)
    positive = totals > 0.0
    result[positive] = np.log(totals[positive]) + log_scale
    return result


def evaluate(
    *,
    histogram: dict[int, int],
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
    state_factor = np.empty(STATE_COUNT, dtype=np.float64)
    state_factor[0] = 1.0 / theta
    for state in range(1, STATE_COUNT):
        exponent = u * state.bit_count()
        state_factor[state] = math.exp(-exponent) / -math.expm1(-exponent)

    log_path_counts = fixed_packet_path_log_counts(state_factor, weight_tilt)
    log_x = math.log(weight_tilt)
    profile_terms = []
    for weight, count in histogram.items():
        complement = 128 - weight
        term = (
            math.log(count)
            - compressed.kernel.log_binom(128, weight)
            - compressed.kernel.log_binom(128, complement)
            - VARIABLE_WEIGHT * log_x
        )
        profile_terms.append((weight, term))
    log_outer_factor = float(logsumexp([term for _, term in profile_terms]))

    interleaving_log_normalization = compressed.kernel.log_binom(
        STAGING_PACKETS, FIXED_PACKETS
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
    starts = (
        np.log(np.asarray([68.8, 24.9, 1.0])),
        np.log(np.asarray([90.0, 30.0, 0.8])),
    )
    results = []

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

    for start in starts:
        results.append(
            minimize(
                objective,
                start,
                method="Nelder-Mead",
                options={"xatol": 2e-5, "fatol": 2e-7, "maxiter": maxiter},
            )
        )
    result = min(results, key=lambda item: float(item.fun))
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
            "optimizer_evaluations_selected_start": int(result.nfev),
            "optimizer_total_evaluations": int(sum(item.nfev for item in results)),
        }
    )
    return selected


def direct_small_counts(
    state_factor: np.ndarray,
    weight_tilt: float,
    fixed_packets: int,
    variable_packets: int,
) -> np.ndarray:
    """Independently enumerate a small fixed/variable packet instance."""
    total_packets = fixed_packets + variable_packets
    totals = np.zeros(total_packets + 1, dtype=np.float64)
    for fixed_positions in itertools.combinations(range(total_packets), fixed_packets):
        fixed_set = set(fixed_positions)
        for variable_values in itertools.product(range(STATE_COUNT), repeat=variable_packets):
            iterator = iter(variable_values)
            state = 0
            support = 0
            contribution = 1.0
            for slot in range(total_packets):
                value = FIXED_VALUE if slot in fixed_set else next(iterator)
                if value:
                    support += 1
                    state ^= value
                    contribution *= (
                        state_factor[state] * weight_tilt ** value.bit_count()
                        if slot not in fixed_set
                        else state_factor[state]
                    )
            totals[support] += contribution
    return totals


def validate_operator() -> dict[str, str]:
    state_factor = np.asarray(
        [1.7] + [0.55 + 0.03 * state.bit_count() for state in range(1, STATE_COUNT)],
        dtype=np.float64,
    )
    weight_tilt = 0.73
    fixed_packets = 2
    variable_packets = 3
    observed_logs = fixed_packet_path_log_counts(
        state_factor,
        weight_tilt,
        fixed_packets=fixed_packets,
        variable_packets=variable_packets,
    )
    expected = direct_small_counts(
        state_factor, weight_tilt, fixed_packets, variable_packets
    )
    for support, exact in enumerate(expected):
        observed = 0.0 if not math.isfinite(observed_logs[support]) else math.exp(
            observed_logs[support]
        )
        if not math.isclose(observed, exact, rel_tol=2e-13, abs_tol=2e-13):
            raise RuntimeError(f"small direct audit failed at support {support}")

    neutral = fixed_packet_path_log_counts(
        np.ones(STATE_COUNT, dtype=np.float64), 1.0
    )
    observed_total_log = float(logsumexp(neutral))
    expected_total_log = (
        compressed.kernel.log_binom(STAGING_PACKETS, FIXED_PACKETS)
        + VARIABLE_PACKETS * math.log(STATE_COUNT)
    )
    if not math.isclose(
        observed_total_log, expected_total_log, rel_tol=0.0, abs_tol=2e-11
    ):
        raise RuntimeError("full neutral-mass audit failed")
    return {
        "small_direct_enumeration": "PASS (2 fixed, 3 variable packets)",
        "full_neutral_mass": "PASS (C(96,32) * 16^64)",
    }


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
    if args.distance <= 0:
        raise ValueError("distance must be positive")
    if args.packet_positions < STAGING_PACKETS:
        raise ValueError("packet count is smaller than the endpoint superblock")

    validation = validate_operator()
    histogram = load_histogram(args.scan)
    result = optimize(
        histogram, args.distance, args.packet_positions, args.optimizer_maxiter
    )
    payload = {
        "schema": "riffle-shiftalpha64-endpoint-fixed-packets-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "source_scan": str(args.scan),
        "endpoint_words_in_finite_support_scan": sum(histogram.values()),
        "preserved_structure": {
            "fixed_packets": FIXED_PACKETS,
            "fixed_packet_value": "1111",
            "variable_packets": VARIABLE_PACKETS,
            "variable_total_binary_weight": VARIABLE_WEIGHT,
            "random_interleaving_normalization": "C(96,32)",
        },
        "bound": result,
        "validation": {
            **validation,
            "exact_histogram_total": "PASS",
            "selected_parameters_valid_if_optimizer_stops": True,
        },
        "scope": (
            "Rigorous transfer for outer weight-three words on finite supports "
            "that contain the all-one BCH word. The two complementary blocks "
            "are pooled before their common total-weight coefficient bound."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
