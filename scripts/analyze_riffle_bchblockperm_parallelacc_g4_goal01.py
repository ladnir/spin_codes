#!/usr/bin/env python3
"""Analyze one permuted BCH block before the global packet accumulator.

The block permutation maps every BCH word of binary weight h to a uniform
weight-h vector in 128 coordinates.  The script averages the joint four-lane
gap moment over that exact slice and then optimizes its two Chernoff
parameters.  It prints a JSON receipt and does not modify the workspace.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp


BLOCK_BITS = 128
PACKET_BITS = 4
PACKETS_PER_BLOCK = BLOCK_BITS // PACKET_BITS
STATE_COUNT = 1 << PACKET_BITS


def log_binom(n: int, k: int) -> float:
    if not 0 <= k <= n:
        return -math.inf
    return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))


def packet_support_count(input_weight: int, support: int) -> int:
    """Count weight-h block vectors that occupy exactly H packets."""
    coefficient = 0
    for selected in range(support + 1):
        available = PACKET_BITS * selected
        if available < input_weight:
            continue
        coefficient += (
            (-1) ** (support - selected)
            * math.comb(support, selected)
            * math.comb(available, input_weight)
        )
    return math.comb(PACKETS_PER_BLOCK, support) * coefficient


def packet_support_distribution(input_weight: int) -> dict[int, float]:
    denominator = math.comb(BLOCK_BITS, input_weight)
    return {
        support: count / denominator
        for support in range(1, PACKETS_PER_BLOCK + 1)
        if (count := packet_support_count(input_weight, support))
    }


def weighted_path_counts(input_weight: int, state_factor: np.ndarray) -> np.ndarray:
    """Sum a state-path product over every weight-h vector, separated by H."""
    shape = (input_weight + 1, PACKETS_PER_BLOCK + 1, STATE_COUNT)
    current = np.zeros(shape, dtype=np.float64)
    current[0, 0, 0] = 1.0
    states = np.arange(STATE_COUNT)

    for _slot in range(PACKETS_PER_BLOCK):
        following = np.zeros_like(current)
        for value in range(STATE_COUNT):
            weight = value.bit_count()
            support_increment = int(value != 0)
            if weight > input_weight:
                continue
            source = current[
                : input_weight + 1 - weight,
                : PACKETS_PER_BLOCK + 1 - support_increment,
                states ^ value,
            ]
            if value:
                source = source * state_factor[states]
            following[
                weight:,
                support_increment:,
                states,
            ] += source
        current = following

    return np.sum(current[input_weight, :, :], axis=1)


def zero_return_statistics(
    input_weight: int, distance: int, packet_positions: int
) -> dict[str, object]:
    """Count zero prefix states exactly over the uniform block slice."""
    max_returns = PACKETS_PER_BLOCK // 2
    shape = (
        input_weight + 1,
        PACKETS_PER_BLOCK + 1,
        STATE_COUNT,
        max_returns + 1,
    )
    current = np.zeros(shape, dtype=np.float64)
    current[0, 0, 0, 0] = 1.0

    for _slot in range(PACKETS_PER_BLOCK):
        following = np.zeros_like(current)
        for value in range(STATE_COUNT):
            weight = value.bit_count()
            support_increment = int(value != 0)
            if weight > input_weight:
                continue
            for next_state in range(STATE_COUNT):
                return_increment = int(value != 0 and next_state == 0)
                following[
                    weight:,
                    support_increment:,
                    next_state,
                    return_increment:,
                ] += current[
                    : input_weight + 1 - weight,
                    : PACKETS_PER_BLOCK + 1 - support_increment,
                    next_state ^ value,
                    : max_returns + 1 - return_increment,
                ]
        current = following

    denominator = math.comb(BLOCK_BITS, input_weight)
    return_counts = np.sum(current[input_weight, :, :, :], axis=(0, 1))
    probabilities = return_counts / denominator
    terminal_zero = float(
        np.sum(current[input_weight, :, 0, :]) / denominator
    )
    joint_counts = np.sum(current[input_weight, :, :, :], axis=1)
    exact_zero_gap_bound = 0.0
    contributions = []
    for support in range(1, PACKETS_PER_BLOCK + 1):
        for zero_returns in range(max_returns + 1):
            slice_probability = float(joint_counts[support, zero_returns] / denominator)
            if slice_probability <= 0.0:
                continue
            positive_gaps = support - zero_returns
            if distance < positive_gaps:
                conditional_bound = 0.0
            else:
                log_bound = (
                    log_binom(distance, positive_gaps)
                    + log_binom(
                        packet_positions - support + zero_returns,
                        zero_returns,
                    )
                    - log_binom(packet_positions, support)
                )
                conditional_bound = min(1.0, math.exp(log_bound))
            contribution = slice_probability * conditional_bound
            exact_zero_gap_bound += contribution
            contributions.append(
                {
                    "packet_support": support,
                    "zero_prefix_states": zero_returns,
                    "slice_probability": slice_probability,
                    "conditional_gap_bound": conditional_bound,
                    "weighted_contribution": contribution,
                }
            )
    return {
        "terminal_zero_probability": terminal_zero,
        "mean_zero_prefix_states": float(
            sum(index * probability for index, probability in enumerate(probabilities))
        ),
        "probability_any_zero_prefix_state": float(1.0 - probabilities[0]),
        "distribution": [
            [index, float(probability)]
            for index, probability in enumerate(probabilities)
            if probability > 0.0
        ],
        "exact_zero_gap_log2_bound": (
            math.log2(exact_zero_gap_bound) if exact_zero_gap_bound else -math.inf
        ),
        "exact_zero_gap_effective_root": exact_zero_gap_bound
        ** (1.0 / input_weight),
        "top_exact_zero_gap_contributions": sorted(
            contributions,
            key=lambda item: item["weighted_contribution"],
            reverse=True,
        )[:12],
    }


def expected_path_products(
    input_weight: int, u: float, scaled_v: float, packet_positions: int
) -> np.ndarray:
    """Average the ordered-state product for every nonzero-packet support H."""
    v_gap = scaled_v / packet_positions
    state_factor = np.empty(STATE_COUNT, dtype=np.float64)
    for state in range(STATE_COUNT):
        weight = state.bit_count()
        exponent = v_gap + weight * u
        denominator = -math.expm1(-exponent)
        state_factor[state] = math.exp(-weight * u) / denominator

    normalization = math.comb(BLOCK_BITS, input_weight)
    return weighted_path_counts(input_weight, state_factor) / normalization


def optimize_block_moment(
    input_weight: int, distance: int, packet_positions: int
) -> dict[str, object]:
    def objective(point: np.ndarray) -> float:
        u = float(point[0])
        scaled_v = float(point[1])
        averages = expected_path_products(
            input_weight, u, scaled_v, packet_positions
        )
        v_gap = scaled_v / packet_positions
        terms = []
        for support, average in enumerate(averages):
            if average <= 0.0:
                continue
            terms.append(
                (packet_positions - support) * v_gap
                - log_binom(packet_positions, support)
                + math.log(average)
            )
        if not terms:
            return math.inf
        value = distance * u + float(logsumexp(terms))
        return value if math.isfinite(value) else 1e300

    start_u = max(1e-3, math.log(packet_positions / max(1, distance)) / 4.0)
    result = minimize(
        objective,
        np.asarray([start_u, float(max(1, input_weight))]),
        method="Nelder-Mead",
        bounds=((1e-9, 30.0), (1e-6, float(packet_positions))),
        options={"xatol": 2e-8, "fatol": 2e-10, "maxiter": 600},
    )
    u = float(result.x[0])
    scaled_v = float(result.x[1])
    log_bound = min(0.0, objective(np.asarray([u, scaled_v])))
    return {
        "log2_bound": log_bound / math.log(2.0),
        "effective_root_per_input_bit": math.exp(log_bound / input_weight),
        "z": math.exp(-u),
        "tau": math.exp(-scaled_v / packet_positions),
        "scaled_gap_parameter": scaled_v,
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "optimizer_evaluations": int(result.nfev),
    }


def scalar_accumulator_bad_probability(
    binary_length: int, input_weight: int, distance: int
) -> float:
    numerator = 0
    left = input_weight // 2
    right = (input_weight + 1) // 2 - 1
    for output_weight in range(1, distance + 1):
        if left <= binary_length - output_weight and right <= output_weight - 1:
            numerator += math.comb(binary_length - output_weight, left) * math.comb(
                output_weight - 1, right
            )
    return numerator / math.comb(binary_length, input_weight)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-positions", type=int, default=524_352)
    parser.add_argument("--input-weight", type=int, default=22)
    parser.add_argument("--distance", type=int, default=40)
    args = parser.parse_args()

    if not 1 <= args.input_weight <= BLOCK_BITS:
        raise ValueError("input weight must lie in [1,128]")
    if args.packet_positions < PACKETS_PER_BLOCK:
        raise ValueError("global packet count must contain one BCH block")

    support_distribution = packet_support_distribution(args.input_weight)
    support_total = sum(support_distribution.values())
    if not math.isclose(support_total, 1.0, rel_tol=2e-13, abs_tol=2e-13):
        raise RuntimeError("packet-support distribution failed normalization")
    unweighted_counts = weighted_path_counts(
        args.input_weight, np.ones(STATE_COUNT, dtype=np.float64)
    )
    denominator = math.comb(BLOCK_BITS, args.input_weight)
    for support, probability in support_distribution.items():
        dynamic_probability = unweighted_counts[support] / denominator
        if not math.isclose(
            dynamic_probability, probability, rel_tol=2e-12, abs_tol=2e-14
        ):
            raise RuntimeError("state-path DP failed the exact packet-support gate")

    moment = optimize_block_moment(
        args.input_weight, args.distance, args.packet_positions
    )
    zero_returns = zero_return_statistics(
        args.input_weight, args.distance, args.packet_positions
    )
    binary_length = PACKET_BITS * args.packet_positions
    scalar_probability = scalar_accumulator_bad_probability(
        binary_length, args.input_weight, args.distance
    )
    theta = 4.0 * math.e * args.distance / args.packet_positions
    one_lane_probability = min(1.0, theta ** (args.input_weight / (2 * PACKET_BITS)))

    ordered_support = sorted(
        support_distribution.items(), key=lambda item: item[1], reverse=True
    )
    payload = {
        "schema": "riffle-bchblockperm-parallelacc-g4-goal01-v1",
        "evidence_label": "EXACT_BLOCK_SLICE_WITH_TWO_VALID_GAP_BOUNDS",
        "parameters": {
            "bch_block_bits": BLOCK_BITS,
            "packet_bits": PACKET_BITS,
            "packets_per_bch_block": PACKETS_PER_BLOCK,
            "global_packet_positions": args.packet_positions,
            "global_binary_length": binary_length,
            "bch_word_weight": args.input_weight,
            "failure_weight_inclusive": args.distance,
        },
        "block_permutation": {
            "law": "uniform bit permutation within the 128-bit BCH block",
            "mean_nonzero_packets": sum(
                support * probability
                for support, probability in support_distribution.items()
            ),
            "minimum_nonzero_packets": min(support_distribution),
            "minimum_support_probability": support_distribution[
                min(support_distribution)
            ],
            "mode_nonzero_packets": ordered_support[0][0],
            "mode_probability": ordered_support[0][1],
            "top_support_probabilities": ordered_support[:10],
            "exact_support_gate": "PASS",
        },
        "joint_lane_moment": moment,
        "zero_return_statistics": zero_returns,
        "comparison": {
            "old_one_lane_log2_bound": math.log2(one_lane_probability),
            "old_one_lane_effective_root": one_lane_probability
            ** (1.0 / args.input_weight),
            "full_global_bit_permutation_exact_log2_probability": (
                math.log2(scalar_probability) if scalar_probability else None
            ),
            "full_global_bit_permutation_effective_root": (
                scalar_probability ** (1.0 / args.input_weight)
                if scalar_probability
                else 0.0
            ),
        },
        "scope": (
            "The result covers one BCH block of fixed binary weight. It does not "
            "compose the double-parity outer or prove a full-code distance bound."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
