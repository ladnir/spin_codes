#!/usr/bin/env python3
"""Exact sparse two-active diagnostic for packet-g=4 FieldCheckpoint.

The calculation separates the three possible relationships between two
active outer blocks: the same fixed packet group, different groups with the
same packet lane, and different groups with different packet lanes.  It also
computes the uniform-bit-shuffle baseline from the existing transfer code.
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
from analyze_riffle_fieldcheckpoint_twoactive_w38 import (
    category_log_moments,
    intersection_log_probabilities,
    region_matrices_up_to_two,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import LOG2


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_packetshuffle_fieldcheckpoint_g4/"
    "receipts/two_active_w38_packet_law.json"
)


def fixed_support_epoch_matrix(
    state_bits: int,
    visits: int,
    support: tuple[tuple[int, int], ...],
    z: float,
) -> np.ndarray:
    """Return the exact state-class transfer for a fixed epoch support.

    A support entry is ``(lane, visit)``.  The starting live state is uniform
    over all nonzero state words, exactly as it is after a random nonzero
    field multiplier.
    """
    lane_supports = [0] * state_bits
    for lane, visit in support:
        lane_supports[lane] ^= 1 << visit

    zero_weight = 0
    endpoint_zero_weight = 0
    endpoint_delta = 0
    all_state_moment = 1.0
    for lane_support in lane_supports:
        bit = 0
        lane_zero_weight = 0
        for visit in range(visits):
            bit ^= (lane_support >> visit) & 1
            lane_zero_weight += bit
        delta = lane_support.bit_count() & 1
        lane_one_weight = visits - lane_zero_weight
        zero_weight += lane_zero_weight
        endpoint_zero_weight += lane_one_weight if delta else lane_zero_weight
        endpoint_delta |= delta
        all_state_moment *= z**lane_zero_weight + z**lane_one_weight

    live_states = math.ldexp(1.0, state_bits) - 1.0
    zero_moment = z**zero_weight
    live_total = (all_state_moment - zero_moment) / live_states
    matrix = np.zeros((2, 2), dtype=np.float64)
    matrix[0, endpoint_delta] = zero_moment
    if endpoint_delta:
        matrix[1, 0] = z**endpoint_zero_weight / live_states
    matrix[1, 1] = live_total - matrix[1, 0]
    if matrix[1, 1] < 0.0 and matrix[1, 1] > -2e-13:
        matrix[1, 1] = 0.0
    return matrix


def _powers(matrix: np.ndarray, maximum: int) -> list[np.ndarray]:
    values = [np.eye(2)]
    for _ in range(maximum):
        values.append(values[-1] @ matrix)
    return values


def packet_region_matrices(
    state_bits: int,
    epoch_bits: int,
    epochs_per_region: int,
    packet_bits: int,
    relation: str,
    z: float,
) -> list[np.ndarray]:
    """Return R0, R1, R2 for one packet-pair relationship."""
    if epoch_bits % state_bits:
        raise ValueError("state width must divide epoch length")
    if state_bits % packet_bits:
        raise ValueError("packet width must divide state width")
    visits = epoch_bits // state_bits
    lanes_per_residue = state_bits // packet_bits
    slots_per_epoch = visits * lanes_per_residue
    slots_per_region = epochs_per_region * slots_per_epoch

    zero_epoch, averaged_one = epoch_matrices(state_bits, epoch_bits, z)
    zero_powers = _powers(zero_epoch, epochs_per_region)
    one_by_visit = [
        fixed_support_epoch_matrix(state_bits, visits, ((0, visit),), z)
        for visit in range(visits)
    ]
    one_average_check = sum(one_by_visit, np.zeros((2, 2))) / visits
    if float(np.max(np.abs(one_average_check - averaged_one))) > 2e-13:
        raise AssertionError("fixed-lane one-impulse law is not the baseline law")

    inactive = zero_powers[epochs_per_region]
    active = np.zeros((2, 2), dtype=np.float64)
    for epoch in range(epochs_per_region):
        for visit in range(visits):
            active += (
                zero_powers[epoch]
                @ one_by_visit[visit]
                @ zero_powers[epochs_per_region - 1 - epoch]
            )
    active /= epochs_per_region * visits

    two_diff_lane = {
        (first, second): fixed_support_epoch_matrix(
            state_bits,
            visits,
            ((0, first), (1, second)),
            z,
        )
        for first in range(visits)
        for second in range(visits)
    }
    two_same_lane = {
        (first, second): fixed_support_epoch_matrix(
            state_bits,
            visits,
            ((0, first), (0, second)),
            z,
        )
        for first in range(visits)
        for second in range(visits)
        if first != second
    }

    overlap = np.zeros((2, 2), dtype=np.float64)
    if relation == "same_group":
        # Both block bits occupy the same packet slot and distinct packet lanes.
        for epoch in range(epochs_per_region):
            for visit in range(visits):
                overlap += (
                    zero_powers[epoch]
                    @ two_diff_lane[visit, visit]
                    @ zero_powers[epochs_per_region - 1 - epoch]
                )
        overlap /= epochs_per_region * visits
        return [inactive, active, overlap]

    if relation not in ("different_group_same_lane", "different_group_different_lane"):
        raise ValueError(f"unknown relation: {relation}")

    # The two packet slots are a uniform ordered distinct pair.  First handle
    # different epochs.  The factor two selects which labeled block is earlier.
    for first_epoch in range(epochs_per_region):
        for second_epoch in range(first_epoch + 1, epochs_per_region):
            gap = second_epoch - first_epoch - 1
            tail = epochs_per_region - 1 - second_epoch
            for first_visit in range(visits):
                for second_visit in range(visits):
                    overlap += (
                        2.0
                        * lanes_per_residue**2
                        * zero_powers[first_epoch]
                        @ one_by_visit[first_visit]
                        @ zero_powers[gap]
                        @ one_by_visit[second_visit]
                        @ zero_powers[tail]
                    )

    # Same-epoch ordered slot pairs.  Equal local slots are forbidden.
    for epoch in range(epochs_per_region):
        for first_visit in range(visits):
            for second_visit in range(visits):
                if first_visit == second_visit:
                    # Equal visits require distinct k values, hence distinct
                    # accumulator lanes for either residue relationship.
                    overlap += (
                        lanes_per_residue * (lanes_per_residue - 1)
                        * zero_powers[epoch]
                        @ two_diff_lane[first_visit, second_visit]
                        @ zero_powers[epochs_per_region - 1 - epoch]
                    )
                    continue

                if relation == "different_group_same_lane":
                    overlap += (
                        lanes_per_residue
                        * zero_powers[epoch]
                        @ two_same_lane[first_visit, second_visit]
                        @ zero_powers[epochs_per_region - 1 - epoch]
                    )
                    different_lane_count = lanes_per_residue * (lanes_per_residue - 1)
                else:
                    different_lane_count = lanes_per_residue**2
                overlap += (
                    different_lane_count
                    * zero_powers[epoch]
                    @ two_diff_lane[first_visit, second_visit]
                    @ zero_powers[epochs_per_region - 1 - epoch]
                )

    overlap /= slots_per_region * (slots_per_region - 1)
    return [inactive, active, overlap]


def conditional_w38_log_moment(
    regions: list[np.ndarray],
    outer_bits: int,
    outer_weight: int,
    intersection_logs: np.ndarray,
) -> tuple[float, int]:
    categories = category_log_moments(
        regions, outer_bits, 2 * outer_weight, outer_weight
    )
    terms = []
    dominant_intersection = -1
    dominant_term = -math.inf
    for intersection in range(outer_weight + 1):
        if not math.isfinite(float(intersection_logs[intersection])):
            continue
        singletons = 2 * outer_weight - 2 * intersection
        term = float(intersection_logs[intersection]) + float(
            categories[singletons, intersection]
        )
        terms.append(term)
        if term > dominant_term:
            dominant_term = term
            dominant_intersection = intersection
    return float(logsumexp(np.asarray(terms))), dominant_intersection


def self_test() -> dict[str, float]:
    maximum_fixed_support_error = 0.0
    for z in (0.41, 0.73, 1.0):
        baseline = region_matrices_up_to_two(3, 6, 1, z)[2]
        averaged = np.zeros((2, 2), dtype=np.float64)
        positions = [(lane, visit) for visit in range(2) for lane in range(3)]
        count = 0
        for first in range(len(positions)):
            for second in range(first + 1, len(positions)):
                averaged += fixed_support_epoch_matrix(
                    3, 2, (positions[first], positions[second]), z
                )
                count += 1
        averaged /= count
        maximum_fixed_support_error = max(
            maximum_fixed_support_error,
            float(np.max(np.abs(averaged - baseline))),
        )
    if maximum_fixed_support_error > 3e-13:
        raise AssertionError("fixed-support matrices disagree with occupancy baseline")

    maximum_packet_stochastic_error = 0.0
    for relation in (
        "same_group",
        "different_group_same_lane",
        "different_group_different_lane",
    ):
        regions = packet_region_matrices(8, 16, 3, 4, relation, 1.0)
        maximum_packet_stochastic_error = max(
            maximum_packet_stochastic_error,
            max(
                float(np.max(np.abs(matrix.sum(axis=1) - 1.0)))
                for matrix in regions
            ),
        )
    if maximum_packet_stochastic_error > 3e-13:
        raise AssertionError("packet region transfer is not stochastic at z=1")
    return {
        "maximum_fixed_support_occupancy_error": maximum_fixed_support_error,
        "maximum_packet_region_stochastic_error": maximum_packet_stochastic_error,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    if outer_blocks % args.packet_bits:
        raise ValueError("packet width must divide the outer-block count")
    groups = outer_blocks // args.packet_bits
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    spectrum_log = float(spectrum[args.outer_weight])
    intersection_logs = intersection_log_probabilities(
        args.outer_bits, args.outer_weight
    )

    pair_counts = {
        "same_group": groups * math.comb(args.packet_bits, 2),
        "different_group_same_lane": args.packet_bits * math.comb(groups, 2),
    }
    pair_counts["different_group_different_lane"] = (
        math.comb(outer_blocks, 2) - sum(pair_counts.values())
    )
    relations = list(pair_counts)
    best = {relation: math.inf for relation in relations}
    best_tilt = {relation: math.nan for relation in relations}
    best_intersection = {relation: -1 for relation in relations}
    baseline_best = math.inf
    baseline_tilt = math.nan
    baseline_intersection = -1

    grid_count = int(
        math.floor((args.grid_max - args.grid_min) / args.grid_step + 0.5)
    ) + 1
    for grid_index in range(grid_count):
        log_surprisal = args.grid_min + grid_index * args.grid_step
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        baseline_regions = region_matrices_up_to_two(
            args.state_bits, args.epoch_bits, args.epochs_per_region, z
        )
        baseline_moment, baseline_i = conditional_w38_log_moment(
            baseline_regions,
            args.outer_bits,
            args.outer_weight,
            intersection_logs,
        )
        baseline_candidate = baseline_moment + distance * surprisal
        if baseline_candidate < baseline_best:
            baseline_best = baseline_candidate
            baseline_tilt = log_surprisal
            baseline_intersection = baseline_i

        for relation in relations:
            regions = packet_region_matrices(
                args.state_bits,
                args.epoch_bits,
                args.epochs_per_region,
                args.packet_bits,
                relation,
                z,
            )
            log_moment, intersection = conditional_w38_log_moment(
                regions,
                args.outer_bits,
                args.outer_weight,
                intersection_logs,
            )
            candidate = log_moment + distance * surprisal
            if candidate < best[relation]:
                best[relation] = candidate
                best_tilt[relation] = log_surprisal
                best_intersection[relation] = intersection
        print(
            f"grid,{grid_index + 1},{grid_count},log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    outer_word_log = 2.0 * spectrum_log
    rows = []
    packet_terms = []
    for relation in relations:
        inner_log = min(0.0, best[relation])
        outer_log = math.log(pair_counts[relation]) + outer_word_log
        pointwise = outer_log + inner_log
        packet_terms.append(pointwise)
        rows.append(
            {
                "relation": relation,
                "outer_block_pair_count": pair_counts[relation],
                "outer_log2_multiplicity": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": pointwise / LOG2,
                "lambda_bits_lower_for_class": -pointwise / LOG2,
                "best_log_surprisal": best_tilt[relation],
                "dominant_outer_support_intersection": best_intersection[relation],
            }
        )

    packet_total = float(logsumexp(np.asarray(packet_terms)))
    baseline_inner = min(0.0, baseline_best)
    baseline_total = (
        math.log(math.comb(outer_blocks, 2)) + outer_word_log + baseline_inner
    )
    return {
        "schema": "riffle-packet4-fieldcheckpoint-two-active-w38-v1",
        "candidate": "BCHPerm-TransposePacketShuffle-g4-FieldCheckpoint",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "relative_distance": args.relative_distance,
            "distance": distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "outer_weight": args.outer_weight,
            "packet_bits": args.packet_bits,
            "groups": groups,
            "state_bits": args.state_bits,
            "epoch_bits": args.epoch_bits,
            "epochs_per_region": args.epochs_per_region,
        },
        "self_test": self_test(),
        "one_active_transfer_result": (
            "Exact equality with bit shuffle: restricting the single impulse to "
            "one packet residue preserves the uniform visit law, and lane labels "
            "do not enter the two-state transfer."
        ),
        "packet_relation_rows": rows,
        "packet_two_active_log2_upper": packet_total / LOG2,
        "packet_two_active_lambda_bits_lower_float": -packet_total / LOG2,
        "uniform_bitshuffle_baseline": {
            "inner_log2_upper": baseline_inner / LOG2,
            "two_active_log2_upper": baseline_total / LOG2,
            "two_active_lambda_bits_lower_float": -baseline_total / LOG2,
            "best_log_surprisal": baseline_tilt,
            "dominant_outer_support_intersection": baseline_intersection,
        },
        "scope": (
            "Exact packet-position averaging and exact two-state inner transfer "
            "for the modeled outer weight-38 shell with exactly two active outer "
            "blocks. Floating-point Chernoff arithmetic is not outward rounded."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--outer-weight", type=int, default=38)
    parser.add_argument("--packet-bits", type=int, default=4)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--epoch-bits", type=int, default=256)
    parser.add_argument("--epochs-per-region", type=int, default=32)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--grid-min", type=float, default=-10.0)
    parser.add_argument("--grid-max", type=float, default=-5.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "packet_two_active_w38_lambda_bits,"
        f"{payload['packet_two_active_lambda_bits_lower_float']:.12f}"
    )
    print(
        "bitshuffle_two_active_w38_lambda_bits,"
        f"{payload['uniform_bitshuffle_baseline']['two_active_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
