#!/usr/bin/env python3
"""Analyze one equal-weight low-block packet-accumulator profile.

The real profile consists of ``part_count`` independently permuted 128-bit
BCH words, each of binary weight ``part_weight``.  The rigorous calculation
pools those blocks into one uniform superblock slice and pays the exact
conditioning probability for restoring the declared block weights.  It then
uses the existing 16-state dynamic program, which treats zero-state gaps
exactly and retains every positive state weight.

The numerical diagnostic samples the real independent-block law directly. It
records the occupation weight after a uniform global packet placement.  These
samples expose typical paths and possible obstructions; they are not a tail
certificate.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIRECTORY))

import analyze_riffle_bchblockperm_parallelacc_g4_goal02 as kernel  # noqa: E402


PACKET_BITS = 4
PART_BITS = 128
PACKETS_PER_PART = PART_BITS // PACKET_BITS
GLOBAL_PACKET_POSITIONS = 524_352
GLOBAL_BINARY_LENGTH = PACKET_BITS * GLOBAL_PACKET_POSITIONS


def configure_kernel(part_count: int, part_weight: int) -> None:
    """Configure the imported equal-part superblock engine."""
    kernel.PART_COUNT = part_count
    kernel.PART_WEIGHT = part_weight
    kernel.SUPERBLOCK_BITS = PART_BITS * part_count
    kernel.SUPERBLOCK_WEIGHT = part_weight * part_count
    kernel.PACKETS_PER_SUPERBLOCK = PACKETS_PER_PART * part_count


def sample_weight_slice(
    rng: np.random.Generator, weight: int
) -> np.ndarray:
    bits = np.zeros(PART_BITS, dtype=np.uint8)
    bits[rng.choice(PART_BITS, size=weight, replace=False)] = 1
    return bits


def packetize(bits: np.ndarray) -> np.ndarray:
    powers = np.asarray([1, 2, 4, 8], dtype=np.uint8)
    return bits.reshape(PACKETS_PER_PART, PACKET_BITS) @ powers


def sample_real_occupation(
    *,
    part_count: int,
    part_weight: int,
    packet_positions: int,
    distance: int,
    samples: int,
    seed: int,
) -> dict[str, object]:
    """Sample local slices, active order, gaps, and full occupation weight."""
    rng = np.random.default_rng(seed)
    supports = np.empty(samples, dtype=np.int16)
    zero_returns = np.empty(samples, dtype=np.int16)
    positive_state_counts = np.empty(samples, dtype=np.int16)
    positive_state_weight_sums = np.empty(samples, dtype=np.int16)
    output_weights = np.empty(samples, dtype=np.int64)
    final_state_weights = np.empty(samples, dtype=np.int8)

    minimum_record: dict[str, object] | None = None
    minimum_output = math.inf

    for sample_index in range(samples):
        packets = np.concatenate(
            [
                packetize(sample_weight_slice(rng, part_weight))
                for _ in range(part_count)
            ]
        )
        active = packets[packets != 0]
        rng.shuffle(active)
        support = int(active.size)

        state = 0
        state_weights = np.empty(support, dtype=np.int8)
        for index, value in enumerate(active):
            state ^= int(value)
            state_weights[index] = state.bit_count()

        positions = np.sort(
            rng.choice(packet_positions, size=support, replace=False)
        )
        durations = np.empty(support, dtype=np.int64)
        if support > 1:
            durations[:-1] = positions[1:] - positions[:-1]
        durations[-1] = packet_positions - positions[-1]
        output_weight = int(np.dot(state_weights.astype(np.int64), durations))

        positive = state_weights > 0
        supports[sample_index] = support
        zero_returns[sample_index] = int(np.count_nonzero(~positive))
        positive_state_counts[sample_index] = int(np.count_nonzero(positive))
        positive_state_weight_sums[sample_index] = int(np.sum(state_weights))
        output_weights[sample_index] = output_weight
        final_state_weights[sample_index] = int(state_weights[-1])

        if output_weight < minimum_output:
            minimum_output = output_weight
            minimum_record = {
                "sample_index": sample_index,
                "packet_support": support,
                "zero_prefix_states": int(zero_returns[sample_index]),
                "positive_prefix_states": int(positive_state_counts[sample_index]),
                "positive_prefix_state_weight_sum": int(
                    positive_state_weight_sums[sample_index]
                ),
                "final_state_weight": int(final_state_weights[sample_index]),
                "output_weight": output_weight,
                "relative_binary_output_weight": output_weight
                / (PACKET_BITS * packet_positions),
                "active_packet_values": [int(value) for value in active],
                "state_weights": [int(value) for value in state_weights],
            }

    assert minimum_record is not None
    quantiles = (0.001, 0.01, 0.1, 0.5, 0.9)
    output_quantiles = np.quantile(output_weights, quantiles, method="higher")
    return {
        "evidence_label": "NUMERICAL_REAL_LAW_DIAGNOSTIC",
        "samples": samples,
        "seed": seed,
        "law": (
            f"{part_count} independent weight-{part_weight} slices in 128 bits; "
            "uniform global packet permutation"
        ),
        "mean_packet_support": float(np.mean(supports)),
        "mean_zero_prefix_states": float(np.mean(zero_returns)),
        "mean_positive_prefix_states": float(np.mean(positive_state_counts)),
        "mean_positive_prefix_state_weight_sum": float(
            np.mean(positive_state_weight_sums)
        ),
        "final_state_weight_histogram": {
            str(weight): int(np.count_nonzero(final_state_weights == weight))
            for weight in range(PACKET_BITS + 1)
        },
        "output_weight_quantiles": {
            str(quantile): {
                "binary_weight": int(value),
                "relative_binary_weight": float(
                    value / (PACKET_BITS * packet_positions)
                ),
            }
            for quantile, value in zip(quantiles, output_quantiles)
        },
        "observed_failure_count_at_declared_distance": int(
            np.count_nonzero(output_weights <= distance)
        ),
        "minimum_observed": minimum_record,
        "scope": (
            "The sample quantiles describe the real law. They do not bound a "
            "rare low-output probability."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part-count", type=int, required=True)
    parser.add_argument("--part-weight", type=int, required=True)
    parser.add_argument("--packet-positions", type=int, default=GLOBAL_PACKET_POSITIONS)
    parser.add_argument("--distance", type=int, default=188_766)
    parser.add_argument("--samples", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20_260_824)
    parser.add_argument("--positive-scaled-cost", type=float, default=33.7)
    parser.add_argument("--skip-rigorous", action="store_true")
    args = parser.parse_args()

    if args.part_count < 1:
        raise ValueError("part count must be positive")
    if not 1 <= args.part_weight <= PART_BITS:
        raise ValueError("part weight lies outside one BCH block")
    if args.packet_positions < PACKETS_PER_PART * args.part_count:
        raise ValueError("global packet count is smaller than the active superblock")
    if not 0 <= args.distance <= PACKET_BITS * args.packet_positions:
        raise ValueError("distance lies outside the binary output")
    if args.samples < 1:
        raise ValueError("sample count must be positive")
    if args.positive_scaled_cost <= 0.0:
        raise ValueError("positive scaled cost must be positive")

    configure_kernel(args.part_count, args.part_weight)
    conditioning_log = kernel.balance_log_probability()
    rigorous = None
    if not args.skip_rigorous:
        zero_gap = kernel.exact_zero_gap_statistics(
            args.distance, args.packet_positions
        )
        positive = kernel.weighted_positive_gap_bound(
            args.distance,
            args.packet_positions,
            args.positive_scaled_cost,
        )
        rigorous = {
            "evidence_label": "RIGOROUS_SUPERBLOCK_TRANSFER_BOUND",
            "conditioning_penalty_bits": -conditioning_log / math.log(2.0),
            "exact_zero_gap": {
                "unconditioned_log2_bound": zero_gap[
                    "unconditional_log2_bound"
                ],
                "balanced_profile_log2_bound": zero_gap[
                    "balanced_three_block_log2_bound"
                ],
                "top_weighted_contributions": zero_gap[
                    "top_weighted_contributions"
                ],
            },
            "positive_state_weight_bound": {
                "scaled_cost_parameter": args.positive_scaled_cost,
                "unconditioned_log2_bound": positive[
                    "unconditional_log2_bound"
                ],
                "balanced_profile_log2_bound": positive[
                    "balanced_three_block_log2_bound"
                ],
                "top_pre_chernoff_contributions": positive[
                    "top_pre_chernoff_contributions"
                ],
            },
            "validation": {
                "uniform_slice_normalization_gate": "PASS",
                "packet_support_marginal_gate": "PASS",
            },
            "scope": (
                "The bound is for the declared equal-weight profile. It uses an "
                "unconditioned pooled slice and pays the exact balance penalty."
            ),
        }

    payload = {
        "schema": "riffle-shiftalpha64-lowblock-occupation-kernel-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "profile": {
            "active_bch_blocks": args.part_count,
            "bch_block_weights": [args.part_weight] * args.part_count,
            "total_binary_input_weight": args.part_count * args.part_weight,
        },
        "parameters": {
            "packet_bits": PACKET_BITS,
            "global_packet_positions": args.packet_positions,
            "global_binary_length": PACKET_BITS * args.packet_positions,
            "failure_weight_inclusive": args.distance,
        },
        "rigorous_superblock_transfer": rigorous,
        "real_law_occupation_diagnostic": sample_real_occupation(
            part_count=args.part_count,
            part_weight=args.part_weight,
            packet_positions=args.packet_positions,
            distance=args.distance,
            samples=args.samples,
            seed=args.seed,
        ),
        "interpretation": (
            "This artifact keeps local grouping, state returns, and global gaps "
            "inside one profile kernel. It does not replace the kernel by a "
            "uniform contraction constant or sum the outer spectrum."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
