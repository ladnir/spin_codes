#!/usr/bin/env python3
"""Full regular diagnostic for maximally packed four-block groups at s=256.

The 256-bit checkpoint state is visited once per epoch.  Therefore, an epoch
with A fair candidate bits has an exact two-state transfer that depends only
on A, not on candidate positions.  This diagnostic packs a active blocks into
floor(a/4) full groups and at most one partial group.  Packet permutations are
averaged exactly across all 32 epochs in each region.

The outer multiplicity deliberately counts every a-block subset.  A proof that
maximal group packing dominates every other group profile would therefore turn
this complete floating-point sum into a regular-word certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp

from analyze_riffle_fieldcheckpoint_regular_bulk_logdp import (
    log_matrix_power_moments_batch,
    regular_region_log_matrices,
)
from analyze_riffle_fieldcheckpoint_regular_envelope import (
    spectrum_density_envelope_log,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import (
    LOG2,
    log_choose,
    log_matmul_batch,
    log_two_power_minus_one,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_packetshuffle_fieldcheckpoint_g4_s256/"
    "receipts/packed_groups_full_regular_delta09.json"
)


def fair_candidate_epoch_log_matrix(state_bits: int, candidates: int, z: float) -> np.ndarray:
    """Return the exact one-visit transfer for A fair candidate coordinates."""
    live_states = math.ldexp(1.0, state_bits) - 1.0
    zero_probability = math.ldexp(1.0, -candidates)
    zero_start_moment = ((1.0 + z) * 0.5) ** candidates
    all_state_moment = (1.0 + z) ** state_bits
    matrix = np.asarray(
        (
            (zero_probability, zero_start_moment - zero_probability),
            (
                (1.0 - zero_probability) / live_states,
                (
                    all_state_moment
                    - zero_start_moment
                    - 1.0
                    + zero_probability
                )
                / live_states,
            ),
        ),
        dtype=np.float64,
    )
    matrix[(matrix < 0.0) & (matrix > -3e-13)] = 0.0
    if np.any(matrix < 0.0):
        raise ArithmeticError("fair-candidate epoch transfer is negative")
    with np.errstate(divide="ignore"):
        return np.log(matrix)


def log_choose_vector(total: int) -> np.ndarray:
    values = np.arange(total + 1, dtype=np.float64)
    return (
        gammaln(total + 1.0)
        - gammaln(values + 1.0)
        - gammaln(total - values + 1.0)
    )


def marked_partial_region_logs(
    *,
    base_epoch_logs: np.ndarray,
    marked_epoch_logs: np.ndarray,
    packet_slots_per_epoch: int,
    epochs_per_region: int,
    maximum_full_groups: int,
) -> np.ndarray:
    """Average h full groups and one labeled partial group in a region."""
    negative = -math.inf
    no_mark = np.full((maximum_full_groups + 1, 2, 2), negative)
    have_mark = np.full_like(no_mark, negative)
    no_mark[0] = np.asarray(((0.0, negative), (negative, 0.0)))
    choose_full = log_choose_vector(packet_slots_per_epoch)
    choose_with_mark = (
        math.log(packet_slots_per_epoch)
        + log_choose_vector(packet_slots_per_epoch - 1)
    )

    for completed_epochs in range(epochs_per_region):
        maximum_before = min(
            completed_epochs * packet_slots_per_epoch,
            maximum_full_groups,
        )
        maximum_after = min(
            (completed_epochs + 1) * packet_slots_per_epoch,
            maximum_full_groups,
        )
        next_no_mark = np.full_like(no_mark, negative)
        next_have_mark = np.full_like(have_mark, negative)
        for here in range(packet_slots_per_epoch + 1):
            source_count = min(maximum_before, maximum_after - here) + 1
            if source_count <= 0:
                break
            destination = slice(here, here + source_count)
            no_product = log_matmul_batch(
                no_mark[:source_count], base_epoch_logs[here]
            ) + choose_full[here]
            next_no_mark[destination] = np.logaddexp(
                next_no_mark[destination], no_product
            )
            have_product = log_matmul_batch(
                have_mark[:source_count], base_epoch_logs[here]
            ) + choose_full[here]
            next_have_mark[destination] = np.logaddexp(
                next_have_mark[destination], have_product
            )
            if here < packet_slots_per_epoch:
                marked_product = log_matmul_batch(
                    no_mark[:source_count], marked_epoch_logs[here]
                ) + choose_with_mark[here]
                next_have_mark[destination] = np.logaddexp(
                    next_have_mark[destination], marked_product
                )
        no_mark = next_no_mark
        have_mark = next_have_mark

    region_slots = packet_slots_per_epoch * epochs_per_region
    for full_groups in range(maximum_full_groups + 1):
        normalization = math.log(region_slots) + log_choose(
            region_slots - 1, full_groups
        )
        have_mark[full_groups] -= normalization
    return have_mark


def self_test(args: argparse.Namespace) -> dict[str, float]:
    stochastic_error = 0.0
    for candidates in range(args.epoch_bits + 1):
        matrix = np.exp(
            fair_candidate_epoch_log_matrix(args.state_bits, candidates, 1.0)
        )
        stochastic_error = max(
            stochastic_error,
            float(np.max(np.abs(matrix.sum(axis=1) - 1.0))),
        )
    if stochastic_error > 5e-13:
        raise AssertionError("one-visit epoch transfer is not stochastic")

    toy_state_bits = 6
    toy_packet_bits = 2
    toy_slots_per_epoch = 3
    toy_epochs = 2
    toy_z = 0.73
    toy_all = np.stack(
        [
            fair_candidate_epoch_log_matrix(toy_state_bits, count, toy_z)
            for count in range(toy_state_bits + 1)
        ]
    )
    toy_marked = marked_partial_region_logs(
        base_epoch_logs=toy_all[::toy_packet_bits],
        marked_epoch_logs=np.stack(
            [toy_all[toy_packet_bits * full + 1] for full in range(toy_slots_per_epoch)]
        ),
        packet_slots_per_epoch=toy_slots_per_epoch,
        epochs_per_region=toy_epochs,
        maximum_full_groups=5,
    )
    marked_error = 0.0
    region_slots = toy_slots_per_epoch * toy_epochs
    for full_groups in range(6):
        brute = np.zeros((2, 2), dtype=np.float64)
        cases = 0
        for partial_position in range(region_slots):
            remaining = [
                position
                for position in range(region_slots)
                if position != partial_position
            ]
            for full_positions in itertools.combinations(remaining, full_groups):
                full_set = set(full_positions)
                product = np.eye(2)
                for epoch in range(toy_epochs):
                    start = epoch * toy_slots_per_epoch
                    stop = start + toy_slots_per_epoch
                    full_here = sum(position in full_set for position in range(start, stop))
                    partial_here = start <= partial_position < stop
                    candidates = toy_packet_bits * full_here + int(partial_here)
                    product = product @ np.exp(toy_all[candidates])
                brute += product
                cases += 1
        brute /= cases
        marked_error = max(
            marked_error,
            float(np.max(np.abs(np.exp(toy_marked[full_groups]) - brute))),
        )
    if marked_error > 2e-12:
        raise AssertionError("marked partial-group recurrence failed")
    return {
        "maximum_epoch_stochastic_error": stochastic_error,
        "maximum_marked_region_bruteforce_error": marked_error,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    if args.state_bits != args.epoch_bits:
        raise ValueError("this checker requires exactly one visit per state lane")
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    groups = outer_blocks // args.packet_bits
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
    if packet_slots_per_epoch * args.epochs_per_region != groups:
        raise ValueError("epoch geometry does not fill a packet region")
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(dimension) + log_eta

    best = np.full(outer_blocks + 1, math.inf)
    best_tilt = np.full(outer_blocks + 1, math.nan)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        all_epoch_logs = np.stack(
            [
                fair_candidate_epoch_log_matrix(args.state_bits, count, z)
                for count in range(args.epoch_bits + 1)
            ]
        )
        base_epoch_logs = all_epoch_logs[:: args.packet_bits]
        full_regions = regular_region_log_matrices(
            base_epoch_logs,
            packet_slots_per_epoch,
            args.epochs_per_region,
            groups,
        )
        full_moments = log_matrix_power_moments_batch(
            full_regions, args.outer_bits
        )
        for full_groups in range(1, groups + 1):
            active_blocks = args.packet_bits * full_groups
            candidate = full_moments[full_groups] + distance * surprisal
            if candidate < best[active_blocks]:
                best[active_blocks] = candidate
                best_tilt[active_blocks] = log_surprisal

        for partial_rank in range(1, args.packet_bits):
            marked_epoch_logs = np.stack(
                [
                    all_epoch_logs[
                        args.packet_bits * full_count + partial_rank
                    ]
                    for full_count in range(packet_slots_per_epoch)
                ]
            )
            marked_regions = marked_partial_region_logs(
                base_epoch_logs=base_epoch_logs,
                marked_epoch_logs=marked_epoch_logs,
                packet_slots_per_epoch=packet_slots_per_epoch,
                epochs_per_region=args.epochs_per_region,
                maximum_full_groups=groups - 1,
            )
            marked_moments = log_matrix_power_moments_batch(
                marked_regions, args.outer_bits
            )
            for full_groups in range(groups):
                active_blocks = args.packet_bits * full_groups + partial_rank
                candidate = marked_moments[full_groups] + distance * surprisal
                if candidate < best[active_blocks]:
                    best[active_blocks] = candidate
                    best_tilt[active_blocks] = log_surprisal
        print(
            f"tilt,{tilt_index + 1},{len(args.log_surprisals)},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    contributions = []
    for active_blocks in range(1, outer_blocks + 1):
        outer_log = (
            log_choose(outer_blocks, active_blocks)
            + active_blocks * regular_log_mass
        )
        inner_log = min(0.0, float(best[active_blocks]))
        contribution = outer_log + inner_log
        contributions.append(contribution)
        rows.append(
            {
                "active_regular_outer_blocks": active_blocks,
                "full_packet_groups": active_blocks // args.packet_bits,
                "partial_group_rank": active_blocks % args.packet_bits,
                "best_log_surprisal": float(best_tilt[active_blocks]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
            }
        )
    total_log = float(logsumexp(np.asarray(contributions)))
    return {
        "schema": "riffle-packet4-fieldcheckpoint-s256-packed-groups-regular-v1",
        "candidate": "BCHPerm-TransposePacketShuffle-FieldCheckpoint-s256-g4",
        "method": {
            "outer_reduction": "modeled regular-spectrum density envelope",
            "profile": "maximally packed four-block groups",
            "epoch_transfer": "exact one-visit fair-candidate transfer",
            "region_recurrence": (
                "exact packet permutation with full groups and at most one "
                "labeled partial group"
            ),
            "outer_placement_count": "all subsets at each active-block count",
            "arithmetic": "nearest binary64 log semiring",
            "outward_rounded": False,
        },
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "groups": groups,
            "packet_bits": args.packet_bits,
            "state_bits": args.state_bits,
            "epoch_bits": args.epoch_bits,
            "epochs_per_region": args.epochs_per_region,
            "log_surprisals": args.log_surprisals,
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
            "special_support_excluded": "the unique all-one outer word",
        },
        "self_test": self_test(args),
        "packed_group_log2_upper": total_log / LOG2,
        "packed_group_lambda_bits_lower_float": -total_log / LOG2,
        "dominant_rows": sorted(
            rows,
            key=lambda row: float(row["pointwise_log2_upper"]),
            reverse=True,
        )[:30],
        "occupation_rows": rows,
        "scope": (
            "Complete floating-point sum for maximally packed regular group "
            "profiles. A group-packing domination lemma, all-one cases, and "
            "outward rounding remain open."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--packet-bits", type=int, default=4)
    parser.add_argument("--state-bits", type=int, default=256)
    parser.add_argument("--epoch-bits", type=int, default=256)
    parser.add_argument("--epochs-per-region", type=int, default=32)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0,
                 -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0,
                 0.5, 0.75, 1.0),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "packed_group_lambda_bits,"
        f"{payload['packed_group_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
