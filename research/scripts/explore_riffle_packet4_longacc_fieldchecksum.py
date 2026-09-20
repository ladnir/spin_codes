#!/usr/bin/env python3
"""Regular-word diagnostic for Packet4 LongAcc-FieldChecksum s=64.

One epoch contains 64 packets.  The four packet lanes form four binary words
of length 64.  The incoming 64-bit state is added to every lane word, and
each lane passes through a length-64 accumulator.  Four fresh coefficients
in GF(2^64) compress the four accumulated words to the next 64-bit state.

For every fixed nonzero four-word output, the random field checksum is a
uniform 64-bit word.  A refreshed incoming state therefore gives a uniform
translate of the fourfold binary repetition code.  The resulting live-state
moment is at most ((1+z^4)/2)^64 for every fixed epoch input.

This script combines that live-state envelope with a coarse zero-state fact:
any nonempty epoch emits at least one bit because the accumulator is
invertible.  The packet-count recurrence and outer reduction match the
existing packet4 s=256 calculation.  The result is a floating-point regular
word diagnostic; exceptional all-one outer words and outward rounding remain
outside its scope.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

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
from certify_riffle_packet4_fieldcheckpoint_s256_support import (
    packed_profile_region_logs,
)


DEFAULT_OUTPUT = Path(
    "constructions/"
    "riffle_bchperm_transpose_packetshuffle_longacc_fieldchecksum_g4_s64/"
    "receipts/regular_coarse_delta09.json"
)


def longacc_epoch_log_matrices(
    *,
    packet_slots_per_epoch: int,
    state_bits: int,
    packet_bits: int,
    z: float,
) -> np.ndarray:
    """Return the count-indexed two-state epoch envelope.

    State zero denotes a deterministic zero state.  State fresh denotes an
    independent uniform state in F_2^64, including zero.  The row index is
    the incoming class, and the column index is the outgoing class.
    """
    if state_bits != packet_slots_per_epoch:
        raise ValueError("the repeated-state model requires 64 packet slots")
    if packet_bits != 4:
        raise ValueError("the fourfold repetition moment requires g=4")

    reset_mass = math.ldexp(1.0, -state_bits)
    live_moment = ((1.0 + z**packet_bits) * 0.5) ** state_bits
    matrices = []
    for active_packets in range(packet_slots_per_epoch + 1):
        if active_packets == 0:
            zero_row = (1.0, 0.0)
        else:
            # A nonzero epoch input has a nonzero accumulated output.  Its
            # random field checksum is a fresh uniform 64-bit state.
            zero_row = (0.0, z)
        matrix = np.asarray(
            (
                zero_row,
                # For fixed input, at most one of the 2^64 state values can
                # cancel all four output words.  The second entry keeps the
                # complete coset-moment envelope, so this row overcounts by
                # at most 2^-64 at z=1.
                (reset_mass, live_moment),
            ),
            dtype=np.float64,
        )
        with np.errstate(divide="ignore"):
            matrices.append(np.log(matrix))
    return np.stack(matrices)


def deterministic_support_region_envelope(
    *,
    packet_slots_per_epoch: int,
    epochs_per_region: int,
    groups: int,
    state_bits: int,
    packet_bits: int,
    z: float,
) -> np.ndarray:
    """Return the entrywise suffix envelope over packet counts."""
    epoch_logs = longacc_epoch_log_matrices(
        packet_slots_per_epoch=packet_slots_per_epoch,
        state_bits=state_bits,
        packet_bits=packet_bits,
        z=z,
    )
    exact_by_count = regular_region_log_matrices(
        epoch_logs,
        packet_slots_per_epoch,
        epochs_per_region,
        groups,
    )
    return np.maximum.accumulate(exact_by_count[::-1], axis=0)[::-1]


def full_rank_epoch_log_matrices(
    *, packet_slots_per_epoch: int, state_bits: int, packet_bits: int, z: float
) -> np.ndarray:
    """Return exact epoch matrices for k rank-four packets and 64-k zeros."""
    states = 1 << packet_bits
    complement = states - 1
    weights = np.asarray(
        [z ** state.bit_count() for state in range(states)], dtype=np.float64
    )

    def moments(fresh: bool) -> np.ndarray:
        # The active-count index counts choices of packet positions.  Random
        # candidate bits are averaged inside each active transition.
        current = np.zeros((packet_slots_per_epoch + 1, states))
        current[0, 0] = 1.0
        for position in range(packet_slots_per_epoch):
            following = np.zeros_like(current)
            maximum_used = position
            if fresh:
                for state in range(states):
                    following[: maximum_used + 1, state] += (
                        0.5 * current[: maximum_used + 1, state] * weights[state]
                    )
                    following[: maximum_used + 1, state ^ complement] += (
                        0.5
                        * current[: maximum_used + 1, state]
                        * weights[state ^ complement]
                    )
            else:
                following[: maximum_used + 1] += (
                    current[: maximum_used + 1] * weights[None, :]
                )

            active_mass = (
                np.sum(current[: maximum_used + 1], axis=1) / states
            )
            following[1 : maximum_used + 2] += (
                active_mass[:, None] * weights[None, :]
            )
            current = following

        result = np.sum(current, axis=1)
        for active_packets in range(packet_slots_per_epoch + 1):
            result[active_packets] /= math.comb(
                packet_slots_per_epoch, active_packets
            )
        return result

    zero_moments = moments(False)
    fresh_moments = moments(True)
    matrices = []
    for active_packets in range(packet_slots_per_epoch + 1):
        zero_probability = math.ldexp(1.0, -packet_bits * active_packets)
        fresh_zero_probability = math.ldexp(
            1.0,
            -(
                state_bits
                + (packet_bits - 1) * active_packets
            ),
        )
        matrix = np.asarray(
            (
                (
                    zero_probability,
                    max(0.0, zero_moments[active_packets] - zero_probability),
                ),
                (
                    fresh_zero_probability,
                    max(
                        0.0,
                        fresh_moments[active_packets] - fresh_zero_probability,
                    ),
                ),
            ),
            dtype=np.float64,
        )
        with np.errstate(divide="ignore"):
            matrices.append(np.log(matrix))
    return np.stack(matrices)


def fixed_rank_epoch_log_matrices(
    *, packet_slots_per_epoch: int, packet_bits: int, rank: int, z: float
) -> np.ndarray:
    """Return exact matrices for k packets carrying a fixed rank-r lane set."""
    if not 1 <= rank <= packet_bits:
        raise ValueError("rank must be between one and the packet width")
    states = 1 << packet_bits
    complement = states - 1
    indices = np.arange(states)
    weights = np.asarray(
        [z ** state.bit_count() for state in range(states)], dtype=np.float64
    )

    def moments(fresh: bool) -> np.ndarray:
        current = np.zeros((packet_slots_per_epoch + 1, states))
        current[0, 0] = 1.0
        inactive_offsets = (0, complement) if fresh else (0,)
        active_offsets = set(range(1 << rank))
        if fresh:
            active_offsets |= {value ^ complement for value in active_offsets}
        for position in range(packet_slots_per_epoch):
            following = np.zeros_like(current)
            maximum_used = position
            for offset in inactive_offsets:
                following[: maximum_used + 1] += (
                    current[: maximum_used + 1, indices ^ offset]
                    * weights[None, :]
                    / len(inactive_offsets)
                )
            for offset in active_offsets:
                following[1 : maximum_used + 2] += (
                    current[: maximum_used + 1, indices ^ offset]
                    * weights[None, :]
                    / len(active_offsets)
                )
            current = following
        result = np.sum(current, axis=1)
        for active_packets in range(packet_slots_per_epoch + 1):
            result[active_packets] /= math.comb(
                packet_slots_per_epoch, active_packets
            )
        return result

    zero_moments = moments(False)
    fresh_moments = moments(True)
    matrices = []
    for active_packets in range(packet_slots_per_epoch + 1):
        zero_probability = math.ldexp(1.0, -rank * active_packets)
        fresh_zero_probability = math.ldexp(
            1.0,
            -(
                packet_slots_per_epoch
                + min(rank, packet_bits - 1) * active_packets
            ),
        )
        matrix = np.asarray(
            (
                (
                    zero_probability,
                    max(0.0, zero_moments[active_packets] - zero_probability),
                ),
                (
                    fresh_zero_probability,
                    max(
                        0.0,
                        fresh_moments[active_packets]
                        - fresh_zero_probability,
                    ),
                ),
            ),
            dtype=np.float64,
        )
        with np.errstate(divide="ignore"):
            matrices.append(np.log(matrix))
    return np.stack(matrices)


def layered_epoch_log_matrices(
    *,
    packet_slots_per_epoch: int,
    packet_bits: int,
    base_rank: int,
    z: float,
) -> np.ndarray:
    """Return exact matrices with rank b or b+1 at every packet position."""
    if not 0 <= base_rank < packet_bits:
        raise ValueError("base rank must be between zero and g-1")
    states = 1 << packet_bits
    complement = states - 1
    weights = np.asarray(
        [z ** state.bit_count() for state in range(states)], dtype=np.float64
    )

    def offsets(rank: int, fresh: bool) -> tuple[int, ...]:
        values = set(range(1 << rank))
        if fresh:
            values |= {value ^ complement for value in values}
        return tuple(values)

    def moments(fresh: bool) -> np.ndarray:
        low_offsets = offsets(base_rank, fresh)
        high_offsets = offsets(base_rank + 1, fresh)
        current = np.zeros((packet_slots_per_epoch + 1, states))
        current[0, 0] = 1.0
        for position in range(packet_slots_per_epoch):
            following = np.zeros_like(current)
            maximum_high = position
            for offset in low_offsets:
                destinations = np.arange(states) ^ offset
                following[: maximum_high + 1] += (
                    current[: maximum_high + 1, destinations]
                    * weights[None, :]
                    / len(low_offsets)
                )
            for offset in high_offsets:
                destinations = np.arange(states) ^ offset
                following[1 : maximum_high + 2] += (
                    current[: maximum_high + 1, destinations]
                    * weights[None, :]
                    / len(high_offsets)
                )
            current = following
        result = np.sum(current, axis=1)
        for high_count in range(packet_slots_per_epoch + 1):
            result[high_count] /= math.comb(
                packet_slots_per_epoch, high_count
            )
        return result

    zero_moments = moments(False)
    fresh_moments = moments(True)
    matrices = []
    for high_count in range(packet_slots_per_epoch + 1):
        candidate_bits = (
            base_rank * packet_slots_per_epoch + high_count
        )
        zero_probability = math.ldexp(1.0, -candidate_bits)
        if base_rank < packet_bits - 1:
            fresh_zero_exponent = packet_slots_per_epoch + candidate_bits
        else:
            fresh_zero_exponent = packet_bits * packet_slots_per_epoch
        fresh_zero_probability = math.ldexp(1.0, -fresh_zero_exponent)
        matrix = np.asarray(
            (
                (
                    zero_probability,
                    max(0.0, zero_moments[high_count] - zero_probability),
                ),
                (
                    fresh_zero_probability,
                    max(
                        0.0,
                        fresh_moments[high_count]
                        - fresh_zero_probability,
                    ),
                ),
            )
        )
        with np.errstate(divide="ignore"):
            matrices.append(np.log(matrix))
    return np.stack(matrices)


def exhaustive_full_rank_epoch_self_test() -> dict[str, object]:
    """Compare the rank-four epoch DP with exhaustive toy enumeration."""
    length = 3
    packet_bits = 4
    packet_mask = (1 << packet_bits) - 1
    tested_z = (0.37, 1.0)
    maximum_error = 0.0
    cases = 0
    for z in tested_z:
        full_rank_logs = full_rank_epoch_log_matrices(
            packet_slots_per_epoch=length,
            state_bits=length,
            packet_bits=packet_bits,
            z=z,
        )
        generic_logs = fixed_rank_epoch_log_matrices(
            packet_slots_per_epoch=length,
            packet_bits=packet_bits,
            rank=packet_bits,
            z=z,
        )
        generic_error = float(
            np.max(np.abs(np.exp(full_rank_logs) - np.exp(generic_logs)))
        )
        maximum_error = max(maximum_error, generic_error)
        if generic_error > 2e-13:
            raise AssertionError(
                f"generic rank-four epoch mismatch at z={z}: {generic_error}"
            )
        actual = np.exp(full_rank_logs)
        for active_count in range(length + 1):
            expected = np.zeros((2, 2), dtype=np.float64)
            subsets = tuple(combinations(range(length), active_count))
            for fresh in (False, True):
                q_words = range(1 << length) if fresh else (0,)
                denominator = (
                    len(subsets)
                    * (1 << (packet_bits * active_count))
                    * len(q_words)
                )
                for active_positions in subsets:
                    for packet_values in product(
                        range(1 << packet_bits), repeat=active_count
                    ):
                        inputs = dict(zip(active_positions, packet_values))
                        for q_word in q_words:
                            accumulator = 0
                            output_weight = 0
                            for position in range(length):
                                repeated_q = (
                                    packet_mask
                                    if (q_word >> position) & 1
                                    else 0
                                )
                                accumulator ^= inputs.get(position, 0) ^ repeated_q
                                output_weight += accumulator.bit_count()
                            outgoing = int(output_weight != 0)
                            expected[int(fresh), outgoing] += z**output_weight
                expected[int(fresh)] /= denominator
            error = float(np.max(np.abs(expected - actual[active_count])))
            maximum_error = max(maximum_error, error)
            cases += 1
            if error > 2e-13:
                raise AssertionError(
                    f"full-rank epoch mismatch at z={z}, k={active_count}: "
                    f"{error}"
                )
    return {
        "status": "passed",
        "toy_packet_slots": length,
        "z_values": tested_z,
        "cases": cases,
        "maximum_absolute_error": maximum_error,
    }


def aligned_rank_one_diagnostic(
    args: argparse.Namespace,
    *, outer_blocks: int,
    groups: int,
    distance: int,
    regular_log_mass: float,
) -> dict[str, object]:
    """Evaluate active blocks in distinct packets but one common lane."""
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
    best = np.full(groups + 1, math.inf)
    best_tilt = np.full(groups + 1, math.nan)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        epoch_logs = fixed_rank_epoch_log_matrices(
            packet_slots_per_epoch=packet_slots_per_epoch,
            packet_bits=args.packet_bits,
            rank=1,
            z=z,
        )
        region_logs = regular_region_log_matrices(
            epoch_logs,
            packet_slots_per_epoch,
            args.epochs_per_region,
            groups,
        )
        moments = log_matrix_power_moments_batch(region_logs, args.outer_bits)
        candidates = moments[1:] + distance * surprisal
        improved = candidates < best[1:]
        best[1:][improved] = candidates[improved]
        best_tilt[1:][improved] = log_surprisal
        print(
            f"aligned_rank1_tilt,{tilt_index + 1},"
            f"{len(args.log_surprisals)},log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    for active_blocks in range(1, groups + 1):
        outer_log = (
            log_choose(outer_blocks, active_blocks)
            + active_blocks * regular_log_mass
        )
        inner_log = min(0.0, float(best[active_blocks]))
        rows.append(
            {
                "active_regular_outer_blocks": active_blocks,
                "best_log_surprisal": float(best_tilt[active_blocks]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": (outer_log + inner_log) / LOG2,
                "margin_bits": -(outer_log + inner_log) / LOG2,
            }
        )
    dominant = sorted(rows, key=lambda row: row["pointwise_log2_upper"], reverse=True)
    return {
        "status": (
            "diagnostic_only; active blocks occupy distinct packet groups "
            "in one common lane"
        ),
        "minimum_pointwise_margin_bits": dominant[0]["margin_bits"],
        "dominant_rows": dominant[:30],
        "rows": rows,
    }


def uniform_rank_group_diagnostic(
    args: argparse.Namespace,
    *,
    rank: int,
    outer_blocks: int,
    groups: int,
    distance: int,
    regular_log_mass: float,
) -> dict[str, object]:
    """Evaluate profiles whose occupied packet groups all have one rank."""
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
    maximum_groups = min(groups, outer_blocks // rank)
    best = np.full(maximum_groups + 1, math.inf)
    best_tilt = np.full(maximum_groups + 1, math.nan)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        epoch_logs = fixed_rank_epoch_log_matrices(
            packet_slots_per_epoch=packet_slots_per_epoch,
            packet_bits=args.packet_bits,
            rank=rank,
            z=z,
        )
        region_logs = regular_region_log_matrices(
            epoch_logs,
            packet_slots_per_epoch,
            args.epochs_per_region,
            maximum_groups,
        )
        moments = log_matrix_power_moments_batch(region_logs, args.outer_bits)
        candidates = moments[1:] + distance * surprisal
        improved = candidates < best[1:]
        best[1:][improved] = candidates[improved]
        best_tilt[1:][improved] = log_surprisal
        print(
            f"uniform_rank{rank}_tilt,{tilt_index + 1},"
            f"{len(args.log_surprisals)},log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    for occupied_groups in range(1, maximum_groups + 1):
        active_blocks = rank * occupied_groups
        outer_log = (
            log_choose(outer_blocks, active_blocks)
            + active_blocks * regular_log_mass
        )
        inner_log = min(0.0, float(best[occupied_groups]))
        rows.append(
            {
                "active_regular_outer_blocks": active_blocks,
                "occupied_packet_groups": occupied_groups,
                "rank_per_group": rank,
                "best_log_surprisal": float(best_tilt[occupied_groups]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": (outer_log + inner_log) / LOG2,
                "margin_bits": -(outer_log + inner_log) / LOG2,
            }
        )
    dominant = sorted(rows, key=lambda row: row["pointwise_log2_upper"], reverse=True)
    return {
        "status": f"diagnostic_only; every occupied group has rank {rank}",
        "minimum_pointwise_margin_bits": dominant[0]["margin_bits"],
        "dominant_rows": dominant[:30],
        "rows": rows,
    }


def maximally_split_layered_diagnostic(
    args: argparse.Namespace,
    *,
    outer_blocks: int,
    groups: int,
    distance: int,
    regular_log_mass: float,
) -> dict[str, object]:
    """Evaluate the profile that fills all packet groups one lane at a time."""
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
    best = np.full(outer_blocks + 1, math.inf)
    best_tilt = np.full(outer_blocks + 1, math.nan)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        for base_rank in range(args.packet_bits):
            epoch_logs = layered_epoch_log_matrices(
                packet_slots_per_epoch=packet_slots_per_epoch,
                packet_bits=args.packet_bits,
                base_rank=base_rank,
                z=z,
            )
            region_logs = regular_region_log_matrices(
                epoch_logs,
                packet_slots_per_epoch,
                args.epochs_per_region,
                groups,
            )
            moments = log_matrix_power_moments_batch(
                region_logs, args.outer_bits
            )
            active_counts = base_rank * groups + np.arange(groups + 1)
            candidates = moments + distance * surprisal
            improved = candidates < best[active_counts]
            best[active_counts[improved]] = candidates[improved]
            best_tilt[active_counts[improved]] = log_surprisal
        print(
            f"maximally_split_tilt,{tilt_index + 1},"
            f"{len(args.log_surprisals)},log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    for active_blocks in range(1, outer_blocks + 1):
        base_rank, higher_rank_groups = divmod(active_blocks, groups)
        if base_rank == args.packet_bits:
            base_rank = args.packet_bits - 1
            higher_rank_groups = groups
        outer_log = (
            log_choose(outer_blocks, active_blocks)
            + active_blocks * regular_log_mass
        )
        inner_log = min(0.0, float(best[active_blocks]))
        rows.append(
            {
                "active_regular_outer_blocks": active_blocks,
                "base_rank_in_every_group": base_rank,
                "groups_with_one_extra_lane": higher_rank_groups,
                "best_log_surprisal": float(best_tilt[active_blocks]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": (outer_log + inner_log) / LOG2,
                "margin_bits": -(outer_log + inner_log) / LOG2,
            }
        )
    dominant = sorted(rows, key=lambda row: row["pointwise_log2_upper"], reverse=True)
    total_log = float(
        logsumexp(
            np.asarray([row["pointwise_log2_upper"] * LOG2 for row in rows])
        )
    )
    return {
        "status": (
            "exact for the maximally split nested-lane profile under the "
            "regular outer spectrum-density model"
        ),
        "pointwise_minimum_margin_bits": dominant[0]["margin_bits"],
        "summed_lambda_bits": -total_log / LOG2,
        "dominant_rows": dominant[:30],
        "rows": rows,
    }


def two_extreme_sparse_envelope_diagnostic(
    args: argparse.Namespace,
    *,
    maximum_active_blocks: int,
    outer_blocks: int,
    distance: int,
    regular_log_mass: float,
) -> dict[str, object]:
    """Take max(packed, split) before optimizing the tilt for sparse counts."""
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
    maximum_full_groups = maximum_active_blocks // args.packet_bits
    best = np.full(maximum_active_blocks + 1, math.inf)
    best_tilt = np.full(maximum_active_blocks + 1, math.nan)
    worst_geometry = [""] * (maximum_active_blocks + 1)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        full_epoch = full_rank_epoch_log_matrices(
            packet_slots_per_epoch=packet_slots_per_epoch,
            state_bits=args.state_bits,
            packet_bits=args.packet_bits,
            z=z,
        )
        full_regions = regular_region_log_matrices(
            full_epoch,
            packet_slots_per_epoch,
            args.epochs_per_region,
            maximum_full_groups,
        )
        packed = np.full(maximum_active_blocks + 1, -math.inf)
        full_moments = log_matrix_power_moments_batch(
            full_regions, args.outer_bits
        )
        for full_count in range(1, maximum_full_groups + 1):
            packed[args.packet_bits * full_count] = full_moments[full_count]
        for rank in range(1, args.packet_bits):
            marked_epoch = marked_partial_epoch_log_matrices(
                packet_slots_per_epoch=packet_slots_per_epoch,
                packet_bits=args.packet_bits,
                rank=rank,
                z=z,
            )
            marked_regions = marked_partial_region_log_matrices(
                full_epoch_logs=full_epoch,
                marked_epoch_logs=marked_epoch,
                packet_slots_per_epoch=packet_slots_per_epoch,
                epochs_per_region=args.epochs_per_region,
                maximum_full_groups=maximum_full_groups,
            )
            marked_moments = log_matrix_power_moments_batch(
                marked_regions, args.outer_bits
            )
            for full_count, moment in enumerate(marked_moments):
                active_blocks = args.packet_bits * full_count + rank
                if active_blocks <= maximum_active_blocks:
                    packed[active_blocks] = moment

        split_epoch = layered_epoch_log_matrices(
            packet_slots_per_epoch=packet_slots_per_epoch,
            packet_bits=args.packet_bits,
            base_rank=0,
            z=z,
        )
        split_regions = regular_region_log_matrices(
            split_epoch,
            packet_slots_per_epoch,
            args.epochs_per_region,
            maximum_active_blocks,
        )
        split = log_matrix_power_moments_batch(
            split_regions, args.outer_bits
        )
        for active_blocks in range(1, maximum_active_blocks + 1):
            if split[active_blocks] >= packed[active_blocks]:
                envelope = split[active_blocks]
                geometry = "maximally_split_aligned"
            else:
                envelope = packed[active_blocks]
                geometry = "maximally_packed"
            candidate = envelope + distance * surprisal
            if candidate < best[active_blocks]:
                best[active_blocks] = candidate
                best_tilt[active_blocks] = log_surprisal
                worst_geometry[active_blocks] = geometry
        print(
            f"two_extreme_sparse_tilt,{tilt_index + 1},"
            f"{len(args.log_surprisals)},log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    for active_blocks in range(1, maximum_active_blocks + 1):
        outer_log = (
            log_choose(outer_blocks, active_blocks)
            + active_blocks * regular_log_mass
        )
        inner_log = min(0.0, float(best[active_blocks]))
        rows.append(
            {
                "active_regular_outer_blocks": active_blocks,
                "best_log_surprisal": float(best_tilt[active_blocks]),
                "active_worst_geometry": worst_geometry[active_blocks],
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": (outer_log + inner_log) / LOG2,
                "margin_bits": -(outer_log + inner_log) / LOG2,
            }
        )
    total_log = float(
        logsumexp(
            np.asarray([row["pointwise_log2_upper"] * LOG2 for row in rows])
        )
    )
    return {
        "status": (
            "exact max-before-min calculation for the two proposed "
            "extremes; dominance over every mixed profile remains conjectural"
        ),
        "maximum_active_blocks": maximum_active_blocks,
        "summed_lambda_bits": -total_log / LOG2,
        "rows": rows,
    }


def single_partial_epoch_log_matrix(
    *, packet_slots_per_epoch: int, packet_bits: int, rank: int, z: float
) -> np.ndarray:
    """Return the exact epoch matrix for one uniformly placed rank-r packet."""
    if not 1 <= rank <= packet_bits:
        raise ValueError("rank must be between one and the packet width")
    states = 1 << packet_bits
    complement = states - 1
    weights = np.asarray(
        [z ** state.bit_count() for state in range(states)], dtype=np.float64
    )

    def moment(fresh: bool) -> float:
        # The first index records whether the unique active position has
        # already been used.  Active-position choices are counted, then
        # divided by the epoch length at the end.
        current = np.zeros((2, states), dtype=np.float64)
        current[0, 0] = 1.0
        inactive_offsets = (0, complement) if fresh else (0,)
        active_offsets = set(range(1 << rank))
        if fresh:
            active_offsets |= {value ^ complement for value in active_offsets}
        inactive_probability = 1.0 / len(inactive_offsets)
        active_probability = 1.0 / len(active_offsets)

        for _ in range(packet_slots_per_epoch):
            following = np.zeros_like(current)
            for used in (0, 1):
                for state in range(states):
                    mass = current[used, state]
                    if mass == 0.0:
                        continue
                    for offset in inactive_offsets:
                        destination = state ^ offset
                        following[used, destination] += (
                            mass
                            * inactive_probability
                            * weights[destination]
                        )
                    if used == 0:
                        for offset in active_offsets:
                            destination = state ^ offset
                            following[1, destination] += (
                                mass
                                * active_probability
                                * weights[destination]
                            )
            current = following
        return float(np.sum(current[1]) / packet_slots_per_epoch)

    zero_moment = moment(False)
    fresh_moment = moment(True)
    zero_probability = math.ldexp(1.0, -rank)
    fresh_zero_probability = math.ldexp(
        1.0, -(packet_slots_per_epoch + min(rank, packet_bits - 1))
    )
    matrix = np.asarray(
        (
            (zero_probability, max(0.0, zero_moment - zero_probability)),
            (
                fresh_zero_probability,
                max(0.0, fresh_moment - fresh_zero_probability),
            ),
        ),
        dtype=np.float64,
    )
    with np.errstate(divide="ignore"):
        return np.log(matrix)


def marked_partial_epoch_log_matrices(
    *, packet_slots_per_epoch: int, packet_bits: int, rank: int, z: float
) -> np.ndarray:
    """Return exact matrices for k full groups and one rank-r marked group."""
    if not 1 <= rank < packet_bits:
        raise ValueError("the marked group must have partial rank")
    states = 1 << packet_bits
    complement = states - 1
    indices = np.arange(states)
    weights = np.asarray(
        [z ** state.bit_count() for state in range(states)], dtype=np.float64
    )

    def moments(fresh: bool) -> np.ndarray:
        current = np.zeros((packet_slots_per_epoch + 1, 2, states))
        current[0, 0, 0] = 1.0
        inactive_offsets = (0, complement) if fresh else (0,)
        partial_offsets = set(range(1 << rank))
        if fresh:
            partial_offsets |= {
                value ^ complement for value in partial_offsets
            }
        for position in range(packet_slots_per_epoch):
            following = np.zeros_like(current)
            maximum_full = position
            for offset in inactive_offsets:
                following[: maximum_full + 1] += (
                    current[: maximum_full + 1, :, indices ^ offset]
                    * weights[None, None, :]
                    / len(inactive_offsets)
                )

            # A full-rank packet makes the next four-bit accumulator state
            # uniform, independently of its previous value.
            active_mass = np.sum(
                current[: maximum_full + 1], axis=2
            ) / states
            following[1 : maximum_full + 2] += (
                active_mass[:, :, None] * weights[None, None, :]
            )

            # The marked partial packet can be used at exactly one position.
            for offset in partial_offsets:
                following[: maximum_full + 1, 1] += (
                    current[: maximum_full + 1, 0, indices ^ offset]
                    * weights[None, :]
                    / len(partial_offsets)
                )
            current = following

        result = np.sum(current[:, 1], axis=1)
        for full_groups in range(packet_slots_per_epoch):
            placements = (
                math.comb(packet_slots_per_epoch, full_groups)
                * (packet_slots_per_epoch - full_groups)
            )
            result[full_groups] /= placements
        result[-1] = 0.0
        return result

    zero_moments = moments(False)
    fresh_moments = moments(True)
    matrices = []
    for full_groups in range(packet_slots_per_epoch + 1):
        if full_groups == packet_slots_per_epoch:
            matrices.append(np.full((2, 2), -math.inf))
            continue
        zero_probability = math.ldexp(
            1.0, -(packet_bits * full_groups + rank)
        )
        fresh_zero_probability = math.ldexp(
            1.0,
            -(
                packet_slots_per_epoch
                + (packet_bits - 1) * full_groups
                + rank
            ),
        )
        matrix = np.asarray(
            (
                (
                    zero_probability,
                    max(0.0, zero_moments[full_groups] - zero_probability),
                ),
                (
                    fresh_zero_probability,
                    max(
                        0.0,
                        fresh_moments[full_groups]
                        - fresh_zero_probability,
                    ),
                ),
            )
        )
        with np.errstate(divide="ignore"):
            matrices.append(np.log(matrix))
    return np.stack(matrices)


def marked_partial_region_log_matrices(
    *,
    full_epoch_logs: np.ndarray,
    marked_epoch_logs: np.ndarray,
    packet_slots_per_epoch: int,
    epochs_per_region: int,
    maximum_full_groups: int,
) -> np.ndarray:
    """Average q full groups and one marked partial group over a region."""
    identity = np.asarray(((0.0, -math.inf), (-math.inf, 0.0)))
    no_mark = np.full((maximum_full_groups + 1, 2, 2), -math.inf)
    with_mark = np.full_like(no_mark, -math.inf)
    no_mark[0] = identity

    epoch_log_choose = np.asarray(
        [log_choose(packet_slots_per_epoch, k) for k in range(packet_slots_per_epoch + 1)]
    )
    marked_log_placements = np.full(packet_slots_per_epoch + 1, -math.inf)
    for k in range(packet_slots_per_epoch):
        marked_log_placements[k] = (
            epoch_log_choose[k] + math.log(packet_slots_per_epoch - k)
        )

    current_maximum = 0
    for _ in range(epochs_per_region):
        next_maximum = min(
            maximum_full_groups, current_maximum + packet_slots_per_epoch
        )
        next_no_mark = np.full((maximum_full_groups + 1, 2, 2), -math.inf)
        next_with_mark = np.full_like(next_no_mark, -math.inf)
        for epoch_full in range(packet_slots_per_epoch + 1):
            source_count = min(
                current_maximum, next_maximum - epoch_full
            ) + 1
            if source_count <= 0:
                break
            destination = slice(epoch_full, epoch_full + source_count)
            no_then_empty = log_matmul_batch(
                no_mark[:source_count], full_epoch_logs[epoch_full]
            ) + epoch_log_choose[epoch_full]
            next_no_mark[destination] = np.logaddexp(
                next_no_mark[destination], no_then_empty
            )
            marked_then_empty = log_matmul_batch(
                with_mark[:source_count], full_epoch_logs[epoch_full]
            ) + epoch_log_choose[epoch_full]
            next_with_mark[destination] = np.logaddexp(
                next_with_mark[destination], marked_then_empty
            )
            if epoch_full < packet_slots_per_epoch:
                no_then_marked = log_matmul_batch(
                    no_mark[:source_count], marked_epoch_logs[epoch_full]
                ) + marked_log_placements[epoch_full]
                next_with_mark[destination] = np.logaddexp(
                    next_with_mark[destination], no_then_marked
                )
        no_mark = next_no_mark
        with_mark = next_with_mark
        current_maximum = next_maximum

    region_slots = packet_slots_per_epoch * epochs_per_region
    for full_groups in range(maximum_full_groups + 1):
        normalizer = (
            log_choose(region_slots, full_groups)
            + math.log(region_slots - full_groups)
        )
        with_mark[full_groups] -= normalizer
    return with_mark


@lru_cache(maxsize=None)
def exact_profile_epoch_matrix(
    *, packet_slots_per_epoch: int, packet_bits: int, masks: tuple[int, ...], z: float
) -> np.ndarray:
    """Exactly average a small fixed mask profile over distinct packet slots."""
    group_count = len(masks)
    if group_count > packet_slots_per_epoch:
        raise ValueError("more active groups than packet slots")
    states = 1 << packet_bits
    complement = states - 1
    all_used = (1 << group_count) - 1

    def moment(fresh: bool, moment_z: float) -> float:
        weights = np.asarray(
            [moment_z ** state.bit_count() for state in range(states)],
            dtype=np.float64,
        )
        current = np.zeros((1 << group_count, states), dtype=np.float64)
        current[0, 0] = 1.0
        for position in range(packet_slots_per_epoch):
            following = np.zeros_like(current)
            remaining_positions = packet_slots_per_epoch - position
            for used in range(1 << group_count):
                used_count = used.bit_count()
                zero_slots = (
                    packet_slots_per_epoch
                    - group_count
                    - position
                    + used_count
                )
                if zero_slots > 0:
                    q_offsets = (0, complement) if fresh else (0,)
                    probability = zero_slots / remaining_positions / len(q_offsets)
                    for state in range(states):
                        mass = current[used, state]
                        for offset in q_offsets:
                            destination = state ^ offset
                            following[used, destination] += (
                                mass * probability * weights[destination]
                            )
                for group, mask in enumerate(masks):
                    if (used >> group) & 1:
                        continue
                    packet_offsets = tuple(
                        value for value in range(states) if value & ~mask == 0
                    )
                    offsets = set(packet_offsets)
                    if fresh:
                        offsets |= {value ^ complement for value in offsets}
                    probability = 1.0 / remaining_positions / len(offsets)
                    destination_used = used | (1 << group)
                    for state in range(states):
                        mass = current[used, state]
                        for offset in offsets:
                            destination = state ^ offset
                            following[destination_used, destination] += (
                                mass * probability * weights[destination]
                            )
            current = following
        return float(np.sum(current[all_used]))

    matrix = np.zeros((2, 2), dtype=np.float64)
    for fresh in (False, True):
        total = moment(fresh, z)
        zero = moment(fresh, 0.0)
        matrix[int(fresh), 0] = zero
        matrix[int(fresh), 1] = max(0.0, total - zero)
    return matrix


def exact_profile_region_matrix(
    *,
    packet_slots_per_epoch: int,
    epochs_per_region: int,
    packet_bits: int,
    masks: tuple[int, ...],
    z: float,
) -> np.ndarray:
    """Exactly average a small labeled group profile over a complete region."""
    group_count = len(masks)
    region_slots = packet_slots_per_epoch * epochs_per_region
    if group_count > region_slots:
        raise ValueError("more active groups than region slots")
    subset_matrices = []
    subset_log_placements = []
    for subset in range(1 << group_count):
        subset_masks = tuple(
            masks[group]
            for group in range(group_count)
            if (subset >> group) & 1
        )
        count = len(subset_masks)
        if count > packet_slots_per_epoch:
            subset_matrices.append(None)
            subset_log_placements.append(-math.inf)
            continue
        subset_matrices.append(
            exact_profile_epoch_matrix(
                packet_slots_per_epoch=packet_slots_per_epoch,
                packet_bits=packet_bits,
                masks=subset_masks,
                z=z,
            )
        )
        subset_log_placements.append(
            math.lgamma(packet_slots_per_epoch + 1)
            - math.lgamma(packet_slots_per_epoch - count + 1)
        )

    with np.errstate(divide="ignore"):
        subset_logs = [
            None if matrix is None else np.log(matrix)
            for matrix in subset_matrices
        ]
    identity = np.asarray(((0.0, -math.inf), (-math.inf, 0.0)))
    current = np.full((1 << group_count, 2, 2), -math.inf)
    current[0] = identity
    all_used = (1 << group_count) - 1
    for _ in range(epochs_per_region):
        following = np.full_like(current, -math.inf)
        for used in range(1 << group_count):
            remaining = all_used ^ used
            subset = remaining
            while True:
                epoch_log = subset_logs[subset]
                if epoch_log is not None:
                    product_log = log_matmul_2x2(current[used], epoch_log)
                    destination = used | subset
                    following[destination] = np.logaddexp(
                        following[destination],
                        product_log + subset_log_placements[subset],
                    )
                if subset == 0:
                    break
                subset = (subset - 1) & remaining
        current = following
    log_normalizer = (
        math.lgamma(region_slots + 1)
        - math.lgamma(region_slots - group_count + 1)
    )
    return np.exp(current[all_used] - log_normalizer)


def integer_partitions(total: int, maximum: int) -> list[tuple[int, ...]]:
    if total == 0:
        return [()]
    result = []
    for first in range(min(total, maximum), 0, -1):
        for suffix in integer_partitions(total - first, first):
            result.append((first,) + suffix)
    return result


def canonical_mask_profile(
    masks: tuple[int, ...], packet_bits: int = 4
) -> tuple[int, ...]:
    """Canonicalize a profile under lane and packet-group permutations."""
    candidates = []
    for lane_permutation in permutations(range(packet_bits)):
        transformed = []
        for mask in masks:
            mapped = 0
            for source_lane, destination_lane in enumerate(lane_permutation):
                if (mask >> source_lane) & 1:
                    mapped |= 1 << destination_lane
            transformed.append(mapped)
        candidates.append(tuple(sorted(transformed)))
    return min(candidates)


def canonical_profiles_of_total_rank(
    total_rank: int, packet_bits: int = 4
) -> list[tuple[int, ...]]:
    """Enumerate mask multisets modulo lane permutations."""
    canonical_profiles: set[tuple[int, ...]] = set()

    def visit(remaining_rank: int, minimum_mask: int, prefix: tuple[int, ...]) -> None:
        if remaining_rank == 0:
            canonical_profiles.add(
                canonical_mask_profile(prefix, packet_bits)
            )
            return
        for mask in range(minimum_mask, 1 << packet_bits):
            rank = mask.bit_count()
            if rank <= remaining_rank:
                visit(remaining_rank - rank, mask, prefix + (mask,))

    visit(total_rank, 1, ())
    return sorted(canonical_profiles)


def packing_merge_toy_audit() -> dict[str, object]:
    """Test packed-profile entrywise dominance for every profile through weight 4."""
    packet_bits = 4
    packet_slots = 4
    masks_by_rank = {
        rank: tuple(
            mask
            for mask in range(1, 1 << packet_bits)
            if mask.bit_count() == rank
        )
        for rank in range(1, packet_bits + 1)
    }
    maximum_violation = 0.0
    worst_case: dict[str, object] | None = None
    profiles_checked = 0
    for z in (0.2, 0.5, 0.8, 1.0):
        for active_bits in range(1, packet_bits + 1):
            packed_masks = (
                (1 << active_bits) - 1,
            )
            packed = exact_profile_epoch_matrix(
                packet_slots_per_epoch=packet_slots,
                packet_bits=packet_bits,
                masks=packed_masks,
                z=z,
            )
            for ranks in integer_partitions(active_bits, packet_bits):
                if len(ranks) > packet_slots:
                    continue
                for masks in product(*(masks_by_rank[rank] for rank in ranks)):
                    candidate = exact_profile_epoch_matrix(
                        packet_slots_per_epoch=packet_slots,
                        packet_bits=packet_bits,
                        masks=tuple(masks),
                        z=z,
                    )
                    violation = float(np.max(candidate - packed))
                    profiles_checked += 1
                    if violation > maximum_violation:
                        maximum_violation = violation
                        worst_case = {
                            "z": z,
                            "active_bits": active_bits,
                            "ranks": ranks,
                            "masks": masks,
                            "candidate": candidate.tolist(),
                            "packed": packed.tolist(),
                        }
    return {
        "claim_tested": (
            "the single packed group entrywise dominates every fixed mask "
            "profile of the same total rank"
        ),
        "packet_slots": packet_slots,
        "maximum_active_bits": packet_bits,
        "profiles_checked": profiles_checked,
        "maximum_entrywise_violation": maximum_violation,
        "status": "passed" if maximum_violation <= 2e-13 else "refuted",
        "worst_case": worst_case,
    }


def region_packing_merge_toy_audit() -> dict[str, object]:
    """Test both matrix and repeated scalar orders after toy regions."""
    packet_bits = 4
    packet_slots = 4
    epochs = 4
    masks_by_rank = {
        rank: tuple(
            mask
            for mask in range(1, 1 << packet_bits)
            if mask.bit_count() == rank
        )
        for rank in range(1, packet_bits + 1)
    }
    maximum_violation = 0.0
    worst_case: dict[str, object] | None = None
    maximum_scalar_log2_violation = -math.inf
    worst_scalar_case: dict[str, object] | None = None
    maximum_two_extreme_log2_violation = -math.inf
    worst_two_extreme_case: dict[str, object] | None = None
    repeated_regions = 256
    profiles_checked = 0
    for z in (0.2, 0.5, 0.8, 1.0):
        for active_bits in range(1, packet_bits + 1):
            packed_masks = ((1 << active_bits) - 1,)
            packed = exact_profile_region_matrix(
                packet_slots_per_epoch=packet_slots,
                epochs_per_region=epochs,
                packet_bits=packet_bits,
                masks=packed_masks,
                z=z,
            )
            split = exact_profile_region_matrix(
                packet_slots_per_epoch=packet_slots,
                epochs_per_region=epochs,
                packet_bits=packet_bits,
                masks=(1,) * active_bits,
                z=z,
            )
            with np.errstate(divide="ignore"):
                packed_reference_scalar = float(
                    log_matrix_power_moments_batch(
                        np.log(packed)[None], repeated_regions
                    )[0]
                )
                split_reference_scalar = float(
                    log_matrix_power_moments_batch(
                        np.log(split)[None], repeated_regions
                    )[0]
                )
            extreme_reference = max(
                packed_reference_scalar, split_reference_scalar
            )
            for ranks in integer_partitions(active_bits, packet_bits):
                for masks in product(*(masks_by_rank[rank] for rank in ranks)):
                    candidate = exact_profile_region_matrix(
                        packet_slots_per_epoch=packet_slots,
                        epochs_per_region=epochs,
                        packet_bits=packet_bits,
                        masks=tuple(masks),
                        z=z,
                    )
                    violation = float(np.max(candidate - packed))
                    with np.errstate(divide="ignore"):
                        candidate_log = np.log(candidate)
                        packed_log = np.log(packed)
                    candidate_scalar = float(
                        log_matrix_power_moments_batch(
                            candidate_log[None], repeated_regions
                        )[0]
                    )
                    packed_scalar = float(
                        log_matrix_power_moments_batch(
                            packed_log[None], repeated_regions
                        )[0]
                    )
                    scalar_log2_violation = (
                        candidate_scalar - packed_scalar
                    ) / LOG2
                    two_extreme_log2_violation = (
                        candidate_scalar - extreme_reference
                    ) / LOG2
                    profiles_checked += 1
                    if violation > maximum_violation:
                        maximum_violation = violation
                        worst_case = {
                            "z": z,
                            "active_bits": active_bits,
                            "ranks": ranks,
                            "masks": masks,
                            "candidate": candidate.tolist(),
                            "packed": packed.tolist(),
                        }
                    if scalar_log2_violation > maximum_scalar_log2_violation:
                        maximum_scalar_log2_violation = scalar_log2_violation
                        worst_scalar_case = {
                            "z": z,
                            "active_bits": active_bits,
                            "ranks": ranks,
                            "masks": masks,
                            "candidate_log2_moment": candidate_scalar / LOG2,
                            "packed_log2_moment": packed_scalar / LOG2,
                        }
                    if (
                        two_extreme_log2_violation
                        > maximum_two_extreme_log2_violation
                    ):
                        maximum_two_extreme_log2_violation = (
                            two_extreme_log2_violation
                        )
                        worst_two_extreme_case = {
                            "z": z,
                            "active_bits": active_bits,
                            "ranks": ranks,
                            "masks": masks,
                            "candidate_log2_moment": candidate_scalar / LOG2,
                            "packed_log2_moment": (
                                packed_reference_scalar / LOG2
                            ),
                            "fully_split_log2_moment": (
                                split_reference_scalar / LOG2
                            ),
                        }
    return {
        "claim_tested": (
            "the packed profile entrywise dominates every fixed mask "
            "profile after a complete four-epoch toy region"
        ),
        "packet_slots_per_epoch": packet_slots,
        "epochs_per_region": epochs,
        "maximum_active_bits": packet_bits,
        "profiles_checked": profiles_checked,
        "maximum_entrywise_violation": maximum_violation,
        "entrywise_status": (
            "passed" if maximum_violation <= 2e-13 else "refuted"
        ),
        "worst_entrywise_case": worst_case,
        "repeated_regions": repeated_regions,
        "maximum_repeated_scalar_log2_violation": (
            maximum_scalar_log2_violation
        ),
        "repeated_scalar_status": (
            "passed"
            if maximum_scalar_log2_violation <= 2e-12
            else "refuted"
        ),
        "worst_repeated_scalar_case": worst_scalar_case,
        "maximum_two_extreme_log2_violation": (
            maximum_two_extreme_log2_violation
        ),
        "two_extreme_status": (
            "passed"
            if maximum_two_extreme_log2_violation <= 2e-12
            else "refuted"
        ),
        "worst_two_extreme_case": worst_two_extreme_case,
    }


def spectral_potential_toy_audit() -> dict[str, object]:
    """Search for one positive potential covering every total-rank-four profile."""
    packet_bits = 4
    packet_slots = 4
    epochs = 4
    repeated_regions = 256
    active_bits = 4
    masks_by_rank = {
        rank: tuple(
            mask
            for mask in range(1, 1 << packet_bits)
            if mask.bit_count() == rank
        )
        for rank in range(1, packet_bits + 1)
    }
    profiles = []
    for ranks in integer_partitions(active_bits, packet_bits):
        for masks in product(*(masks_by_rank[rank] for rank in ranks)):
            profiles.append((ranks, tuple(masks)))

    rows = []
    for log_surprisal in (-8.5, -8.0, -7.5, -7.0, -6.5, math.log(-math.log(0.8))):
        z = math.exp(-math.exp(log_surprisal))
        matrices = np.stack(
            [
                exact_profile_region_matrix(
                    packet_slots_per_epoch=packet_slots,
                    epochs_per_region=epochs,
                    packet_bits=packet_bits,
                    masks=masks,
                    z=z,
                )
                for _ranks, masks in profiles
            ]
        )

        def collatz_envelope(log_ratio: float) -> float:
            ratio = math.exp(log_ratio)
            first = matrices[:, 0, 0] + matrices[:, 0, 1] * ratio
            second = (
                matrices[:, 1, 0] + matrices[:, 1, 1] * ratio
            ) / ratio
            return float(max(np.max(first), np.max(second)))

        optimum = minimize_scalar(
            collatz_envelope,
            bounds=(-16.0, 16.0),
            method="bounded",
            options={"xatol": 1e-13},
        )
        potential_ratio = math.exp(float(optimum.x))
        common_lambda = float(optimum.fun)
        first_values = (
            matrices[:, 0, 0] + matrices[:, 0, 1] * potential_ratio
        )
        second_values = (
            matrices[:, 1, 0]
            + matrices[:, 1, 1] * potential_ratio
        ) / potential_ratio
        active_profile_index = int(
            np.argmax(np.maximum(first_values, second_values))
        )
        active_row = (
            "zero"
            if first_values[active_profile_index]
            >= second_values[active_profile_index]
            else "fresh"
        )

        eigenvalues = np.linalg.eigvals(matrices)
        spectral_radii = np.max(np.abs(eigenvalues), axis=1)
        maximum_spectral_index = int(np.argmax(spectral_radii))
        packed_index = profiles.index(((4,), (15,)))
        split_index = profiles.index(((1, 1, 1, 1), (1, 1, 1, 1)))
        extreme_spectral_radius = float(
            max(spectral_radii[packed_index], spectral_radii[split_index])
        )

        with np.errstate(divide="ignore"):
            exact_log_moments = log_matrix_power_moments_batch(
                np.log(matrices), repeated_regions
            )
        worst_scalar_index = int(np.argmax(exact_log_moments))
        terminal_factor = max(1.0, 1.0 / potential_ratio)
        potential_log2_bound = (
            repeated_regions * math.log2(common_lambda)
            + math.log2(terminal_factor)
        )
        rows.append(
            {
                "log_surprisal": log_surprisal,
                "z": z,
                "profiles": len(profiles),
                "potential_ratio_v1_over_v0": potential_ratio,
                "common_collatz_lambda": common_lambda,
                "common_log2_bound": potential_log2_bound,
                "active_potential_profile": {
                    "ranks": profiles[active_profile_index][0],
                    "masks": profiles[active_profile_index][1],
                    "row": active_row,
                },
                "maximum_spectral_radius": float(
                    spectral_radii[maximum_spectral_index]
                ),
                "maximum_spectral_profile": {
                    "ranks": profiles[maximum_spectral_index][0],
                    "masks": profiles[maximum_spectral_index][1],
                },
                "extreme_spectral_radius": extreme_spectral_radius,
                "spectral_extreme_gap": float(
                    spectral_radii[maximum_spectral_index]
                    - extreme_spectral_radius
                ),
                "exact_worst_log2_moment": float(
                    exact_log_moments[worst_scalar_index] / LOG2
                ),
                "exact_worst_profile": {
                    "ranks": profiles[worst_scalar_index][0],
                    "masks": profiles[worst_scalar_index][1],
                },
                "potential_overhead_bits": float(
                    potential_log2_bound
                    - exact_log_moments[worst_scalar_index] / LOG2
                ),
            }
        )
    return {
        "status": "diagnostic finite-profile potential search",
        "geometry": {
            "packet_slots_per_epoch": packet_slots,
            "epochs_per_region": epochs,
            "repeated_regions": repeated_regions,
            "total_rank": active_bits,
        },
        "rows": rows,
    }


def full_geometry_sparse_potential_audit(
    minimum_total_rank: int = 1, maximum_total_rank: int = 4
) -> dict[str, object]:
    """Exhaust a range of total ranks in the actual geometry."""
    packet_bits = 4
    packet_slots = 64
    epochs = 32
    repeated_regions = 256
    tilt_by_rank = {
        1: -8.5,
        2: -8.0,
        3: -7.5,
        4: -7.0,
        5: -7.0,
        6: -7.0,
        7: -6.5,
        8: -6.5,
        9: -6.5,
        10: -6.0,
        11: -6.0,
        12: -6.0,
    }
    rows = []
    for active_bits in range(minimum_total_rank, maximum_total_rank + 1):
        log_surprisal = tilt_by_rank[active_bits]
        z = math.exp(-math.exp(log_surprisal))
        canonical_profiles = canonical_profiles_of_total_rank(
            active_bits, packet_bits
        )
        profiles = [
            (
                tuple(
                    sorted(
                        (mask.bit_count() for mask in canonical), reverse=True
                    )
                ),
                canonical,
            )
            for canonical in canonical_profiles
        ]
        matrices = np.stack(
            [
                exact_profile_region_matrix(
                    packet_slots_per_epoch=packet_slots,
                    epochs_per_region=epochs,
                    packet_bits=packet_bits,
                    masks=masks,
                    z=z,
                )
                for _ranks, masks in profiles
            ]
        )

        def collatz_envelope(log_ratio: float) -> float:
            ratio = math.exp(log_ratio)
            return float(
                max(
                    np.max(
                        matrices[:, 0, 0]
                        + matrices[:, 0, 1] * ratio
                    ),
                    np.max(
                        (
                            matrices[:, 1, 0]
                            + matrices[:, 1, 1] * ratio
                        )
                        / ratio
                    ),
                )
            )

        def finite_horizon_objective(log_ratio: float) -> float:
            multiplier = collatz_envelope(log_ratio)
            terminal_log_factor = max(0.0, -log_ratio)
            return repeated_regions * math.log(multiplier) + terminal_log_factor

        optimum = minimize_scalar(
            finite_horizon_objective,
            bounds=(-80.0, 16.0),
            method="bounded",
            options={"xatol": 1e-13},
        )
        ratio = math.exp(float(optimum.x))
        common_lambda = collatz_envelope(float(optimum.x))
        first = matrices[:, 0, 0] + matrices[:, 0, 1] * ratio
        second = (
            matrices[:, 1, 0] + matrices[:, 1, 1] * ratio
        ) / ratio
        active_index = int(np.argmax(np.maximum(first, second)))
        with np.errstate(divide="ignore"):
            exact_log_moments = log_matrix_power_moments_batch(
                np.log(matrices), repeated_regions
            )
        worst_index = int(np.argmax(exact_log_moments))
        full_groups, partial_rank = divmod(active_bits, packet_bits)
        packed_profile = (15,) * full_groups
        if partial_rank:
            packed_profile += ((1 << partial_rank) - 1,)
        packed_masks = canonical_mask_profile(packed_profile, packet_bits)
        packed_index = next(
            index
            for index, (_ranks, masks) in enumerate(profiles)
            if masks == packed_masks
        )
        potential_bound = float(optimum.fun / LOG2)
        rows.append(
            {
                "active_bits": active_bits,
                "log_surprisal": log_surprisal,
                "z": z,
                "symmetry_classes": len(profiles),
                "potential_ratio_v1_over_v0": ratio,
                "common_collatz_lambda": common_lambda,
                "potential_objective": "256 log(lambda) + log(max(1,1/r))",
                "active_potential_profile": {
                    "ranks": profiles[active_index][0],
                    "masks": profiles[active_index][1],
                    "row": "zero" if first[active_index] >= second[active_index] else "fresh",
                },
                "packed_is_active_potential_profile": active_index == packed_index,
                "exact_worst_profile": {
                    "ranks": profiles[worst_index][0],
                    "masks": profiles[worst_index][1],
                },
                "packed_is_exact_worst_profile": worst_index == packed_index,
                "exact_worst_log2_moment": float(
                    exact_log_moments[worst_index] / LOG2
                ),
                "potential_log2_bound": potential_bound,
                "potential_overhead_bits": float(
                    potential_bound - exact_log_moments[worst_index] / LOG2
                ),
            }
        )
    return {
        "status": "exhaustive floating-point potential audit",
        "geometry": {
            "packet_slots_per_epoch": packet_slots,
            "epochs_per_region": epochs,
            "repeated_regions": repeated_regions,
            "minimum_total_rank": minimum_total_rank,
            "maximum_total_rank": maximum_total_rank,
        },
        "rows": rows,
    }


def log_matmul_2x2(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.empty((2, 2), dtype=np.float64)
    for row in range(2):
        for column in range(2):
            result[row, column] = np.logaddexp(
                left[row, 0] + right[0, column],
                left[row, 1] + right[1, column],
            )
    return result


def one_partial_packet_region_log_matrix(
    *, empty_epoch: np.ndarray, partial_epoch: np.ndarray, epochs: int
) -> np.ndarray:
    """Average one partial packet over all epoch positions in a region."""
    identity = np.asarray(((0.0, -math.inf), (-math.inf, 0.0)))
    powers = [identity]
    for _ in range(epochs):
        powers.append(log_matmul_2x2(powers[-1], empty_epoch))
    placements = []
    for epoch in range(epochs):
        placements.append(
            log_matmul_2x2(
                log_matmul_2x2(powers[epoch], partial_epoch),
                powers[epochs - epoch - 1],
            )
        )
    return logsumexp(np.stack(placements), axis=0) - math.log(epochs)


def sparse_partial_diagnostic(
    args: argparse.Namespace,
    *, outer_blocks: int,
    distance: int,
    regular_log_mass: float,
) -> dict[str, object]:
    """Evaluate packed profiles with only one rank-one to rank-three group."""
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
    best = np.full(args.packet_bits, math.inf)
    best_tilt = np.full(args.packet_bits, math.nan)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        empty_epoch = full_rank_epoch_log_matrices(
            packet_slots_per_epoch=packet_slots_per_epoch,
            state_bits=args.state_bits,
            packet_bits=args.packet_bits,
            z=z,
        )[0]
        for rank in range(1, args.packet_bits):
            partial_epoch = single_partial_epoch_log_matrix(
                packet_slots_per_epoch=packet_slots_per_epoch,
                packet_bits=args.packet_bits,
                rank=rank,
                z=z,
            )
            region = one_partial_packet_region_log_matrix(
                empty_epoch=empty_epoch,
                partial_epoch=partial_epoch,
                epochs=args.epochs_per_region,
            )
            moment = float(
                log_matrix_power_moments_batch(
                    region[None, :, :], args.outer_bits
                )[0]
            )
            candidate = moment + distance * surprisal
            if candidate < best[rank]:
                best[rank] = candidate
                best_tilt[rank] = log_surprisal
        print(
            f"partial_tilt,{tilt_index + 1},{len(args.log_surprisals)},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    for rank in range(1, args.packet_bits):
        outer_log = log_choose(outer_blocks, rank) + rank * regular_log_mass
        inner_log = min(0.0, float(best[rank]))
        rows.append(
            {
                "active_regular_outer_blocks": rank,
                "partial_group_rank": rank,
                "best_log_surprisal": float(best_tilt[rank]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": (outer_log + inner_log) / LOG2,
                "margin_bits": -(outer_log + inner_log) / LOG2,
            }
        )
    return {
        "status": "diagnostic_only; uses the packed partial-group profile",
        "minimum_pointwise_margin_bits": min(
            float(row["margin_bits"]) for row in rows
        ),
        "rows": rows,
    }


def full_group_diagnostic(
    args: argparse.Namespace,
    *,
    outer_blocks: int,
    groups: int,
    distance: int,
    regular_log_mass: float,
) -> dict[str, object]:
    """Evaluate packed profiles containing only rank-four packet groups."""
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
    best = np.full(groups + 1, math.inf)
    best_tilt = np.full(groups + 1, math.nan)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        epoch_logs = full_rank_epoch_log_matrices(
            packet_slots_per_epoch=packet_slots_per_epoch,
            state_bits=args.state_bits,
            packet_bits=args.packet_bits,
            z=z,
        )
        region_logs = regular_region_log_matrices(
            epoch_logs,
            packet_slots_per_epoch,
            args.epochs_per_region,
            groups,
        )
        moments = log_matrix_power_moments_batch(region_logs, args.outer_bits)
        candidates = moments[1:] + distance * surprisal
        improved = candidates < best[1:]
        best[1:][improved] = candidates[improved]
        best_tilt[1:][improved] = log_surprisal
        print(
            f"full_group_tilt,{tilt_index + 1},{len(args.log_surprisals)},"
            f"log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    for full_groups in range(1, groups + 1):
        active_blocks = args.packet_bits * full_groups
        outer_log = (
            log_choose(outer_blocks, active_blocks)
            + active_blocks * regular_log_mass
        )
        inner_log = min(0.0, float(best[full_groups]))
        rows.append(
            {
                "active_regular_outer_blocks": active_blocks,
                "full_packet_groups": full_groups,
                "best_log_surprisal": float(best_tilt[full_groups]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": (outer_log + inner_log) / LOG2,
            }
        )
    dominant = sorted(
        rows,
        key=lambda row: float(row["pointwise_log2_upper"]),
        reverse=True,
    )
    return {
        "status": "diagnostic_only; no packing theorem is claimed",
        "maximum_pointwise_log2_upper": dominant[0]["pointwise_log2_upper"],
        "minimum_pointwise_margin_bits": -dominant[0][
            "pointwise_log2_upper"
        ],
        "dominant_rows": dominant[:30],
        "rows": rows,
    }


def packed_mixed_profile_diagnostic(
    args: argparse.Namespace,
    *,
    outer_blocks: int,
    groups: int,
    distance: int,
    regular_log_mass: float,
) -> dict[str, object]:
    """Evaluate every packed profile 4q+r with its exact rank information."""
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
    best = np.full(outer_blocks + 1, math.inf)
    best_tilt = np.full(outer_blocks + 1, math.nan)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        full_epoch = full_rank_epoch_log_matrices(
            packet_slots_per_epoch=packet_slots_per_epoch,
            state_bits=args.state_bits,
            packet_bits=args.packet_bits,
            z=z,
        )
        full_regions = regular_region_log_matrices(
            full_epoch,
            packet_slots_per_epoch,
            args.epochs_per_region,
            groups,
        )
        full_moments = log_matrix_power_moments_batch(
            full_regions, args.outer_bits
        )
        full_candidates = full_moments[1:] + distance * surprisal
        full_indices = args.packet_bits * np.arange(1, groups + 1)
        improved = full_candidates < best[full_indices]
        best[full_indices[improved]] = full_candidates[improved]
        best_tilt[full_indices[improved]] = log_surprisal

        for rank in range(1, args.packet_bits):
            marked_epoch = marked_partial_epoch_log_matrices(
                packet_slots_per_epoch=packet_slots_per_epoch,
                packet_bits=args.packet_bits,
                rank=rank,
                z=z,
            )
            marked_regions = marked_partial_region_log_matrices(
                full_epoch_logs=full_epoch,
                marked_epoch_logs=marked_epoch,
                packet_slots_per_epoch=packet_slots_per_epoch,
                epochs_per_region=args.epochs_per_region,
                maximum_full_groups=groups - 1,
            )
            marked_moments = log_matrix_power_moments_batch(
                marked_regions, args.outer_bits
            )
            candidates = marked_moments + distance * surprisal
            indices_for_rank = (
                args.packet_bits * np.arange(groups) + rank
            )
            improved = candidates < best[indices_for_rank]
            best[indices_for_rank[improved]] = candidates[improved]
            best_tilt[indices_for_rank[improved]] = log_surprisal
        print(
            f"packed_mixed_tilt,{tilt_index + 1},"
            f"{len(args.log_surprisals)},log_surprisal,{log_surprisal:.6f}",
            flush=True,
        )

    rows = []
    for active_blocks in range(1, outer_blocks + 1):
        full_groups, partial_rank = divmod(active_blocks, args.packet_bits)
        outer_log = (
            log_choose(outer_blocks, active_blocks)
            + active_blocks * regular_log_mass
        )
        inner_log = min(0.0, float(best[active_blocks]))
        rows.append(
            {
                "active_regular_outer_blocks": active_blocks,
                "full_packet_groups": full_groups,
                "partial_group_rank": partial_rank,
                "best_log_surprisal": float(best_tilt[active_blocks]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": (outer_log + inner_log) / LOG2,
                "margin_bits": -(outer_log + inner_log) / LOG2,
            }
        )
    dominant = sorted(rows, key=lambda row: row["pointwise_log2_upper"], reverse=True)
    total_log = float(
        logsumexp(
            np.asarray([row["pointwise_log2_upper"] * LOG2 for row in rows])
        )
    )
    return {
        "status": (
            "exact for every packed rank profile under the regular outer "
            "spectrum-density model; packing extremality is not yet proved"
        ),
        "pointwise_minimum_margin_bits": dominant[0]["margin_bits"],
        "summed_lambda_bits": -total_log / LOG2,
        "dominant_rows": dominant[:30],
        "rows": rows,
    }


def accumulator_weight_counts(length: int, toggles: int) -> dict[int, int]:
    """Count output weights for fixed-weight inputs to one accumulator."""
    states: dict[tuple[int, int, int], int] = {(0, 0, 0): 1}
    for _ in range(length):
        following: dict[tuple[int, int, int], int] = defaultdict(int)
        for (used, state, weight), count in states.items():
            following[(used, state, weight + state)] += count
            if used < toggles:
                next_state = state ^ 1
                following[(used + 1, next_state, weight + next_state)] += count
        states = following
    result: dict[int, int] = defaultdict(int)
    for (used, _state, weight), count in states.items():
        if used == toggles:
            result[weight] += count
    return dict(result)


def convolve_counts(
    left: dict[int, int], right: dict[int, int]
) -> dict[int, int]:
    result: dict[int, int] = defaultdict(int)
    for left_weight, left_count in left.items():
        for right_weight, right_count in right.items():
            result[left_weight + right_weight] += left_count * right_count
    return dict(result)


def negative_log2_probability(numerator: int, denominator: int) -> float:
    if numerator == 0:
        return math.inf
    return math.log2(denominator) - math.log2(numerator)


def local_diagnostics(
    *, epoch_bits: int, packet_bits: int, relative_distance: float
) -> dict[str, object]:
    packet_slots = epoch_bits // packet_bits
    local_distance = math.floor(relative_distance * epoch_bits)
    live_low_count = sum(math.comb(packet_slots, weight) for weight in range(6))
    live_denominator = 1 << packet_slots
    live_distribution = {
        packet_bits * weight: math.comb(packet_slots, weight)
        for weight in range(packet_slots + 1)
    }

    rows = []
    for toggles in (1, 2, 4, 8, 16, 24, 32, 48, 64):
        one_lane = accumulator_weight_counts(packet_slots, toggles)
        one_lane_denominator = math.comb(packet_slots, toggles)
        one_lane_low = sum(
            count for weight, count in one_lane.items() if weight <= local_distance
        )

        lane_toggles = [
            toggles // packet_bits + int(lane < toggles % packet_bits)
            for lane in range(packet_bits)
        ]
        balanced = {0: 1}
        balanced_denominator = 1
        for lane_count in lane_toggles:
            balanced = convolve_counts(
                balanced, accumulator_weight_counts(packet_slots, lane_count)
            )
            balanced_denominator *= math.comb(packet_slots, lane_count)
        balanced_low = sum(
            count for weight, count in balanced.items() if weight <= local_distance
        )

        two_epoch_threshold = math.floor(relative_distance * 2 * epoch_bits)
        two_epoch_low = sum(
            first_count * live_count
            for first_weight, first_count in one_lane.items()
            for live_weight, live_count in live_distribution.items()
            if first_weight + live_weight <= two_epoch_threshold
        )
        three_epoch_threshold = math.floor(relative_distance * 3 * epoch_bits)
        two_live = convolve_counts(live_distribution, live_distribution)
        three_epoch_low = sum(
            first_count * live_count
            for first_weight, first_count in one_lane.items()
            for live_weight, live_count in two_live.items()
            if first_weight + live_weight <= three_epoch_threshold
        )

        rows.append(
            {
                "toggles": toggles,
                "balanced_lane_toggles": lane_toggles,
                "one_lane_low_weight_bits": negative_log2_probability(
                    one_lane_low, one_lane_denominator
                ),
                "balanced_low_weight_bits": negative_log2_probability(
                    balanced_low, balanced_denominator
                ),
                "one_lane_then_one_fresh_epoch_bits": negative_log2_probability(
                    two_epoch_low,
                    one_lane_denominator * live_denominator,
                ),
                "one_lane_then_two_fresh_epochs_bits": negative_log2_probability(
                    three_epoch_low,
                    one_lane_denominator * live_denominator**2,
                ),
            }
        )

    return {
        "epoch_distance": local_distance,
        "fresh_state_low_weight_bits": negative_log2_probability(
            live_low_count, live_denominator
        ),
        "fresh_state_exact_zero_bits": float(packet_slots),
        "canonical_zero_start_rows": rows,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    groups = outer_blocks // args.packet_bits
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
    if packet_slots_per_epoch != args.state_bits:
        raise ValueError("the current construction repeats one state bit per packet")
    if packet_slots_per_epoch * args.epochs_per_region != groups:
        raise ValueError("epoch geometry does not fill one packet region")

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
        support_envelope = deterministic_support_region_envelope(
            packet_slots_per_epoch=packet_slots_per_epoch,
            epochs_per_region=args.epochs_per_region,
            groups=groups,
            state_bits=args.state_bits,
            packet_bits=args.packet_bits,
            z=z,
        )
        packed_regions = packed_profile_region_logs(
            support_envelope,
            groups=groups,
            packet_bits=args.packet_bits,
        )
        moments = log_matrix_power_moments_batch(packed_regions, args.outer_bits)
        candidates = moments[1:] + distance * surprisal
        improved = candidates < best[1:]
        best[1:][improved] = candidates[improved]
        best_tilt[1:][improved] = log_surprisal
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
    full_groups = full_group_diagnostic(
        args,
        outer_blocks=outer_blocks,
        groups=groups,
        distance=distance,
        regular_log_mass=regular_log_mass,
    )
    packed_mixed = packed_mixed_profile_diagnostic(
        args,
        outer_blocks=outer_blocks,
        groups=groups,
        distance=distance,
        regular_log_mass=regular_log_mass,
    )
    sparse_partials = sparse_partial_diagnostic(
        args,
        outer_blocks=outer_blocks,
        distance=distance,
        regular_log_mass=regular_log_mass,
    )
    aligned_rank_one = aligned_rank_one_diagnostic(
        args,
        outer_blocks=outer_blocks,
        groups=groups,
        distance=distance,
        regular_log_mass=regular_log_mass,
    )
    uniform_rank_profiles = {
        str(rank): uniform_rank_group_diagnostic(
            args,
            rank=rank,
            outer_blocks=outer_blocks,
            groups=groups,
            distance=distance,
            regular_log_mass=regular_log_mass,
        )
        for rank in (2, 3)
    }
    maximally_split = maximally_split_layered_diagnostic(
        args,
        outer_blocks=outer_blocks,
        groups=groups,
        distance=distance,
        regular_log_mass=regular_log_mass,
    )
    return {
        "schema": "riffle-packet4-longacc-fieldchecksum-s64-geometry-v2",
        "candidate": (
            "BCHPerm-TransposePacketShuffle-LongAcc-FieldChecksum-g4-s64"
        ),
        "method": {
            "outer_reduction": "modeled regular-spectrum density envelope",
            "profile": "proved packed nonzero-packet-count upper bound",
            "zero_state_epoch": (
                "one output bit for every nonempty epoch; all other output "
                "weight is discarded"
            ),
            "fresh_state_epoch": "fourfold-repetition coset moment envelope",
            "checksum": (
                "four independent uniform GF(2^64) coefficients per epoch"
            ),
            "region_recurrence": (
                "exact conditional packet-count matrices followed by the "
                "entrywise suffix envelope max_{K>=H} R_K"
            ),
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
        "local_diagnostics": local_diagnostics(
            epoch_bits=args.epoch_bits,
            packet_bits=args.packet_bits,
            relative_distance=args.relative_distance,
        ),
        "exact_epoch_self_test": exhaustive_full_rank_epoch_self_test(),
        "packing_merge_toy_audit": packing_merge_toy_audit(),
        "region_packing_merge_toy_audit": region_packing_merge_toy_audit(),
        "packed_sparse_partial_diagnostic": sparse_partials,
        "packed_full_rank_group_diagnostic": full_groups,
        "packed_mixed_profile_diagnostic": packed_mixed,
        "aligned_rank_one_diagnostic": aligned_rank_one,
        "uniform_rank_group_diagnostics": uniform_rank_profiles,
        "maximally_split_layered_diagnostic": maximally_split,
        "regular_log2_upper": total_log / LOG2,
        "regular_lambda_bits_lower_float": -total_log / LOG2,
        "dominant_rows": sorted(
            rows,
            key=lambda row: float(row["pointwise_log2_upper"]),
            reverse=True,
        )[:30],
        "occupation_rows": rows,
        "scope": (
            "Every regular active-block count, packed packet-count profile, "
            "and zero/fresh trajectory under the stated epoch envelope.  "
            "Exceptional all-one outer words and outward rounding remain open."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--packet-bits", type=int, default=4)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--epoch-bits", type=int, default=256)
    parser.add_argument("--epochs-per-region", type=int, default=32)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(
            -10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
            -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0,
            -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5,
            0.75, 1.0,
        ),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "regular_lambda_bits,"
        f"{payload['regular_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
