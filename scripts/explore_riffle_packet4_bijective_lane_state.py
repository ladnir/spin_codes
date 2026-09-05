#!/usr/bin/env python3
"""Profile-independent diagnostic for a coded, bijective lane state.

The four long-accumulator lanes receive a state codeword in every column.
The current experiment compares two state codes:

* k=1: the repetition code with weight enumerator 1 + z^4;
* k=2: the [4,2] code {0000,0011,1101,1110}, whose weight enumerator is
  1 + z^2 + 2 z^3;
* k=3: the even-parity [4,3] code, with enumerator 1 + 6 z^2 + z^4;
* k=4: the full four-bit lane space.

The state update has a fresh nonzero field coefficient.  A cancellation in
the next state is therefore charged at probability at most 2^(-64k).  This
script deliberately discards all output weight on a region containing such
a cancellation.  On a cancellation-free region, the affine-coset moment of
the lane state supplies the live-output bound.

For a zero incoming state and a total candidate dimension a, two independent
bounds are combined: the zero-output/zero-state event is at most the packed
epoch-occupancy activation factor, and the total output moment is at most
((1+z)/2)^a because the long accumulators are invertible.  This yields a
profile-independent diagnostic; it does not assume that packed rank profiles
maximize the complete output moment.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
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
    log_two_power_minus_one,
)
from certify_riffle_packet4_fieldcheckpoint_s256_support import (
    packed_profile_region_logs,
)
from explore_riffle_packet4_longacc_fieldchecksum import (
    marked_partial_region_log_matrices,
)


DEFAULT_OUTPUT = Path(
    "constructions/"
    "riffle_bchperm_transpose_packetshuffle_longacc_bijectivelanestate_g4/"
    "receipts/profile_independent_delta09.json"
)


def state_code_polynomial(state_dimension: int, z: float) -> float:
    """Return the unnormalized maximum affine-coset weight polynomial."""
    if state_dimension == 1:
        return 1.0 + z**4
    if state_dimension == 2:
        # C={0000,0011,1101,1110}.  Its other coset enumerators are
        # 2z+z^2+z^4 and z+2z^2+z^3, both no larger on 0<=z<=1.
        return 1.0 + z**2 + 2.0 * z**3
    if state_dimension == 3:
        # The odd coset has enumerator 4z+4z^3.  Their difference is
        # (1-z)^4, so the even coset is maximal on 0<=z<=1.
        return 1.0 + 6.0 * z**2 + z**4
    if state_dimension == 4:
        return (1.0 + z) ** 4
    raise ValueError("implemented state dimensions are 1 through 4")


def packed_activation_logs(
    *,
    groups: int,
    packet_bits: int,
    packet_slots_per_epoch: int,
    epochs_per_region: int,
    state_bits: int,
) -> np.ndarray:
    """Return log E[2^(-state_bits H)] for every packed candidate rank.

    H is the number of nonempty epochs in one region.  Candidate values and
    packet positions are averaged exactly for the packed 4q+r rank profile.
    Packing is used only for this activation event, for which the existing
    packet-support coupling proves extremality.
    """
    epoch_logs = np.full(
        (packet_slots_per_epoch + 1, 2, 2), -math.inf, dtype=np.float64
    )
    epoch_logs[0, 0, 0] = 0.0
    epoch_logs[1:, 0, 0] = -state_bits * LOG2
    exact_by_nonempty_packets = regular_region_log_matrices(
        epoch_logs,
        packet_slots_per_epoch,
        epochs_per_region,
        groups,
    )
    packed = packed_profile_region_logs(
        exact_by_nonempty_packets, groups=groups, packet_bits=packet_bits
    )
    return packed[:, 0, 0]


def region_log_matrices(
    *,
    z: float,
    activation_logs: np.ndarray,
    state_dimension: int,
    packet_slots_per_epoch: int,
    epochs_per_region: int,
) -> np.ndarray:
    """Return one conservative region matrix for every candidate dimension."""
    state_bits = state_dimension * packet_slots_per_epoch
    candidates = np.arange(activation_logs.shape[0], dtype=np.float64)
    log_candidate_moment = candidates * math.log((1.0 + z) * 0.5)

    code_size = 1 << state_dimension
    log_live_epoch = (
        packet_slots_per_epoch
        * math.log(state_code_polynomial(state_dimension, z) / code_size)
    )
    # Conditioning on a nonzero state costs at most 1/(1-2^-s).
    log_live_epoch -= math.log1p(-math.ldexp(1.0, -state_bits))
    log_live_region = epochs_per_region * log_live_epoch

    # On a region containing a state cancellation, discard every output bit.
    # A union bound covers all epoch locations.
    log_reset_region = (
        math.log(epochs_per_region)
        - math.log(math.ldexp(1.0, state_bits) - 1.0)
    )

    result = np.full(
        (activation_logs.shape[0], 2, 2), -math.inf, dtype=np.float64
    )
    result[:, 0, 0] = np.minimum(activation_logs, log_candidate_moment)
    result[:, 0, 1] = log_candidate_moment
    result[:, 1, 0] = log_reset_region
    result[:, 1, 1] = np.logaddexp(log_live_region, log_reset_region)
    return result


def packed_profile_diagnostic(
    args: argparse.Namespace,
    *,
    outer_blocks: int,
    groups: int,
    packet_slots_per_epoch: int,
    distance: int,
    regular_log_mass: float,
) -> dict[str, object]:
    """Evaluate packed profiles while retaining reset-output weight.

    This calculation is stronger than the profile-independent reset union
    bound, but packing extremality for the complete two-state moment has not
    been proved.  It is therefore evidence about the construction rather
    than a certificate.
    """
    state_bits = args.state_dimension * packet_slots_per_epoch
    reset_mass = math.ldexp(1.0, -state_bits)
    best = np.full(outer_blocks + 1, math.inf)
    best_tilt = np.full(outer_blocks + 1, math.nan)
    for log_surprisal in args.log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        code_size = 1 << args.state_dimension
        live = (
            state_code_polynomial(args.state_dimension, z) / code_size
        ) ** packet_slots_per_epoch
        live /= 1.0 - reset_mass

        epoch = np.zeros(
            (packet_slots_per_epoch + 1, 2, 2), dtype=np.float64
        )
        epoch[0, 0, 0] = 1.0
        # For an actual nonzero packet input, invertibility forces a nonzero
        # epoch output, for which z is a support-only moment bound.
        epoch[1:, 0, 0] = reset_mass * z
        epoch[1:, 0, 1] = z
        # The all-zero output has mass at most 2^-s.  For every other output,
        # a random checksum cancellation contributes another 2^-s factor.
        epoch[:, 1, 0] = reset_mass * (1.0 + live)
        epoch[:, 1, 1] = live
        with np.errstate(divide="ignore"):
            epoch_logs = np.log(epoch)
        exact_regions = regular_region_log_matrices(
            epoch_logs,
            packet_slots_per_epoch,
            args.epochs_per_region,
            groups,
        )
        packed_regions = packed_profile_region_logs(
            exact_regions, groups=groups, packet_bits=args.packet_bits
        )
        # Independently of the support recurrence, a zero-start region is an
        # injective linear image of its candidate variables before the state
        # checksum.  Its complete output moment is therefore at most B(z)^a.
        # Intersect this bound with each zero-row transition.
        log_b = math.log((1.0 + z) * 0.5)
        candidate_dimensions = np.arange(
            packed_regions.shape[0], dtype=np.float64
        )
        packed_regions[:, 0, :] = np.minimum(
            packed_regions[:, 0, :],
            candidate_dimensions[:, None] * log_b,
        )
        moments = log_matrix_power_moments_batch(
            packed_regions, args.outer_bits
        )
        candidates = moments[1:] + distance * surprisal
        improved = candidates < best[1:]
        best[1:][improved] = candidates[improved]
        best_tilt[1:][improved] = log_surprisal

    rows: list[dict[str, float | int]] = []
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
                "best_log_surprisal": float(best_tilt[active_blocks]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "margin_bits": -contribution / LOG2,
            }
        )
    dominant = sorted(rows, key=lambda row: float(row["margin_bits"]))
    return {
        "status": (
            "diagnostic_only; exact packed-profile recurrence under the "
            "stated epoch envelopes, without a complete packing theorem"
        ),
        "pointwise_minimum_margin_bits": float(dominant[0]["margin_bits"]),
        "summed_lambda_bits_lower_float": -float(
            logsumexp(np.asarray(contributions))
        )
        / LOG2,
        "dominant_rows": dominant[:40],
        "rows": rows,
    }


def full_rank_group_diagnostic(
    args: argparse.Namespace,
    *,
    outer_blocks: int,
    groups: int,
    packet_slots_per_epoch: int,
    distance: int,
    regular_log_mass: float,
) -> dict[str, object]:
    """Evaluate profiles made only of complete rank-four packet groups.

    Unlike the support-only packed diagnostic, a zero-state epoch with h
    rank-four packet variables receives its exact B(z)^(4h) information-set
    moment.  This retains candidate output weight after a state reset.
    """
    state_bits = args.state_dimension * packet_slots_per_epoch
    reset_mass = math.ldexp(1.0, -state_bits)
    best = np.full(groups + 1, math.inf)
    best_tilt = np.full(groups + 1, math.nan)
    for log_surprisal in args.log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        b = (1.0 + z) * 0.5
        code_size = 1 << args.state_dimension
        live = (
            state_code_polynomial(args.state_dimension, z) / code_size
        ) ** packet_slots_per_epoch
        live /= 1.0 - reset_mass

        epoch = np.zeros(
            (packet_slots_per_epoch + 1, 2, 2), dtype=np.float64
        )
        for full_rank_packets in range(packet_slots_per_epoch + 1):
            candidate_dimension = args.packet_bits * full_rank_packets
            zero_probability = math.ldexp(1.0, -candidate_dimension)
            zero_moment = b**candidate_dimension
            if args.state_dimension == args.packet_bits:
                epoch[full_rank_packets, 0, 0] = zero_probability
            else:
                epoch[full_rank_packets, 0, 0] = (
                    zero_probability
                    + reset_mass * (zero_moment - zero_probability)
                )
            epoch[full_rank_packets, 0, 1] = zero_moment - zero_probability
        # A live-state affine coset contains the zero output for at most one
        # state.  Every other output has checksum-zero probability 2^-s.
        if args.state_dimension == args.packet_bits:
            epoch[:, 1, 0] = 1.0 / (
                math.ldexp(1.0, state_bits) - 1.0
            )
        else:
            epoch[:, 1, 0] = reset_mass * (1.0 + live)
        epoch[:, 1, 1] = live
        with np.errstate(divide="ignore"):
            epoch_logs = np.log(epoch)
        regions = regular_region_log_matrices(
            epoch_logs,
            packet_slots_per_epoch,
            args.epochs_per_region,
            groups,
        )
        moments = log_matrix_power_moments_batch(regions, args.outer_bits)
        candidates = moments[1:] + distance * surprisal
        improved = candidates < best[1:]
        best[1:][improved] = candidates[improved]
        best_tilt[1:][improved] = log_surprisal

    rows: list[dict[str, float | int]] = []
    contributions = []
    for full_groups in range(1, groups + 1):
        active_blocks = args.packet_bits * full_groups
        outer_log = (
            log_choose(outer_blocks, active_blocks)
            + active_blocks * regular_log_mass
        )
        inner_log = min(0.0, float(best[full_groups]))
        contribution = outer_log + inner_log
        contributions.append(contribution)
        rows.append(
            {
                "active_regular_outer_blocks": active_blocks,
                "full_rank_packet_groups": full_groups,
                "best_log_surprisal": float(best_tilt[full_groups]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "margin_bits": -contribution / LOG2,
            }
        )
    dominant = sorted(rows, key=lambda row: float(row["margin_bits"]))
    return {
        "status": (
            "exact for complete rank-four packet-group profiles under the "
            "state-code coset envelope"
        ),
        "pointwise_minimum_margin_bits": float(dominant[0]["margin_bits"]),
        "summed_lambda_bits_lower_float": -float(
            logsumexp(np.asarray(contributions))
        )
        / LOG2,
        "dominant_rows": dominant[:40],
        "rows": rows,
    }


def candidate_dimension_epoch_logs(
    *,
    packet_slots_per_epoch: int,
    packet_bits: int,
    partial_rank: int,
    state_dimension: int,
    z: float,
) -> np.ndarray:
    """Return matrices indexed by full groups plus one optional partial."""
    state_bits = state_dimension * packet_slots_per_epoch
    reset_mass = math.ldexp(1.0, -state_bits)
    b = (1.0 + z) * 0.5
    code_size = 1 << state_dimension
    live = (
        state_code_polynomial(state_dimension, z) / code_size
    ) ** packet_slots_per_epoch
    live /= 1.0 - reset_mass
    matrices = np.zeros(
        (packet_slots_per_epoch + 1, 2, 2), dtype=np.float64
    )
    for full_rank_packets in range(packet_slots_per_epoch + 1):
        if partial_rank and full_rank_packets == packet_slots_per_epoch:
            continue
        candidate_dimension = (
            packet_bits * full_rank_packets + partial_rank
        )
        zero_probability = math.ldexp(1.0, -candidate_dimension)
        zero_moment = b**candidate_dimension
        if state_dimension == packet_bits:
            # The full-state update is multiplication by a nonzero field
            # scalar.  A nonzero output therefore cannot map to zero.
            matrices[full_rank_packets, 0, 0] = zero_probability
        else:
            matrices[full_rank_packets, 0, 0] = (
                zero_probability
                + reset_mass * (zero_moment - zero_probability)
            )
        matrices[full_rank_packets, 0, 1] = zero_moment - zero_probability
        if state_dimension == packet_bits:
            # A fixed candidate translate cancels at most one of the nonzero
            # states.  The following field multiplication preserves zero.
            matrices[full_rank_packets, 1, 0] = 1.0 / (
                math.ldexp(1.0, state_bits) - 1.0
            )
        else:
            matrices[full_rank_packets, 1, 0] = reset_mass * (1.0 + live)
        matrices[full_rank_packets, 1, 1] = live
    with np.errstate(divide="ignore"):
        return np.log(matrices)


def packed_rank_dimension_diagnostic(
    args: argparse.Namespace,
    *,
    outer_blocks: int,
    groups: int,
    packet_slots_per_epoch: int,
    distance: int,
    regular_log_mass: float,
) -> dict[str, object]:
    """Evaluate all packed 4q+r profiles using candidate dimensions."""
    best = np.full(outer_blocks + 1, math.inf)
    best_tilt = np.full(outer_blocks + 1, math.nan)
    for log_surprisal in args.log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        full_epoch = candidate_dimension_epoch_logs(
            packet_slots_per_epoch=packet_slots_per_epoch,
            packet_bits=args.packet_bits,
            partial_rank=0,
            state_dimension=args.state_dimension,
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
        full_indices = args.packet_bits * np.arange(1, groups + 1)
        full_candidates = full_moments[1:] + distance * surprisal
        improved = full_candidates < best[full_indices]
        best[full_indices[improved]] = full_candidates[improved]
        best_tilt[full_indices[improved]] = log_surprisal

        for rank in range(1, args.packet_bits):
            marked_epoch = candidate_dimension_epoch_logs(
                packet_slots_per_epoch=packet_slots_per_epoch,
                packet_bits=args.packet_bits,
                partial_rank=rank,
                state_dimension=args.state_dimension,
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
            indices = args.packet_bits * np.arange(groups) + rank
            candidates = marked_moments + distance * surprisal
            improved = candidates < best[indices]
            best[indices[improved]] = candidates[improved]
            best_tilt[indices[improved]] = log_surprisal

    rows: list[dict[str, float | int]] = []
    contributions = []
    for active_blocks in range(1, outer_blocks + 1):
        full_groups, partial_rank = divmod(active_blocks, args.packet_bits)
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
                "full_rank_packet_groups": full_groups,
                "partial_group_rank": partial_rank,
                "best_log_surprisal": float(best_tilt[active_blocks]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "margin_bits": -contribution / LOG2,
            }
        )
    dominant = sorted(rows, key=lambda row: float(row["margin_bits"]))
    return {
        "status": (
            "exact for packed 4q+r candidate-dimension profiles under the "
            "state-code coset envelope; packing extremality remains open"
        ),
        "pointwise_minimum_margin_bits": float(dominant[0]["margin_bits"]),
        "summed_lambda_bits_lower_float": -float(
            logsumexp(np.asarray(contributions))
        )
        / LOG2,
        "dominant_rows": dominant[:40],
        "rows": rows,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    groups = outer_blocks // args.packet_bits
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
    if packet_slots_per_epoch * args.epochs_per_region != groups:
        raise ValueError("epoch geometry does not fill one packet region")

    state_bits = args.state_dimension * packet_slots_per_epoch
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(dimension) + log_eta

    activation_logs = packed_activation_logs(
        groups=groups,
        packet_bits=args.packet_bits,
        packet_slots_per_epoch=packet_slots_per_epoch,
        epochs_per_region=args.epochs_per_region,
        state_bits=state_bits,
    )
    best = np.full(outer_blocks + 1, math.inf)
    best_tilt = np.full(outer_blocks + 1, math.nan)
    for log_surprisal in args.log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        regions = region_log_matrices(
            z=z,
            activation_logs=activation_logs,
            state_dimension=args.state_dimension,
            packet_slots_per_epoch=packet_slots_per_epoch,
            epochs_per_region=args.epochs_per_region,
        )
        moments = log_matrix_power_moments_batch(regions, args.outer_bits)
        candidates = moments[1:] + distance * surprisal
        improved = candidates < best[1:]
        best[1:][improved] = candidates[improved]
        best_tilt[1:][improved] = log_surprisal

    rows: list[dict[str, float | int]] = []
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
                "best_log_surprisal": float(best_tilt[active_blocks]),
                "activation_log2": float(activation_logs[active_blocks] / LOG2),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "margin_bits": -contribution / LOG2,
            }
        )
    dominant = sorted(rows, key=lambda row: float(row["margin_bits"]))
    total_log = float(logsumexp(np.asarray(contributions)))
    packed = packed_profile_diagnostic(
        args,
        outer_blocks=outer_blocks,
        groups=groups,
        packet_slots_per_epoch=packet_slots_per_epoch,
        distance=distance,
        regular_log_mass=regular_log_mass,
    )
    full_rank = full_rank_group_diagnostic(
        args,
        outer_blocks=outer_blocks,
        groups=groups,
        packet_slots_per_epoch=packet_slots_per_epoch,
        distance=distance,
        regular_log_mass=regular_log_mass,
    )
    packed_rank_dimension = packed_rank_dimension_diagnostic(
        args,
        outer_blocks=outer_blocks,
        groups=groups,
        packet_slots_per_epoch=packet_slots_per_epoch,
        distance=distance,
        regular_log_mass=regular_log_mass,
    )
    return {
        "schema": "riffle-packet4-bijective-lane-state-v1",
        "candidate": (
            "BCHPerm-TransposePacketShuffle-LongAcc-"
            f"BijectiveLaneState-k{args.state_dimension}-g4"
        ),
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_bits": args.outer_bits,
            "outer_blocks": outer_blocks,
            "packet_bits": args.packet_bits,
            "packet_slots_per_epoch": packet_slots_per_epoch,
            "epochs_per_region": args.epochs_per_region,
            "state_dimension_per_column": args.state_dimension,
            "state_bits": state_bits,
            "log_surprisals": args.log_surprisals,
        },
        "method": {
            "outer": "modeled regular-spectrum density envelope",
            "activation": (
                "proved packed support coupling applied only to "
                "E[2^(-state_bits H)]"
            ),
            "zero_state_output": "invertible-map information-set moment",
            "live_state_output": "maximum affine coset of the lane state code",
            "reset": (
                "union bound over state cancellations; all output on a "
                "reset-containing region is discarded"
            ),
            "arithmetic": "nearest binary64 log semiring",
            "outward_rounded": False,
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
            "special_support_excluded": "the unique all-one outer word",
        },
        "pointwise_minimum_margin_bits": float(dominant[0]["margin_bits"]),
        "summed_lambda_bits_lower_float": -total_log / LOG2,
        "dominant_rows": dominant[:40],
        "rows": rows,
        "packed_profile_diagnostic": packed,
        "full_rank_group_diagnostic": full_rank,
        "packed_rank_dimension_diagnostic": packed_rank_dimension,
        "scope": (
            "Regular outer words under the modeled spectrum.  Exceptional "
            "all-one words, a formal state-update sampler, and outward "
            "rounding remain open."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--modeled-minimum-distance", type=int, default=38)
    parser.add_argument("--packet-bits", type=int, default=4)
    parser.add_argument("--epoch-bits", type=int, default=256)
    parser.add_argument("--epochs-per-region", type=int, default=32)
    parser.add_argument(
        "--state-dimension", type=int, choices=(1, 2, 3, 4), default=2
    )
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
        "pointwise_minimum_margin_bits,"
        f"{payload['pointwise_minimum_margin_bits']:.12f}"
    )
    print(
        "summed_lambda_bits,"
        f"{payload['summed_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
