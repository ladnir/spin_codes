#!/usr/bin/env python3
"""Floating-point regular-word certificate for packet g=4, state s=256.

The epoch transfer replaces each nonzero input packet by one output bit and
discards its additional Hamming weight.  If the state is zero, the relaxed
output weight is therefore the number of nonzero packets in the epoch.  If
the state is live, the emitted/final vector is a translate of a uniform
nonzero field element.  Its zero term is at most 1/(2^256-1), and its nonzero
z-moment is at most ((1+z)^256-1)/(2^256-1), independently of the input.

The proved packet-packing coupling makes the packed profile's number of
nonzero packets stochastically smallest.  Exact region matrices need not be
monotone in that count, so the checker computes every count-conditioned
matrix and takes an entrywise suffix maximum.  This suffix envelope turns the
scalar packing comparison into a valid matrix bound.  The unique all-one
outer word is treated separately by a companion calculation.
"""

from __future__ import annotations

import argparse
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
    log_two_power_minus_one,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_packetshuffle_fieldcheckpoint_g4_s256/"
    "receipts/all_profiles_regular_support_relaxation_delta09.json"
)


def packet_indicator_epoch_log_matrix(
    state_bits: int,
    packet_bits: int,
    full_groups: int,
    partial_rank: int,
    z: float,
) -> np.ndarray:
    """Upper transfer after replacing each nonzero packet by one bit."""
    full_states = math.ldexp(1.0, state_bits)
    candidates = packet_bits * full_groups + partial_rank
    zero_input = math.ldexp(1.0, -candidates)
    full_packet_moment = z + (1.0 - z) * math.ldexp(1.0, -packet_bits)
    partial_packet_moment = (
        z + (1.0 - z) * math.ldexp(1.0, -partial_rank)
        if partial_rank
        else 1.0
    )
    zero_start_moment = full_packet_moment**full_groups * partial_packet_moment
    live_nonzero_moment = math.expm1(state_bits * math.log1p(z)) / full_states
    matrix = np.asarray(
        (
            (zero_input, zero_start_moment - zero_input),
            (1.0 / full_states, live_nonzero_moment),
        ),
        dtype=np.float64,
    )
    with np.errstate(divide="ignore"):
        return np.log(matrix)


def deterministic_support_region_envelope(
    *,
    state_bits: int,
    packet_slots_per_epoch: int,
    epochs_per_region: int,
    groups: int,
    z: float,
) -> np.ndarray:
    """Return entrywise max_{K>=H} R_K for every packet count H."""
    full_states = math.ldexp(1.0, state_bits)
    live_nonzero_moment = math.expm1(state_bits * math.log1p(z)) / full_states
    epoch_logs = []
    for active_packets in range(packet_slots_per_epoch + 1):
        matrix = np.asarray(
            (
                (
                    1.0 if active_packets == 0 else 0.0,
                    0.0 if active_packets == 0 else z**active_packets,
                ),
                (1.0 / full_states, live_nonzero_moment),
            )
        )
        with np.errstate(divide="ignore"):
            epoch_logs.append(np.log(matrix))
    exact_by_count = regular_region_log_matrices(
        np.stack(epoch_logs),
        packet_slots_per_epoch,
        epochs_per_region,
        groups,
    )
    return np.maximum.accumulate(exact_by_count[::-1], axis=0)[::-1]


def packed_profile_region_logs(
    envelope: np.ndarray, *, groups: int, packet_bits: int
) -> np.ndarray:
    """Average the count envelope for every packed rank profile."""
    result = np.full((groups * packet_bits + 1, 2, 2), -math.inf)
    full_nonzero = 1.0 - math.ldexp(1.0, -packet_bits)
    for full_groups in range(groups + 1):
        counts = np.arange(full_groups + 1, dtype=np.float64)
        log_probabilities = (
            gammaln(full_groups + 1.0)
            - gammaln(counts + 1.0)
            - gammaln(full_groups - counts + 1.0)
            + counts * math.log(full_nonzero)
            + (full_groups - counts) * math.log1p(-full_nonzero)
        )
        ordinary = logsumexp(
            log_probabilities[:, None, None] + envelope[: full_groups + 1],
            axis=0,
        )
        result[packet_bits * full_groups] = ordinary
        if full_groups == groups:
            continue
        shifted = logsumexp(
            log_probabilities[:, None, None] + envelope[1 : full_groups + 2],
            axis=0,
        )
        for partial_rank in range(1, packet_bits):
            partial_nonzero = 1.0 - math.ldexp(1.0, -partial_rank)
            result[packet_bits * full_groups + partial_rank] = np.logaddexp(
                math.log1p(-partial_nonzero) + ordinary,
                math.log(partial_nonzero) + shifted,
            )
    return result


def support_monotonicity_self_test() -> dict[str, float | int]:
    """Check the epoch relaxation and the finite packing coupling audit."""
    from certify_riffle_transpose_packetshuffle_sparse import packing_self_test

    maximum_violation = 0.0
    for state_bits in (2, 3, 6):
        for z in (0.13, 0.57, 0.91, 1.0):
            live_states = math.ldexp(1.0, state_bits) - 1.0
            kappa = math.ldexp(1.0, state_bits) / live_states
            for packet_bits in range(1, min(4, state_bits) + 1):
                for full_groups in range(3):
                    for partial_rank in range(packet_bits):
                        candidates = packet_bits * full_groups + partial_rank
                        p = math.ldexp(1.0, -candidates)
                        q = ((1.0 + z) * 0.5) ** candidates
                        exact = np.asarray(
                            (
                                (p, q - p),
                                (
                                    (1.0 - p) / live_states,
                                    (
                                        (1.0 + z) ** state_bits
                                        - q
                                        - 1.0
                                        + p
                                    )
                                    / live_states,
                                ),
                            )
                        )
                        relaxed = np.exp(
                            packet_indicator_epoch_log_matrix(
                                state_bits,
                                packet_bits,
                                full_groups,
                                partial_rank,
                                z,
                            )
                        )
                        relaxed[1] *= kappa
                        maximum_violation = max(
                            maximum_violation,
                            float(np.max(exact - relaxed)),
                        )
    if maximum_violation > 2e-14:
        raise AssertionError("packet-indicator epoch relaxation failed")
    packing = packing_self_test()
    return {
        "maximum_epoch_upper_violation": maximum_violation,
        "packing_occupancy_vectors_checked": packing["occupancy_vectors_checked"],
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    if args.state_bits != args.epoch_bits:
        raise ValueError("this checker requires one visit per state lane")
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
    total_epochs = args.outer_bits * args.epochs_per_region
    log_live_normalization = -math.log1p(-math.ldexp(1.0, -args.state_bits))

    best = np.full(outer_blocks + 1, math.inf)
    best_tilt = np.full(outer_blocks + 1, math.nan)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        support_envelope = deterministic_support_region_envelope(
            state_bits=args.state_bits,
            packet_slots_per_epoch=packet_slots_per_epoch,
            epochs_per_region=args.epochs_per_region,
            groups=groups,
            z=z,
        )
        packed_regions = packed_profile_region_logs(
            support_envelope,
            groups=groups,
            packet_bits=args.packet_bits,
        )
        moments = log_matrix_power_moments_batch(packed_regions, args.outer_bits)
        candidates = (
            moments[1:]
            + total_epochs * log_live_normalization
            + distance * surprisal
        )
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
    return {
        "schema": "riffle-packet4-fieldcheckpoint-s256-regular-support-v1",
        "candidate": "BCHPerm-TransposePacketShuffle-FieldCheckpoint-s256-g4",
        "method": {
            "outer_reduction": "modeled regular-spectrum density envelope",
            "profile": "proved packed nonzero-packet-count upper bound",
            "epoch_transfer": "one-bit-per-nonzero-packet upper transfer",
            "live_normalization": (
                "uniform 256-bit renewal kernel times the global factor "
                "(2^256/(2^256-1))^8192"
            ),
            "region_recurrence": (
                "exact conditional packet-count matrices followed by the "
                "entrywise suffix envelope max_{K>=H} R_K"
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
            "log2_live_normalization_correction": (
                total_epochs * log_live_normalization / LOG2
            ),
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
            "special_support_excluded": "the unique all-one outer word",
        },
        "self_test": support_monotonicity_self_test(),
        "regular_log2_upper": total_log / LOG2,
        "regular_lambda_bits_lower_float": -total_log / LOG2,
        "dominant_rows": sorted(
            rows,
            key=lambda row: float(row["pointwise_log2_upper"]),
            reverse=True,
        )[:30],
        "occupation_rows": rows,
        "scope": (
            "Every regular active-block count, packet-group histogram, and "
            "zero/live trajectory.  The unique all-one outer word and outward "
            "rounding remain open."
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
        default=(-10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
                 -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0,
                 -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5,
                 0.75, 1.0),
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
