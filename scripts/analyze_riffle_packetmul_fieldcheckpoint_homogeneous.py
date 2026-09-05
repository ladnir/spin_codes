#!/usr/bin/env python3
"""Regular-spectrum landscape for homogeneous packet-rank profiles.

For rank k and occupation h, exactly h of the 2048 fixed groups contain k
regular active outer blocks.  The packet permutation places those groups in a
uniform h-subset of the packet slots in every region.  Independent nonzero
GF(16) multipliers randomize every packet before FieldCheckpoint.

The calculation is exact for each homogeneous profile.  Mixed group ranks and
all-one outer words remain separate proof obligations.
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
from analyze_riffle_packetmul_fieldcheckpoint_rank1 import (
    epoch_matrix,
    rank_packet_distribution,
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
    "constructions/riffle_bchperm_transpose_packetmul_fieldcheckpoint_g4/"
    "receipts/homogeneous_regular_delta09.json"
)


def cell_candidate_polynomials(
    *, packet_bits: int, visits: int, rank: int, z: float
) -> tuple[np.ndarray, ...]:
    """Return cell moments indexed by the number of active packet slots."""
    states = 1 << packet_bits
    active_distribution = rank_packet_distribution(packet_bits, rank)
    zero_total = np.zeros(visits + 1, dtype=np.float64)
    zero_to_zero = np.zeros(visits + 1, dtype=np.float64)
    endpoint_start = np.zeros(visits + 1, dtype=np.float64)
    all_starts = np.zeros(visits + 1, dtype=np.float64)

    for active_visits in range(1 << visits):
        count = active_visits.bit_count()
        transfer = np.eye(states, dtype=np.float64)
        for visit in range(visits):
            distribution = (
                active_distribution
                if (active_visits >> visit) & 1
                else np.asarray((1.0,) + (0.0,) * (states - 1))
            )
            step = np.zeros((states, states), dtype=np.float64)
            for before in range(states):
                for mask, probability in enumerate(distribution):
                    after = before ^ mask
                    step[before, after] += probability * z ** after.bit_count()
            transfer = transfer @ step
        zero_total[count] += float(transfer[0].sum())
        zero_to_zero[count] += float(transfer[0, 0])
        endpoint_start[count] += float(transfer[:, 0].sum())
        all_starts[count] += float(transfer.sum())
    return zero_total, zero_to_zero, endpoint_start, all_starts


def _polynomial_power(base: np.ndarray, exponent: int) -> np.ndarray:
    result = np.asarray((1.0,), dtype=np.float64)
    for _ in range(exponent):
        result = np.convolve(result, base)
    return result


def candidate_epoch_logs(
    *, state_bits: int, packet_bits: int, epoch_bits: int, rank: int, z: float
) -> np.ndarray:
    """Return epoch transfers for every active-packet count from zero to 64."""
    visits = epoch_bits // state_bits
    cells = state_bits // packet_bits
    packet_slots = visits * cells
    cell_polynomials = cell_candidate_polynomials(
        packet_bits=packet_bits, visits=visits, rank=rank, z=z
    )
    global_polynomials = [
        _polynomial_power(polynomial, cells)
        for polynomial in cell_polynomials
    ]
    zero_total, zero_to_zero, endpoint_start, all_starts = global_polynomials
    live_states = math.ldexp(1.0, state_bits) - 1.0
    matrices = np.zeros((packet_slots + 1, 2, 2), dtype=np.float64)
    for active_packets in range(packet_slots + 1):
        denominator = math.comb(packet_slots, active_packets)
        m00 = zero_to_zero[active_packets] / denominator
        m01 = (zero_total[active_packets] - zero_to_zero[active_packets]) / denominator
        m10 = (
            endpoint_start[active_packets] - zero_to_zero[active_packets]
        ) / (live_states * denominator)
        live_total = (
            all_starts[active_packets] - zero_total[active_packets]
        ) / (live_states * denominator)
        matrices[active_packets] = ((m00, m01), (m10, live_total - m10))
    matrices[(matrices < 0.0) & (matrices > -4e-13)] = 0.0
    if np.any(matrices < 0.0):
        raise ArithmeticError("homogeneous packet epoch transfer is negative")
    with np.errstate(divide="ignore"):
        return np.log(matrices)


def self_test(args: argparse.Namespace) -> dict[str, float]:
    stochastic_error = 0.0
    boundary_error = 0.0
    packet_slots = args.epoch_bits // args.packet_bits
    for rank in range(1, args.packet_bits + 1):
        logs = candidate_epoch_logs(
            state_bits=args.state_bits,
            packet_bits=args.packet_bits,
            epoch_bits=args.epoch_bits,
            rank=rank,
            z=1.0,
        )
        stochastic_error = max(
            stochastic_error,
            float(np.max(np.abs(np.exp(logsumexp(logs, axis=2)) - 1.0))),
        )
        full = np.exp(
            candidate_epoch_logs(
                state_bits=args.state_bits,
                packet_bits=args.packet_bits,
                epoch_bits=args.epoch_bits,
                rank=rank,
                z=0.73,
            )[packet_slots]
        )
        reference = epoch_matrix(
            state_bits=args.state_bits,
            packet_bits=args.packet_bits,
            epoch_bits=args.epoch_bits,
            rank=rank,
            z=0.73,
        )
        boundary_error = max(
            boundary_error, float(np.max(np.abs(full - reference)))
        )
    if max(stochastic_error, boundary_error) > 5e-12:
        raise AssertionError("homogeneous profile transfer self-test failed")
    return {
        "maximum_epoch_stochastic_error": stochastic_error,
        "maximum_full_occupation_reference_error": boundary_error,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    groups = outer_blocks // args.packet_bits
    packet_slots_per_epoch = args.epoch_bits // args.packet_bits
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

    rank_results = []
    aggregate_terms = []
    for rank in range(1, args.packet_bits + 1):
        best = np.full(groups + 1, math.inf)
        best_tilt = np.full(groups + 1, math.nan)
        for log_surprisal in args.log_surprisals:
            surprisal = math.exp(log_surprisal)
            z = math.exp(-surprisal)
            epoch_logs = candidate_epoch_logs(
                state_bits=args.state_bits,
                packet_bits=args.packet_bits,
                epoch_bits=args.epoch_bits,
                rank=rank,
                z=z,
            )
            region_logs = regular_region_log_matrices(
                epoch_logs,
                packet_slots_per_epoch,
                args.epochs_per_region,
                groups,
            )
            moments = log_matrix_power_moments_batch(region_logs, args.outer_bits)
            candidates = moments + distance * surprisal
            improved = candidates < best
            best[improved] = candidates[improved]
            best_tilt[improved] = log_surprisal

        rows = []
        terms = []
        for active_groups in range(1, groups + 1):
            outer_log = (
                log_choose(groups, active_groups)
                + active_groups * log_choose(args.packet_bits, rank)
                + rank * active_groups * regular_log_mass
            )
            inner_log = min(0.0, float(best[active_groups]))
            contribution = outer_log + inner_log
            terms.append(contribution)
            aggregate_terms.append(contribution)
            rows.append(
                {
                    "active_groups": active_groups,
                    "active_regular_outer_blocks": rank * active_groups,
                    "best_log_surprisal": float(best_tilt[active_groups]),
                    "outer_log2_envelope": outer_log / LOG2,
                    "inner_log2_upper": inner_log / LOG2,
                    "pointwise_log2_upper": contribution / LOG2,
                }
            )
        rank_log = float(logsumexp(np.asarray(terms)))
        rank_results.append(
            {
                "group_rank": rank,
                "homogeneous_log2_upper": rank_log / LOG2,
                "homogeneous_lambda_bits_lower_float": -rank_log / LOG2,
                "dominant_rows": sorted(
                    rows,
                    key=lambda row: float(row["pointwise_log2_upper"]),
                    reverse=True,
                )[:20],
                "occupation_rows": rows,
            }
        )
        print(
            f"rank,{rank},lambda_bits,{-rank_log / LOG2:.12f}", flush=True
        )

    aggregate_log = float(logsumexp(np.asarray(aggregate_terms)))
    return {
        "schema": "riffle-packetmul-fieldcheckpoint-homogeneous-regular-v1",
        "candidate": "BCHPerm-TransposePacketShuffle-PacketMul-FieldCheckpoint-g4",
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
        },
        "self_test": self_test(args),
        "homogeneous_union_log2_upper": aggregate_log / LOG2,
        "homogeneous_union_lambda_bits_lower_float": -aggregate_log / LOG2,
        "rank_results": rank_results,
        "scope": (
            "Complete floating-point sums for the four homogeneous regular "
            "group-rank families. Mixed-rank groups, all-one outer words, and "
            "outward rounding remain open."
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
        default=(-8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 0.75),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "homogeneous_union_lambda_bits,"
        f"{payload['homogeneous_union_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
