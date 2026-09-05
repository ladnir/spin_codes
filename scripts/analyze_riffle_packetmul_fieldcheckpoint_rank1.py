#!/usr/bin/env python3
"""Canonical medium-density test for packet multiplication before FieldCheckpoint.

The active-block set contains one block in every fixed packet group and uses
the same logical packet lane in all groups.  Thus, every region contains 2048
fair rank-one packets.  Before the checkpoint accumulator, each packet is
multiplied by an independent uniform nonzero element of GF(16).

This test exactly averages the packet multipliers, fair outer bits, and field
checkpoint multipliers.  It addresses the canonical aligned-lane obstruction;
it is not a complete construction certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_fieldcheckpoint_regular_envelope import (
    spectrum_density_envelope_log,
)
from analyze_riffle_spectrumperm_bitshuffle_oneblock import (
    modeled_even_floor_spectrum_logs,
)
from analyze_riffle_striped_random_outer import LOG2, log_two_power_minus_one


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_packetmul_fieldcheckpoint_g4/"
    "receipts/canonical_rank1_full_groups_delta09.json"
)


def rank_packet_distribution(packet_bits: int, rank: int) -> np.ndarray:
    """Return the packet mask law after a uniform nonzero field multiplier."""
    if not 1 <= rank <= packet_bits:
        raise ValueError("rank must be between one and the packet width")
    mask_count = 1 << packet_bits
    values = np.empty(mask_count, dtype=np.float64)
    values[0] = math.ldexp(1.0, -rank)
    values[1:] = (1.0 - values[0]) / (mask_count - 1)
    return values


def cell_transfer(
    *, packet_bits: int, visits: int, rank: int, z: float
) -> np.ndarray:
    """Return the exact transfer for one four-lane cell across one epoch."""
    distribution = rank_packet_distribution(packet_bits, rank)
    states = 1 << packet_bits
    transfer = np.zeros((states, states), dtype=np.float64)
    for start in range(states):
        current = np.zeros(states, dtype=np.float64)
        current[start] = 1.0
        for _ in range(visits):
            updated = np.zeros(states, dtype=np.float64)
            for before in range(states):
                if current[before] == 0.0:
                    continue
                for mask, probability in enumerate(distribution):
                    after = before ^ mask
                    updated[after] += (
                        current[before] * probability * z ** after.bit_count()
                    )
            current = updated
        transfer[start] = current
    return transfer


def epoch_matrix(
    *, state_bits: int, packet_bits: int, epoch_bits: int, rank: int, z: float
) -> np.ndarray:
    """Return the two-state epoch transfer when every packet has fixed rank."""
    if state_bits % packet_bits:
        raise ValueError("packet width must divide state width")
    if epoch_bits % state_bits:
        raise ValueError("state width must divide epoch length")
    visits = epoch_bits // state_bits
    cells = state_bits // packet_bits
    cell = cell_transfer(
        packet_bits=packet_bits, visits=visits, rank=rank, z=z
    )
    zero_total_cell = float(cell[0].sum())
    zero_to_zero_cell = float(cell[0, 0])
    endpoint_start_cell = float(cell[:, 0].sum())
    all_starts_cell = float(cell.sum())

    zero_total = zero_total_cell**cells
    zero_to_zero = zero_to_zero_cell**cells
    endpoint_start = endpoint_start_cell**cells
    all_starts = all_starts_cell**cells
    live_states = math.ldexp(1.0, state_bits) - 1.0
    m00 = zero_to_zero
    m01 = zero_total - zero_to_zero
    m10 = (endpoint_start - zero_to_zero) / live_states
    live_total = (all_starts - zero_total) / live_states
    matrix = np.asarray(((m00, m01), (m10, live_total - m10)))
    matrix[(matrix < 0.0) & (matrix > -3e-13)] = 0.0
    if np.any(matrix < 0.0):
        raise ArithmeticError("packet-multiplied epoch matrix has a negative entry")
    return matrix


def matrix_power_log_moment(matrix: np.ndarray, power: int) -> float:
    row = np.asarray((1.0, 0.0), dtype=np.float64)
    accumulated = 0.0
    for _ in range(power):
        row = row @ matrix
        scale = float(np.max(row))
        if scale <= 0.0:
            return -math.inf
        row /= scale
        accumulated += math.log(scale)
    return accumulated + math.log(float(row.sum()))


def self_test(args: argparse.Namespace) -> dict[str, float]:
    probability_error = max(
        abs(float(rank_packet_distribution(args.packet_bits, rank).sum()) - 1.0)
        for rank in range(1, args.packet_bits + 1)
    )
    stochastic_error = max(
        float(
            np.max(
                np.abs(
                    epoch_matrix(
                        state_bits=args.state_bits,
                        packet_bits=args.packet_bits,
                        epoch_bits=args.epoch_bits,
                        rank=rank,
                        z=1.0,
                    ).sum(axis=1)
                    - 1.0
                )
            )
        )
        for rank in range(1, args.packet_bits + 1)
    )
    if max(probability_error, stochastic_error) > 4e-13:
        raise AssertionError("packet multiplier transfer self-test failed")
    return {
        "maximum_packet_probability_error": probability_error,
        "maximum_epoch_stochastic_error": stochastic_error,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    groups = outer_blocks // args.packet_bits
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    epochs = output_bits // args.epoch_bits
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(dimension) + log_eta
    # There are four aligned-lane choices.  Each activates one block in every
    # group, so no block-placement binomial factor appears.
    outer_log = math.log(args.packet_bits) + groups * regular_log_mass

    best = math.inf
    best_tilt = math.nan
    rows = []
    for log_surprisal in args.log_surprisals:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        matrix = epoch_matrix(
            state_bits=args.state_bits,
            packet_bits=args.packet_bits,
            epoch_bits=args.epoch_bits,
            rank=1,
            z=z,
        )
        moment = matrix_power_log_moment(matrix, epochs)
        candidate = moment + distance * surprisal
        rows.append(
            {
                "log_surprisal": log_surprisal,
                "log2_moment": moment / LOG2,
                "inner_log2_upper": candidate / LOG2,
            }
        )
        if candidate < best:
            best = candidate
            best_tilt = log_surprisal
    pointwise = outer_log + min(0.0, best)
    return {
        "schema": "riffle-packetmul-fieldcheckpoint-rank1-full-groups-v1",
        "candidate": (
            "BCHPerm-TransposePacketShuffle-PacketMul-FieldCheckpoint-g4"
        ),
        "canonical_profile": (
            "one regular active block in the same logical lane of every "
            "four-block group"
        ),
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
            "epochs": epochs,
        },
        "envelope": {
            "log2_eta": log_eta / LOG2,
            "maximizing_regular_weight": envelope_weight,
        },
        "self_test": self_test(args),
        "best_log_surprisal": best_tilt,
        "outer_log2_envelope_for_profile": outer_log / LOG2,
        "inner_log2_upper": min(0.0, best) / LOG2,
        "pointwise_log2_upper": pointwise / LOG2,
        "lambda_bits_lower_float": -pointwise / LOG2,
        "tilt_rows": rows,
        "scope": (
            "Exact floating-point diagnostic for the canonical aligned-lane "
            "profile under the modeled regular-spectrum envelope. Other group "
            "profiles, all-one words, and outward rounding remain open."
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
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=(-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 0.75, 1.0),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"lambda_bits,{payload['lambda_bits_lower_float']:.12f}")
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
