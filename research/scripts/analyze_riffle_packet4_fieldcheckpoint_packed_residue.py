#!/usr/bin/env python3
"""Full regular-spectrum diagnostic for packed packet-lane occupations.

The spectrum-density envelope replaces every regular active outer word by a
uniform 256-bit word.  For a fixed active-block count, this diagnostic packs
the active blocks into as few packet residues as possible.  A full residue
contains all 2048 packet slots.  The remaining blocks occupy a uniform subset
of the next residue's packet slots.

This is an exact floating-point calculation for that packed profile.  Using
it as a certificate for arbitrary active-block placements requires a separate
packing-domination lemma.
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


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_packetshuffle_fieldcheckpoint_g4/"
    "receipts/packed_residue_full_regular_delta09.json"
)


def lane_candidate_polynomials(visits: int, z: float) -> tuple[np.ndarray, ...]:
    """Return four lane polynomials indexed by candidate-position count.

    The four arrays respectively sum the moments for a zero starting bit, a
    zero start and zero endpoint, the unique start producing a zero endpoint,
    and both possible starting bits.  Candidate bits are independent and fair.
    """
    zero_total = np.zeros(visits + 1, dtype=np.float64)
    zero_to_zero = np.zeros(visits + 1, dtype=np.float64)
    endpoint_start = np.zeros(visits + 1, dtype=np.float64)
    all_starts = np.zeros(visits + 1, dtype=np.float64)
    for candidates in range(1 << visits):
        count = candidates.bit_count()
        bit_probability = math.ldexp(1.0, -count)
        impulses = candidates
        while True:
            parity = 0
            zero_weight = 0
            for visit in range(visits):
                parity ^= (impulses >> visit) & 1
                zero_weight += parity
            one_weight = visits - zero_weight
            delta = impulses.bit_count() & 1
            zero_moment = z**zero_weight
            zero_total[count] += bit_probability * zero_moment
            if delta == 0:
                zero_to_zero[count] += bit_probability * zero_moment
            endpoint_start[count] += bit_probability * (
                z**one_weight if delta else zero_moment
            )
            all_starts[count] += bit_probability * (
                zero_moment + z**one_weight
            )
            if impulses == 0:
                break
            impulses = (impulses - 1) & candidates
    return zero_total, zero_to_zero, endpoint_start, all_starts


def _polynomial_power(base: np.ndarray, exponent: int) -> np.ndarray:
    result = np.asarray((1.0,), dtype=np.float64)
    for _ in range(exponent):
        result = np.convolve(result, base)
    return result


def packed_candidate_epoch_logs(
    *,
    state_bits: int,
    epoch_bits: int,
    packet_bits: int,
    full_residues: int,
    z: float,
) -> np.ndarray:
    """Return log epoch transfers for 0..slots-per-residue partial candidates."""
    if state_bits % packet_bits:
        raise ValueError("packet width must divide state width")
    if epoch_bits % state_bits:
        raise ValueError("state width must divide epoch length")
    if not 0 <= full_residues < packet_bits:
        raise ValueError("a partial residue must remain")
    visits = epoch_bits // state_bits
    lanes_per_residue = state_bits // packet_bits
    partial_slots = visits * lanes_per_residue
    polynomials = lane_candidate_polynomials(visits, z)

    full_lanes = full_residues * lanes_per_residue
    partial_lanes = lanes_per_residue
    inactive_lanes = state_bits - full_lanes - partial_lanes
    live_states = math.ldexp(1.0, state_bits) - 1.0

    global_polynomials = []
    for lane_polynomial in polynomials:
        full_factor = lane_polynomial[visits] ** full_lanes
        inactive_factor = lane_polynomial[0] ** inactive_lanes
        partial = _polynomial_power(lane_polynomial, partial_lanes)
        global_polynomials.append(partial * full_factor * inactive_factor)
    zero_total, zero_to_zero, endpoint_start, all_starts = global_polynomials

    matrices = np.zeros((partial_slots + 1, 2, 2), dtype=np.float64)
    for candidates in range(partial_slots + 1):
        denominator = math.comb(partial_slots, candidates)
        m00 = zero_to_zero[candidates] / denominator
        m01 = (zero_total[candidates] - zero_to_zero[candidates]) / denominator
        m10 = (
            endpoint_start[candidates] - zero_to_zero[candidates]
        ) / (live_states * denominator)
        live_total = (
            all_starts[candidates] - zero_total[candidates]
        ) / (live_states * denominator)
        m11 = live_total - m10
        matrices[candidates] = ((m00, m01), (m10, m11))
    matrices[(matrices < 0.0) & (matrices > -3e-13)] = 0.0
    if np.any(matrices < 0.0):
        raise ArithmeticError("packed epoch transfer has a negative entry")
    with np.errstate(divide="ignore"):
        return np.log(matrices)


def full_state_epoch_matrix(
    *, state_bits: int, epoch_bits: int, packet_bits: int, z: float
) -> np.ndarray:
    """Return the epoch transfer when all four residues are full."""
    visits = epoch_bits // state_bits
    lane_polynomials = lane_candidate_polynomials(visits, z)
    live_states = math.ldexp(1.0, state_bits) - 1.0
    values = [polynomial[visits] ** state_bits for polynomial in lane_polynomials]
    zero_total, zero_to_zero, endpoint_start, all_starts = values
    m00 = zero_to_zero
    m01 = zero_total - zero_to_zero
    m10 = (endpoint_start - zero_to_zero) / live_states
    live_total = (all_starts - zero_total) / live_states
    matrix = np.asarray(((m00, m01), (m10, live_total - m10)))
    matrix[(matrix < 0.0) & (matrix > -3e-13)] = 0.0
    return matrix


def self_test(args: argparse.Namespace) -> dict[str, float]:
    maximum_stochastic_error = 0.0
    boundary_error = 0.0
    for full_residues in range(args.packet_bits):
        epoch_logs = packed_candidate_epoch_logs(
            state_bits=args.state_bits,
            epoch_bits=args.epoch_bits,
            packet_bits=args.packet_bits,
            full_residues=full_residues,
            z=1.0,
        )
        maximum_stochastic_error = max(
            maximum_stochastic_error,
            float(
                np.max(
                    np.abs(np.exp(logsumexp(epoch_logs, axis=2)) - 1.0)
                )
            ),
        )
        if full_residues + 1 < args.packet_bits:
            next_logs = packed_candidate_epoch_logs(
                state_bits=args.state_bits,
                epoch_bits=args.epoch_bits,
                packet_bits=args.packet_bits,
                full_residues=full_residues + 1,
                z=0.73,
            )
            current_logs = packed_candidate_epoch_logs(
                state_bits=args.state_bits,
                epoch_bits=args.epoch_bits,
                packet_bits=args.packet_bits,
                full_residues=full_residues,
                z=0.73,
            )
            boundary_error = max(
                boundary_error,
                float(np.max(np.abs(np.exp(current_logs[-1]) - np.exp(next_logs[0])))),
            )
    if max(maximum_stochastic_error, boundary_error) > 2e-12:
        raise AssertionError("packed-residue epoch self-test failed")
    return {
        "maximum_epoch_stochastic_error": maximum_stochastic_error,
        "maximum_full_partial_boundary_error": boundary_error,
    }


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    dimension = args.outer_bits // 2
    outer_blocks = args.message_bits // dimension
    if outer_blocks % args.packet_bits:
        raise ValueError("packet width must divide outer-block count")
    slots_per_residue_region = outer_blocks // args.packet_bits
    epoch_slots_per_residue = args.epoch_bits // args.packet_bits
    if slots_per_residue_region != epoch_slots_per_residue * args.epochs_per_region:
        raise ValueError("epoch geometry does not fill a packet residue")
    output_bits = 2 * args.message_bits
    distance = math.floor(args.relative_distance * output_bits)
    spectrum = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.modeled_minimum_distance
    )
    log_eta, envelope_weight = spectrum_density_envelope_log(
        args.outer_bits, dimension, spectrum
    )
    regular_log_mass = log_two_power_minus_one(dimension) + log_eta

    best = np.full(outer_blocks + 1, math.inf, dtype=np.float64)
    best_tilt = np.full(outer_blocks + 1, math.nan, dtype=np.float64)
    for tilt_index, log_surprisal in enumerate(args.log_surprisals):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        for full_residues in range(args.packet_bits):
            epoch_logs = packed_candidate_epoch_logs(
                state_bits=args.state_bits,
                epoch_bits=args.epoch_bits,
                packet_bits=args.packet_bits,
                full_residues=full_residues,
                z=z,
            )
            region_logs = regular_region_log_matrices(
                epoch_logs,
                epoch_slots_per_residue,
                args.epochs_per_region,
                slots_per_residue_region,
            )
            moments = log_matrix_power_moments_batch(region_logs, args.outer_bits)
            candidates = moments + distance * surprisal
            first = full_residues * slots_per_residue_region
            for partial in range(1, slots_per_residue_region + 1):
                active_blocks = first + partial
                if candidates[partial] < best[active_blocks]:
                    best[active_blocks] = candidates[partial]
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
                "full_packet_residues": (active_blocks - 1) // slots_per_residue_region,
                "partial_residue_occupation": (
                    (active_blocks - 1) % slots_per_residue_region + 1
                ),
                "best_log_surprisal": float(best_tilt[active_blocks]),
                "outer_log2_envelope": outer_log / LOG2,
                "inner_log2_upper": inner_log / LOG2,
                "pointwise_log2_upper": contribution / LOG2,
            }
        )
    total_log = float(logsumexp(np.asarray(contributions)))
    dominant = sorted(
        rows, key=lambda row: float(row["pointwise_log2_upper"]), reverse=True
    )[:30]
    return {
        "schema": "riffle-packet4-fieldcheckpoint-packed-residue-regular-v1",
        "candidate": "BCHPerm-TransposePacketShuffle-g4-FieldCheckpoint",
        "method": {
            "outer_reduction": "modeled regular-spectrum density envelope",
            "active_block_profile": (
                "pack into full packet residues, then one partial residue"
            ),
            "epoch_transfer": (
                "exact fair-candidate average over a uniform subset of the "
                "partial residue's 64 packet slots"
            ),
            "region_recurrence": "exact hypergeometric conditioning across 32 epochs",
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
        "packed_profile_log2_upper": total_log / LOG2,
        "packed_profile_lambda_bits_lower_float": -total_log / LOG2,
        "dominant_rows": dominant,
        "occupation_rows": rows,
        "scope": (
            "Complete floating-point sum for packed residue profiles with only "
            "regular outer words. A proof that packed profiles dominate every "
            "packet-group placement and the all-one cases remain open."
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
        "packed_profile_lambda_bits,"
        f"{payload['packed_profile_lambda_bits_lower_float']:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
