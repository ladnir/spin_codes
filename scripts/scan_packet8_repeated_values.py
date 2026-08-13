#!/usr/bin/env python3
"""Exhaustive diagnostic scan of all repeated nonzero packet values.

For one concrete byte ``v``, the repeated-value probe needs

    F_s(x,z) = sum_A x^|A| E[z^wt(Acc(U_A + sigma(S_s)))],

where ``A`` ranges over the subsets of the eight packet slots and ``S_s`` is
a uniform weight-``s`` support.  The original scalar probe constructs the
complete output-weight histogram independently for all 256 masks.  This scan
computes the same coefficients with an eight-byte trellis: a two-state byte
matrix tracks accumulator parity, a ``y`` coefficient tracks state weight,
and an ``x`` coefficient tracks the number of active packet slots.

The byte trellis uses floating evaluations of the output pole, so this remains
a diagnostic.  ``--cross-check`` compares selected rows against the original
exact-histogram implementation before the exhaustive scan.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

from probe_packet8_constant_value_scalar import (
    BITS,
    D,
    PACKET_SLOTS,
    build_distributions,
    log2_binomial,
    optimize_two_state_x,
)


PACKET_BITS = 8
PACKETS_PER_INNER = 8
DEFAULT_POLES = "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98,0.99,0.997"


def byte_transition(packet_value: int, pole: float) -> np.ndarray:
    """Return T[start_parity,end_parity,state_weight] for one byte."""

    transition = np.zeros((2, 2, PACKET_BITS + 1), dtype=np.float64)
    for start in range(2):
        dp = np.zeros((2, PACKET_BITS + 1), dtype=np.float64)
        dp[start, 0] = 1.0
        for position in range(PACKET_BITS):
            current_bit = (packet_value >> position) & 1
            next_dp = np.zeros_like(dp)
            for parity in range(2):
                # The state bit is absent.
                next_parity = parity ^ current_bit
                next_dp[next_parity] += dp[parity] * (
                    pole if next_parity else 1.0
                )
                # The state bit is present and contributes one to its weight.
                next_parity ^= 1
                next_dp[next_parity, 1:] += dp[parity, :-1] * (
                    pole if next_parity else 1.0
                )
            dp = next_dp
        transition[start] = dp
    return transition


def trellis_coefficients(packet_value: int, pole: float) -> np.ndarray:
    """Return the same 65-by-9 coefficient table as the scalar probe."""

    inactive = byte_transition(0, pole)
    active = byte_transition(packet_value, pole)
    # dp[parity, active_packet_count, selected_state_weight]
    dp = np.zeros((2, PACKETS_PER_INNER + 1, BITS + 1), dtype=np.float64)
    dp[0, 0, 0] = 1.0
    for slot in range(PACKETS_PER_INNER):
        next_dp = np.zeros_like(dp)
        old_degree = PACKET_BITS * slot
        for parity in range(2):
            for packet_count in range(slot + 1):
                source = dp[parity, packet_count, : old_degree + 1]
                if not np.any(source):
                    continue
                for next_parity in range(2):
                    zero_row = inactive[parity, next_parity]
                    live_row = active[parity, next_parity]
                    next_dp[
                        next_parity, packet_count, : old_degree + PACKET_BITS + 1
                    ] += np.convolve(source, zero_row)
                    next_dp[
                        next_parity,
                        packet_count + 1,
                        : old_degree + PACKET_BITS + 1,
                    ] += np.convolve(source, live_row)
        dp = next_dp
    numerators = dp.sum(axis=0).T
    denominators = np.array(
        [float(math.comb(BITS, state)) for state in range(BITS + 1)]
    )
    coefficients = numerators / denominators[:, None]
    return coefficients


def legacy_coefficients(packet_value: int, pole: float) -> np.ndarray:
    distributions = build_distributions(packet_value)
    packet_counts = np.array([mask.bit_count() for mask in range(256)])
    pole_powers = np.array([pole**weight for weight in range(BITS + 1)])
    moments = distributions @ pole_powers
    coefficients = np.zeros((BITS + 1, 9), dtype=np.float64)
    for packet_count in range(9):
        coefficients[:, packet_count] = moments[
            :, packet_counts == packet_count
        ].sum(axis=1)
    return coefficients


def cross_check(values: list[int], poles: list[float], tolerance: float) -> None:
    for packet_value in values:
        for pole in poles:
            fast = trellis_coefficients(packet_value, pole)
            legacy = legacy_coefficients(packet_value, pole)
            error = float(np.max(np.abs(fast - legacy)))
            scale = float(np.max(np.abs(legacy)))
            relative = error / max(1.0, scale)
            if relative > tolerance:
                raise SystemExit(
                    "packet8 repeated scan: trellis cross-check failed for "
                    f"0x{packet_value:02x} at pole {pole}: {relative:.3e}"
                )
            print(
                f"cross_check_value=0x{packet_value:02x} pole={pole:.6f} "
                f"relative_error={relative:.3e} PASS"
            )


def scan_value(
    packet_value: int,
    poles: list[float],
    packets: int,
    denominator_log2: float,
    family_dimension: int,
    graph_bits: int,
) -> dict[str, float | int | str]:
    best = None
    for pole in poles:
        coefficients = trellis_coefficients(packet_value, pole)
        (
            cauchy,
            x,
            sequence_log2,
            turnoff_state,
            _turnoff,
            live_state,
        ) = optimize_two_state_x(
            coefficients, packets, packet_value.bit_count()
        )
        probability = cauchy - denominator_log2 - D * math.log2(pole)
        union = probability + family_dimension - graph_bits
        row = (
            union,
            pole,
            x,
            probability,
            sequence_log2,
            turnoff_state,
            live_state,
        )
        if best is None or row < best:
            best = row
    assert best is not None
    return {
        "packet_value": f"0x{packet_value:02x}",
        "packet_weight": packet_value.bit_count(),
        "best_union_log2": best[0],
        "best_pole": best[1],
        "best_x": best[2],
        "probability_log2": best[3],
        "sequence_log2": best[4],
        "turnoff_state": best[5],
        "live_state": best[6],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--active-numerator", type=int, default=21)
    parser.add_argument("--active-denominator", type=int, default=128)
    parser.add_argument("--output-poles", default=DEFAULT_POLES)
    parser.add_argument("--family-dimension", type=int, default=64)
    parser.add_argument("--graph-bits", type=int, default=24)
    parser.add_argument("--output-csv", type=Path)
    parser.add_argument(
        "--cross-check",
        action="store_true",
        help="compare 0x01, 0x03, 0x55, and 0xff against the legacy exact DP",
    )
    args = parser.parse_args()
    if not 0 < args.active_numerator < args.active_denominator:
        raise SystemExit("packet8 repeated scan: invalid active density")
    if PACKET_SLOTS * args.active_numerator % args.active_denominator:
        raise SystemExit("packet8 repeated scan: active packet count is not integral")
    poles = [float(value) for value in args.output_poles.split(",")]
    if any(not 0.0 < pole < 1.0 for pole in poles):
        raise SystemExit("packet8 repeated scan: poles must lie in (0,1)")
    if args.cross_check:
        # Two well-separated poles exercise both long- and short-output rows.
        check_poles = [poles[0], poles[-1]] if len(poles) > 1 else poles
        cross_check([0x01, 0x03, 0x55, 0xFF], check_poles, 3e-13)

    packets = PACKET_SLOTS * args.active_numerator // args.active_denominator
    denominator_log2 = log2_binomial(PACKET_SLOTS, packets)
    rows = []
    for packet_value in range(1, 256):
        row = scan_value(
            packet_value,
            poles,
            packets,
            denominator_log2,
            args.family_dimension,
            args.graph_bits,
        )
        rows.append(row)
        print(
            f"value={row['packet_value']} weight={row['packet_weight']} "
            f"best_union_log2={row['best_union_log2']:.6f} "
            f"pole={row['best_pole']:.6f} x={row['best_x']:.9f}"
        )

    rows.sort(key=lambda row: float(row["best_union_log2"]), reverse=True)
    if args.output_csv is not None:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print("worst_repeated_packet_rows")
    for rank, row in enumerate(rows[:16], 1):
        print(
            f"rank={rank} value={row['packet_value']} weight={row['packet_weight']} "
            f"best_union_log2={row['best_union_log2']:.6f} "
            f"pole={row['best_pole']:.6f} x={row['best_x']:.9f} "
            f"turnoff_state={row['turnoff_state']} live_state={row['live_state']}"
        )
    print("all_255_repeated_values_scanned=PASS")
    print("status=DIAGNOSTIC_BYTE_TRELLIS_FLOATING_POLE_OPTIMIZATION")


if __name__ == "__main__":
    main()
