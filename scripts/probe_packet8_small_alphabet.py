#!/usr/bin/env python3
"""OFF/LIVE packet-fugacity probe for a bounded concrete-byte profile.

The packet multiset contains ``A_v`` copies of each listed nonzero byte and
zeros in all remaining slots.  For positive fugacities ``x_v`` (with the zero
fugacity normalized to one), one byte contributes the parity/state-weight
matrix

    Q = T_0 + sum_v x_v T_v.

The eighth power of this byte matrix sums every ordered eight-packet current
input.  Its state-weight coefficients give the exact local accumulator MGF for
all 65 incoming support weights.  The OFF/LIVE relaxation then retains the
exact zero-state transition and an exact worst-state turnoff polynomial.

Coefficient extraction is upper-bounded by positive multivariate Cauchy and
the global packet orbit is normalized by the exact multinomial coefficient
``M! / product_v A_v!``.  Local combinatorial coefficients are exact before
the selected pole/fugacity evaluation.  Binary64 evaluation and coordinate
optimization make this a diagnostic, not a final certificate.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from probe_packet8_constant_value_scalar import BITS, D, INNER_BLOCKS, PACKET_SLOTS
from scan_packet8_repeated_values import byte_transition


PACKETS_PER_INNER = 8


def parse_counts(specification: str) -> dict[int, int]:
    counts: dict[int, int] = {}
    for item in specification.split(","):
        value_text, count_text = item.split(":", 1)
        value = int(value_text, 0)
        count = int(count_text)
        if not 1 <= value <= 0xFF or count <= 0 or value in counts:
            raise ValueError(f"invalid packet count item: {item}")
        counts[value] = count
    if sum(counts.values()) >= PACKET_SLOTS:
        raise ValueError("nonzero packet counts must leave at least one zero packet")
    return counts


def log2_multinomial(counts: list[int]) -> float:
    result = math.lgamma(sum(counts) + 1)
    result -= sum(math.lgamma(count + 1) for count in counts)
    return result / math.log(2.0)


def local_rows(
    zero_transition: np.ndarray,
    live_transitions: list[np.ndarray],
    fugacities: np.ndarray,
) -> np.ndarray:
    mixed = zero_transition.copy()
    for fugacity, transition in zip(fugacities, live_transitions):
        mixed += fugacity * transition
    dp = np.zeros((2, BITS + 1), dtype=np.float64)
    dp[0, 0] = 1.0
    for slot in range(PACKETS_PER_INNER):
        next_dp = np.zeros_like(dp)
        old_degree = 8 * slot
        for parity in range(2):
            source = dp[parity, : old_degree + 1]
            for next_parity in range(2):
                next_dp[next_parity, : old_degree + 9] += np.convolve(
                    source, mixed[parity, next_parity]
                )
        dp = next_dp
    denominators = np.array(
        [float(math.comb(BITS, state)) for state in range(BITS + 1)]
    )
    return dp.sum(axis=0) / denominators


def turnoff_upper(packet_values: list[int], fugacities: np.ndarray) -> tuple[float, int]:
    one_slot = np.zeros(BITS + 1, dtype=np.float64)
    one_slot[0] = 1.0
    for value, fugacity in zip(packet_values, fugacities):
        one_slot[value.bit_count()] += fugacity
    total = np.array([1.0])
    for _ in range(PACKETS_PER_INNER):
        total = np.convolve(total, one_slot)
    best = 0.0
    best_state = 0
    for state in range(1, BITS + 1):
        candidate = total[state] / math.comb(BITS, state)
        if candidate > best:
            best = float(candidate)
            best_state = state
    return best, best_state


def finite_sequence_log2(
    rows: np.ndarray,
    packet_values: list[int],
    fugacities: np.ndarray,
) -> tuple[float, int, int]:
    off_to_live = max(0.0, float(rows[0]) - 1.0)
    live_state = int(np.argmax(rows[1:])) + 1
    live_to_live = float(rows[live_state])
    turnoff, turnoff_state = turnoff_upper(packet_values, fugacities)
    matrix = np.array([[1.0, off_to_live], [turnoff, live_to_live]])
    matrix_scale = float(np.max(matrix))
    matrix /= matrix_scale
    matrix_log2_scale = math.log2(matrix_scale)
    vector = np.ones(2)
    vector_log2_scale = 0.0
    exponent = INNER_BLOCKS
    while exponent:
        if exponent & 1:
            vector = matrix @ vector
            scale = float(np.max(vector))
            vector /= scale
            vector_log2_scale += matrix_log2_scale + math.log2(scale)
        exponent >>= 1
        if not exponent:
            break
        matrix = matrix @ matrix
        scale = float(np.max(matrix))
        matrix /= scale
        matrix_log2_scale = 2.0 * matrix_log2_scale + math.log2(scale)
    return (
        vector_log2_scale + math.log2(float(vector[0])),
        turnoff_state,
        live_state,
    )


def optimize_fugacities(
    packet_values: list[int],
    packet_counts: np.ndarray,
    pole: float,
    passes: int = 10,
    iterations: int = 48,
) -> tuple[float, np.ndarray, float, int, int]:
    zero_transition = byte_transition(0, pole)
    live_transitions = [byte_transition(value, pole) for value in packet_values]

    def objective(log_fugacities: np.ndarray):
        fugacities = np.exp(log_fugacities)
        rows = local_rows(zero_transition, live_transitions, fugacities)
        sequence_log2, turnoff_state, live_state = finite_sequence_log2(
            rows, packet_values, fugacities
        )
        cauchy = sequence_log2 - float(
            np.dot(packet_counts, log_fugacities / math.log(2.0))
        )
        return cauchy, sequence_log2, turnoff_state, live_state

    point = np.zeros(len(packet_values), dtype=np.float64)
    for _ in range(passes):
        old = point.copy()
        for coordinate in range(len(point)):
            low = -16.0
            high = 16.0
            for _ in range(iterations):
                left = (2.0 * low + high) / 3.0
                right = (low + 2.0 * high) / 3.0
                left_point = point.copy()
                right_point = point.copy()
                left_point[coordinate] = left
                right_point[coordinate] = right
                if objective(left_point)[0] <= objective(right_point)[0]:
                    high = right
                else:
                    low = left
            point[coordinate] = (low + high) / 2.0
        if float(np.max(np.abs(point - old))) < 1e-8:
            break
    cauchy, sequence_log2, turnoff_state, live_state = objective(point)
    return cauchy, np.exp(point), sequence_log2, turnoff_state, live_state


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--counts",
        required=True,
        help="comma-separated nonzero packet counts, e.g. 0x18:21504,0x30:21504",
    )
    parser.add_argument(
        "--output-poles", default="0.6,0.7,0.8,0.9,0.95,0.98,0.99"
    )
    parser.add_argument("--family-dimension", type=int, default=128)
    parser.add_argument("--graph-bits", type=int, default=24)
    args = parser.parse_args()
    try:
        counts = parse_counts(args.counts)
    except ValueError as error:
        raise SystemExit(f"packet8 alphabet: {error}") from error
    poles = [float(value) for value in args.output_poles.split(",")]
    if any(not 0.0 < pole < 1.0 for pole in poles):
        raise SystemExit("packet8 alphabet: output poles must lie in (0,1)")

    packet_values = sorted(counts)
    packet_counts = np.array([counts[value] for value in packet_values], dtype=np.float64)
    zero_count = PACKET_SLOTS - int(packet_counts.sum())
    denominator_log2 = log2_multinomial([zero_count, *map(int, packet_counts)])
    best = None
    print("bounded packet-alphabet OFF/LIVE fugacity probe")
    print(
        f"counts={args.counts} zero_count={zero_count} packet_slots={PACKET_SLOTS} "
        f"d={D} multinomial_log2={denominator_log2:.9f}"
    )
    for pole in poles:
        cauchy, fugacities, sequence_log2, turnoff_state, live_state = (
            optimize_fugacities(packet_values, packet_counts, pole)
        )
        probability = cauchy - denominator_log2 - D * math.log2(pole)
        union = probability + args.family_dimension - args.graph_bits
        row = (
            union,
            pole,
            probability,
            fugacities,
            sequence_log2,
            turnoff_state,
            live_state,
        )
        if best is None or row[0] < best[0]:
            best = row
        print(
            f"pole={pole:.6f} probability_log2={probability:.6f} "
            f"family_graph_union_log2={union:.6f} "
            f"fugacities={','.join(f'{value:.9f}' for value in fugacities)} "
            f"turnoff_state={turnoff_state} live_state={live_state}"
        )
    assert best is not None
    print(
        f"best_union_log2={best[0]:.6f} best_pole={best[1]:.6f} "
        f"best_probability_log2={best[2]:.6f} "
        f"best_fugacities={','.join(f'{value:.9f}' for value in best[3])} "
        f"best_turnoff_state={best[5]} best_live_state={best[6]}"
    )
    print("status=DIAGNOSTIC_EXACT_LOCAL_ALPHABET_FLOATING_OPTIMIZATION")


if __name__ == "__main__":
    main()
