#!/usr/bin/env python3
"""Scalar packet-fugacity probe for a repeated nonzero eight-bit value.

Let exactly p of the M packet slots contain the same nonzero byte v and let
the other packets be zero.  The global packet permutation makes their slots
a uniform p-subset.  At one recursive step, condition on the incoming state
weight s; its random coordinate permutation is a uniform s-subset.

For every s and every subset A of the eight current packet slots, a small
exact DP computes

    E[z^wt(Acc(U_A + sigma(S)))]

without making any assumption about the next BCH state.  If F_s(x,z) sums
that moment with weight x^|A| and Lambda=max_s F_s, tower induction gives

    E[z^total_output | p repeated packets]
      <= Lambda(x,z)^B x^-p / C(M,p).

The maximization over s makes this valid for the actual BCH state evolution;
only the displayed floating evaluation and pole optimization remain
diagnostic.  Exact/outward hardening is straightforward if a row closes.
"""

from __future__ import annotations

import argparse
import math

import numpy as np


BITS = 64
PACKET_BITS = 8
PACKETS_PER_INNER = 8
INNER_BLOCKS = 32768
PACKET_SLOTS = INNER_BLOCKS * PACKETS_PER_INNER
N = BITS * INNER_BLOCKS
D = 9 * N // 100


def accumulator_histograms(current: int) -> np.ndarray:
    """Exact output-weight rows for all uniform state-support weights."""

    # One subset DP records every final selected weight simultaneously.
    # Keeping the transitions as NumPy slice additions avoids a Python loop
    # over the output-weight coordinate.
    dp = np.zeros((BITS + 1, 2, BITS + 1), dtype=np.uint64)
    dp[0, 0, 0] = 1
    for position in range(BITS):
        next_dp = np.zeros_like(dp)
        current_bit = (current >> position) & 1
        selected = slice(0, position + 1)
        selected_plus_one = slice(1, position + 2)
        for parity in range(2):
            next_parity = parity ^ current_bit
            if next_parity:
                next_dp[selected, next_parity, 1:] += dp[selected, parity, :-1]
            else:
                next_dp[selected, next_parity, :] += dp[selected, parity, :]
            next_parity ^= 1
            if next_parity:
                next_dp[selected_plus_one, next_parity, 1:] += dp[selected, parity, :-1]
            else:
                next_dp[selected_plus_one, next_parity, :] += dp[selected, parity, :]
        dp = next_dp
    histograms = dp.sum(axis=1)
    for state_weight, histogram in enumerate(histograms):
        if int(histogram.sum()) != math.comb(BITS, state_weight):
            raise SystemExit("packet8 scalar: accumulator histogram mass mismatch")
    return histograms


def build_distributions(packet_value: int) -> np.ndarray:
    distributions = np.zeros(
        (BITS + 1, 1 << PACKETS_PER_INNER, BITS + 1), dtype=np.float64
    )
    denominators = np.array(
        [float(math.comb(BITS, state)) for state in range(BITS + 1)]
    )
    currents = []
    for mask in range(1 << PACKETS_PER_INNER):
        current = 0
        for slot in range(PACKETS_PER_INNER):
            if (mask >> slot) & 1:
                current |= packet_value << (PACKET_BITS * slot)
        currents.append(current)
    for mask, current in enumerate(currents):
        distributions[:, mask] = accumulator_histograms(current)
    distributions /= denominators[:, None, None]
    if np.max(np.abs(distributions.sum(axis=2) - 1.0)) > 2e-15:
        raise SystemExit("packet8 scalar: normalized rows do not sum to one")
    return distributions


def log2_binomial(n: int, k: int) -> float:
    return (
        math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    ) / math.log(2.0)


def optimize_x(
    coefficients: np.ndarray,
    packets: int,
    iterations: int = 100,
    first_state: int = 0,
):
    def objective(log_x: float):
        x = math.exp(log_x)
        powers = np.array([x**weight for weight in range(9)])
        rows = coefficients[first_state:] @ powers
        lam = float(np.max(rows))
        return (
            INNER_BLOCKS * math.log2(lam) - packets * math.log2(x),
            x,
            lam,
            first_state + int(np.argmax(rows)),
        )

    low = -16.0
    high = 16.0
    for _ in range(iterations):
        left = (2.0 * low + high) / 3.0
        right = (low + 2.0 * high) / 3.0
        if objective(left)[0] <= objective(right)[0]:
            high = right
        else:
            low = left
    return objective((low + high) / 2.0)


def two_state_finite_bound(
    coefficients: np.ndarray, x: float, packet_value_weight: int
):
    powers = np.array([x**weight for weight in range(9)])
    rows = coefficients @ powers
    off_to_live = max(0.0, float(rows[0]) - 1.0)
    live_to_live = float(np.max(rows[1:]))
    turnoff = 0.0
    turnoff_state = 0
    for packet_count in range(1, 9):
        state = packet_value_weight * packet_count
        candidate = (
            math.comb(PACKETS_PER_INNER, packet_count)
            * x**packet_count
            / math.comb(BITS, state)
        )
        if candidate > turnoff:
            turnoff = candidate
            turnoff_state = state
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
    off_log2 = vector_log2_scale + math.log2(float(vector[0]))
    return off_log2, turnoff_state, turnoff, int(np.argmax(rows[1:])) + 1


def optimize_two_state_x(
    coefficients: np.ndarray,
    packets: int,
    packet_value_weight: int,
    iterations: int = 100,
):
    def objective(log_x: float):
        x = math.exp(log_x)
        sequence_log2, turnoff_state, turnoff, live_state = two_state_finite_bound(
            coefficients, x, packet_value_weight
        )
        value = sequence_log2 - packets * math.log2(x)
        return value, x, sequence_log2, turnoff_state, turnoff, live_state

    low = -16.0
    high = 16.0
    for _ in range(iterations):
        left = (2.0 * low + high) / 3.0
        right = (low + 2.0 * high) / 3.0
        if objective(left)[0] <= objective(right)[0]:
            high = right
        else:
            low = left
    return objective((low + high) / 2.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packet-value", type=lambda value: int(value, 0), default=0xFF)
    parser.add_argument("--active-numerator", type=int, default=21)
    parser.add_argument("--active-denominator", type=int, default=128)
    parser.add_argument(
        "--output-poles",
        default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98,0.99,0.997",
    )
    parser.add_argument("--family-dimension", type=int, default=64)
    parser.add_argument("--graph-bits", type=int, default=24)
    args = parser.parse_args()
    if not 1 <= args.packet_value <= 0xFF:
        raise SystemExit("packet8 scalar: packet value must be nonzero")
    if not 0 < args.active_numerator < args.active_denominator:
        raise SystemExit("packet8 scalar: invalid active density")
    if PACKET_SLOTS * args.active_numerator % args.active_denominator:
        raise SystemExit("packet8 scalar: active packet count is not integral")
    packets = PACKET_SLOTS * args.active_numerator // args.active_denominator
    poles = [float(value) for value in args.output_poles.split(",")]
    if any(not 0.0 < pole < 1.0 for pole in poles):
        raise SystemExit("packet8 scalar: output poles must lie in (0,1)")

    distributions = build_distributions(args.packet_value)
    packet_counts = np.array([mask.bit_count() for mask in range(256)])
    denominator_log2 = log2_binomial(PACKET_SLOTS, packets)
    best = None
    best_two = None
    print("repeated-packet scalar fugacity probe")
    print(
        f"packet_value=0x{args.packet_value:02x} active_packets={packets} "
        f"packet_slots={PACKET_SLOTS} d={D}"
    )
    for pole in poles:
        pole_powers = np.array([pole**weight for weight in range(BITS + 1)])
        moments = distributions @ pole_powers
        coefficients = np.zeros((BITS + 1, 9), dtype=np.float64)
        for packet_count in range(9):
            coefficients[:, packet_count] = moments[:, packet_counts == packet_count].sum(axis=1)
        cauchy, x, lam, state = optimize_x(coefficients, packets)
        live_cauchy, live_x, live_lam, live_state = optimize_x(
            coefficients, packets, first_state=1
        )
        probability = cauchy - denominator_log2 - D * math.log2(pole)
        live_probability = (
            live_cauchy - denominator_log2 - D * math.log2(pole)
        )
        (
            two_cauchy,
            two_x,
            two_sequence_log2,
            turnoff_state,
            turnoff,
            two_live_state,
        ) = optimize_two_state_x(
            coefficients, packets, args.packet_value.bit_count()
        )
        two_probability = (
            two_cauchy - denominator_log2 - D * math.log2(pole)
        )
        two_union = two_probability + args.family_dimension - args.graph_bits
        two_row = (
            two_union,
            pole,
            two_x,
            turnoff_state,
            two_live_state,
            two_probability,
        )
        if best_two is None or two_row < best_two:
            best_two = two_row
        union = probability + args.family_dimension - args.graph_bits
        row = (union, pole, x, lam, state, probability)
        if best is None or row < best:
            best = row
        print(
            f"pole={pole:.6f} x={x:.9f} lambda_log2={math.log2(lam):.12f} "
            f"worst_state={state} probability_log2={probability:.6f} "
            f"family_graph_union_log2={union:.6f} "
            f"live_x={live_x:.9f} live_lambda_log2={math.log2(live_lam):.12f} "
            f"live_worst_state={live_state} live_probability_log2={live_probability:.6f} "
            f"two_x={two_x:.9f} two_sequence_log2={two_sequence_log2:.12f} "
            f"two_live_state={two_live_state} turnoff_state={turnoff_state} "
            f"turnoff={turnoff:.3e} two_probability_log2={two_probability:.6f} "
            f"two_family_graph_union_log2={two_union:.6f}"
        )
    assert best is not None
    assert best_two is not None
    print(
        f"best_union_log2={best[0]:.6f} best_pole={best[1]:.6f} "
        f"best_x={best[2]:.9f} best_worst_state={best[4]}"
    )
    print(
        f"best_two_state_union_log2={best_two[0]:.6f} "
        f"best_two_state_pole={best_two[1]:.6f} "
        f"best_two_state_x={best_two[2]:.9f} "
        f"best_two_state_turnoff_state={best_two[3]} "
        f"best_two_state_live_state={best_two[4]}"
    )
    print("status=DIAGNOSTIC_EXACT_LOCAL_DP_FLOATING_POLE_OPTIMIZATION")


if __name__ == "__main__":
    main()
