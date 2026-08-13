#!/usr/bin/env python3
"""OFF/LIVE transfer conditional on an exact packet-weight profile.

Let ``a_j`` be the number of eight-bit packets of weight ``j``.  Conditional
on this profile, the global packet permutation orders the weights uniformly,
and the independent within-packet lane orders make every concrete value of
weight ``j`` uniform among ``C(8,j)`` choices.  Thus the exact normalization is

    M! / product_j a_j! * product_j C(8,j)^a_j.

One fugacity is assigned to each packet-weight class.  The byte transition for
class ``j`` is the sum over all concrete bytes of weight ``j``.  The same
OFF/LIVE construction as the bounded-alphabet probe then gives a conditional
low-output bound without assuming that a 64-bit input is uniform on a Hamming
slice.

The local trellis coefficients are combinatorially exact before evaluation.
Binary64 poles, multivariate Cauchy optimization, and the OFF/LIVE relaxation
make the reported bound diagnostic.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from probe_packet8_constant_value_scalar import BITS, D, PACKET_SLOTS
from probe_packet8_small_alphabet import log2_multinomial
from scan_packet8_repeated_values import byte_transition


PACKET_BITS = 8
CLASSES = tuple(math.comb(PACKET_BITS, weight) for weight in range(PACKET_BITS + 1))


def parse_profile(specification: str) -> list[int]:
    profile = [int(value) for value in specification.split(",")]
    if len(profile) != PACKET_BITS + 1 or any(value < 0 for value in profile):
        raise ValueError("profile must contain nine nonnegative counts")
    if sum(profile) != PACKET_SLOTS:
        raise ValueError(f"profile counts must sum to {PACKET_SLOTS}")
    return profile


def class_transitions(pole: float) -> list[np.ndarray]:
    result = []
    for weight in range(PACKET_BITS + 1):
        total = np.zeros((2, 2, PACKET_BITS + 1), dtype=np.float64)
        for value in range(1 << PACKET_BITS):
            if value.bit_count() == weight:
                total += byte_transition(value, pole)
        result.append(total)
    return result


def class_turnoff_upper(fugacities: np.ndarray) -> tuple[float, int]:
    one_slot = np.zeros(BITS + 1, dtype=np.float64)
    for weight, classes in enumerate(CLASSES):
        one_slot[weight] = classes * fugacities[weight]
    total = np.array([1.0])
    for _ in range(8):
        total = np.convolve(total, one_slot)
    best = 0.0
    best_state = 0
    for state in range(1, BITS + 1):
        candidate = total[state] / math.comb(BITS, state)
        if candidate > best:
            best = float(candidate)
            best_state = state
    return best, best_state


def class_local_rows(
    transitions: list[np.ndarray], fugacities: np.ndarray
) -> np.ndarray:
    mixed = np.zeros_like(transitions[0])
    for fugacity, transition in zip(fugacities, transitions):
        if fugacity:
            mixed += fugacity * transition
    dp = np.zeros((2, BITS + 1), dtype=np.float64)
    dp[0, 0] = 1.0
    for slot in range(8):
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


def class_finite_sequence_log2(
    rows: np.ndarray, fugacities: np.ndarray
) -> tuple[float, int, int]:
    # finite_sequence_log2 needs a concrete-value list only for its turnoff
    # helper.  Reproduce the small two-state wrapper with the class turnoff.
    off_to_off = float(fugacities[0] ** 8)
    off_to_live = max(0.0, float(rows[0]) - off_to_off)
    live_state = int(np.argmax(rows[1:])) + 1
    live_to_live = float(rows[live_state])
    turnoff, turnoff_state = class_turnoff_upper(fugacities)
    matrix = np.array([[off_to_off, off_to_live], [turnoff, live_to_live]])
    matrix_scale = float(np.max(matrix))
    matrix /= matrix_scale
    matrix_log2_scale = math.log2(matrix_scale)
    vector = np.ones(2)
    vector_log2_scale = 0.0
    exponent = 32768
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


def optimize_profile(
    profile: list[int], pole: float, passes: int = 12, iterations: int = 48
):
    transitions = class_transitions(pole)
    active_classes = [weight for weight in range(9) if profile[weight]]
    reference = active_classes[0]
    variable_classes = [weight for weight in active_classes if weight != reference]
    counts = np.array([profile[weight] for weight in variable_classes], dtype=np.float64)

    def objective(point: np.ndarray):
        fugacities = np.zeros(9, dtype=np.float64)
        fugacities[reference] = 1.0
        fugacities[variable_classes] = np.exp(point)
        rows = class_local_rows(transitions, fugacities)
        sequence_log2, turnoff_state, live_state = class_finite_sequence_log2(
            rows, fugacities
        )
        cauchy = sequence_log2 - float(np.dot(counts, point / math.log(2.0)))
        return cauchy, sequence_log2, turnoff_state, live_state, fugacities

    def descend(start: np.ndarray):
        point = start.copy()
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
            if len(point) == 0 or float(np.max(np.abs(point - old))) < 1e-8:
                break
        return objective(point)

    starts = [np.zeros(len(variable_classes), dtype=np.float64)]
    # For the unconstrained sequence enumerator, the saddle point matches the
    # empirical per-concrete-value frequencies.  Normalize it by the selected
    # positive-count reference class; this also handles profiles with a_0=0.
    reference_rate = profile[reference] / CLASSES[reference]
    empirical = np.array(
        [
            math.log(
                max(
                    1e-300,
                    (profile[weight] / CLASSES[weight]) / reference_rate,
                )
            )
            for weight in variable_classes
        ],
        dtype=np.float64,
    )
    starts.append(empirical)
    candidates = [descend(start) for start in starts]
    return min(candidates, key=lambda row: row[0])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--output-poles", default="0.6,0.7,0.8,0.9,0.95,0.98")
    parser.add_argument("--family-bits", type=int, default=1 << 20)
    args = parser.parse_args()
    try:
        profile = parse_profile(args.profile)
    except ValueError as error:
        raise SystemExit(f"packet8 weight profile: {error}") from error
    poles = [float(value) for value in args.output_poles.split(",")]
    if any(not 0.0 < pole < 1.0 for pole in poles):
        raise SystemExit("packet8 weight profile: poles must lie in (0,1)")

    denominator_log2 = log2_multinomial(profile) + sum(
        count * math.log2(classes)
        for count, classes in zip(profile, CLASSES)
        if count
    )
    best = None
    print("packet-weight-profile OFF/LIVE fugacity probe")
    print(
        f"profile={','.join(map(str, profile))} normalization_log2="
        f"{denominator_log2:.9f} d={D} family_bits={args.family_bits}"
    )
    for pole in poles:
        cauchy, sequence_log2, turnoff_state, live_state, fugacities = (
            optimize_profile(profile, pole)
        )
        probability = cauchy - denominator_log2 - D * math.log2(pole)
        union = probability + args.family_bits
        row = (union, pole, probability, fugacities, turnoff_state, live_state)
        if best is None or row[0] < best[0]:
            best = row
        print(
            f"pole={pole:.6f} probability_log2={probability:.6f} "
            f"family_union_log2={union:.6f} "
            f"fugacities={','.join(f'{x:.7g}' for x in fugacities)} "
            f"turnoff_state={turnoff_state} live_state={live_state}"
        )
    assert best is not None
    print(
        f"best_union_log2={best[0]:.6f} best_pole={best[1]:.6f} "
        f"best_probability_log2={best[2]:.6f} "
        f"best_turnoff_state={best[4]} best_live_state={best[5]}"
    )
    print("status=DIAGNOSTIC_EXACT_LOCAL_WEIGHT_PROFILE_FLOATING_OPTIMIZATION")


if __name__ == "__main__":
    main()
