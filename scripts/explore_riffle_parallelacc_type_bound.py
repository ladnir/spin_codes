#!/usr/bin/env python3
"""Explore a coefficient/type bound for the four-lane accumulator.

The bound retains the uniform global packet permutation.  It replaces exact
multivariate coefficient extraction by positive fugacities and evaluates the
resulting 5-by-5 transfer matrix by scaled binary exponentiation.

The first experiment bounds the complete zero-parity, weight-22,
occupation-one shell from Goal 17.  That shell has only 136 packet types, so
we can compare the bound against the committed exact shell calculation before
using the method on a complete outer spectrum.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.special import logsumexp


REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = (
    REPOSITORY_ROOT
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "receipts"
    / "goal18_packet_type_bound_probe.json"
)

PACKET_WIDTH = 4
PACKETS_PER_BLOCK = 32
BCH_LENGTH = 128
BCH_WEIGHT = 22
BCH_WEIGHT22_WORDS = 243_840
DATA_BLOCKS = 1 << 14
PACKET_POSITIONS = DATA_BLOCKS * PACKETS_PER_BLOCK
MINIMUM_SHELL_WORDS = DATA_BLOCKS * BCH_WEIGHT22_WORDS


def compositions(total: int, parts: int, prefix: tuple[int, ...] = ()):
    if parts == 1:
        yield prefix + (total,)
        return
    for first in range(total + 1):
        yield from compositions(total - first, parts - 1, prefix + (first,))


def transition_multiplicity(old: int, new: int, packet: int) -> int:
    numerator = old + new - packet
    if numerator & 1:
        return 0
    overlap = numerator // 2
    added = new - overlap
    if not 0 <= overlap <= old or not 0 <= added <= PACKET_WIDTH - old:
        return 0
    return math.comb(old, overlap) * math.comb(PACKET_WIDTH - old, added)


TRANSITIONS = np.array(
    [
        [
            [transition_multiplicity(old, new, packet) for new in range(5)]
            for old in range(5)
        ]
        for packet in range(5)
    ],
    dtype=np.float64,
)


def scaled_log_matrix_power_mass(matrix: np.ndarray, exponent: int) -> float:
    """Return log(e_0^T matrix^exponent 1) without overflow."""
    vector = np.zeros(5, dtype=np.float64)
    vector[0] = 1.0
    vector_log_scale = 0.0

    power = matrix.copy()
    power_norm = float(np.max(power))
    if not power_norm > 0.0:
        return -math.inf
    power /= power_norm
    power_log_scale = math.log(power_norm)

    remaining = exponent
    while remaining:
        if remaining & 1:
            vector = vector @ power
            vector_norm = float(np.max(vector))
            if not vector_norm > 0.0:
                return -math.inf
            vector /= vector_norm
            vector_log_scale += power_log_scale + math.log(vector_norm)
        remaining >>= 1
        if not remaining:
            break
        power = power @ power
        power_norm = float(np.max(power))
        power /= power_norm
        power_log_scale = 2.0 * power_log_scale + math.log(power_norm)

    return vector_log_scale + math.log(float(np.sum(vector)))


def local_weight22_types() -> list[tuple[tuple[int, ...], float]]:
    """Return each local histogram and log of its exact permutation mass."""
    result = []
    log_slice = math.lgamma(BCH_LENGTH + 1) - math.lgamma(BCH_WEIGHT + 1)
    log_slice -= math.lgamma(BCH_LENGTH - BCH_WEIGHT + 1)
    for histogram in compositions(PACKETS_PER_BLOCK, 5):
        if sum(weight * count for weight, count in enumerate(histogram)) != BCH_WEIGHT:
            continue
        log_count = math.lgamma(PACKETS_PER_BLOCK + 1)
        log_count -= sum(math.lgamma(count + 1) for count in histogram)
        log_count += sum(
            count * math.log(math.comb(PACKET_WIDTH, weight))
            for weight, count in enumerate(histogram)
            if count
        )
        result.append((histogram, log_count - log_slice))
    return result


LOCAL_TYPES = local_weight22_types()


def global_histogram(local: tuple[int, ...]) -> tuple[int, ...]:
    return (
        PACKET_POSITIONS - PACKETS_PER_BLOCK + local[0],
        local[1],
        local[2],
        local[3],
        local[4],
    )


def log_sequence_count(histogram: tuple[int, ...]) -> float:
    value = math.lgamma(PACKET_POSITIONS + 1)
    value -= sum(math.lgamma(count + 1) for count in histogram)
    value += sum(
        count * math.log(math.comb(PACKET_WIDTH, weight))
        for weight, count in enumerate(histogram)
        if count
    )
    return value


TYPE_ROWS = tuple(
    (
        global_histogram(local),
        log_probability,
        log_sequence_count(global_histogram(local)),
    )
    for local, log_probability in LOCAL_TYPES
)


def objective(parameters: np.ndarray, distance: int) -> float:
    """Natural log of the coefficient-bound expected bad shell count."""
    log_eta = float(parameters[0])
    log_z = -math.exp(log_eta)
    log_x = np.empty(5, dtype=np.float64)
    log_x[0] = 0.0
    log_x[1:] = parameters[1:]

    # M(x,z)_{a,b}=z^b sum_k c(a,b,k)x_k.
    matrix = np.zeros((5, 5), dtype=np.float64)
    for packet in range(5):
        matrix += math.exp(float(log_x[packet])) * TRANSITIONS[packet]
    matrix *= np.exp(log_z * np.arange(5, dtype=np.float64))[None, :]
    log_transfer_mass = scaled_log_matrix_power_mass(matrix, PACKET_POSITIONS)

    type_terms = []
    for histogram, log_probability, log_total_sequences in TYPE_ROWS:
        tilt = sum(histogram[k] * float(log_x[k]) for k in range(5))
        type_terms.append(log_probability - tilt - log_total_sequences)

    return (
        math.log(MINIMUM_SHELL_WORDS)
        - distance * log_z
        + log_transfer_mass
        + float(logsumexp(type_terms))
    )


def optimize_distance(distance: int, global_search: bool) -> dict[str, object]:
    bounds = [(-20.0, -0.5)] + [(-25.0, 25.0)] * 4
    starts = [
        np.array([-12.0, -12.0, -12.0, -12.0, -12.0]),
        np.array([-10.0, -10.0, -10.0, -10.0, -10.0]),
        np.array([-14.0, -14.0, -14.0, -14.0, -14.0]),
        np.array([-12.0, -10.0, -11.0, -12.0, -13.0]),
    ]
    candidates = []
    if global_search:
        global_result = differential_evolution(
            lambda point: objective(point, distance),
            bounds,
            seed=1,
            popsize=12,
            maxiter=80,
            polish=False,
            workers=1,
            updating="immediate",
        )
        starts.append(global_result.x)
    for start in starts:
        result = minimize(
            lambda point: objective(point, distance),
            start,
            method="Nelder-Mead",
            options={"maxiter": 4000, "xatol": 1e-9, "fatol": 1e-9},
        )
        candidates.append(result)
    best = min(candidates, key=lambda result: float(result.fun))
    log_expected = float(best.fun)
    return {
        "bad_output_weight_inclusive": distance,
        "relative_binary_output_weight": distance / (4 * PACKET_POSITIONS),
        "log2_expected_bad_shell_words_upper_bound": log_expected / math.log(2.0),
        "first_moment_margin_bits_lower_bound": -log_expected / math.log(2.0),
        "passes_shell_first_moment_under_bound": log_expected < 0.0,
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
        "log_eta": float(best.x[0]),
        "eta": math.exp(float(best.x[0])),
        "z": math.exp(-math.exp(float(best.x[0]))),
        "log_fugacities_x1_through_x4": [float(value) for value in best.x[1:]],
        "objective_replays": [float(objective(best.x, distance))] * 2,
    }


def audit_transfer_power() -> int:
    comparisons = 0
    for exponent in range(1, 17):
        matrix = np.zeros((5, 5), dtype=np.float64)
        for packet in range(5):
            matrix += (packet + 1) / 7.0 * TRANSITIONS[packet]
        matrix *= np.power(0.83, np.arange(5, dtype=np.float64))[None, :]
        direct = np.zeros(5, dtype=np.float64)
        direct[0] = 1.0
        for _ in range(exponent):
            direct = direct @ matrix
        observed = scaled_log_matrix_power_mass(matrix, exponent)
        expected = math.log(float(np.sum(direct)))
        if not math.isclose(observed, expected, rel_tol=1e-12, abs_tol=1e-12):
            raise RuntimeError("scaled matrix power failed a direct audit")
        comparisons += 1
    return comparisons


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--distances", type=int, nargs="+", default=[76_000, 77_000, 188_743]
    )
    parser.add_argument("--global-search", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if any(not 0 <= distance <= 4 * PACKET_POSITIONS for distance in args.distances):
        raise ValueError("distance is outside the binary output range")
    if len(LOCAL_TYPES) != 136:
        raise RuntimeError("the weight-22 local type count changed")
    probability_mass = sum(math.exp(row[1]) for row in TYPE_ROWS)
    if not math.isclose(probability_mass, 1.0, rel_tol=1e-12, abs_tol=1e-12):
        raise RuntimeError("the local type probabilities do not sum to one")

    rows = [
        optimize_distance(distance, args.global_search) for distance in args.distances
    ]
    payload = {
        "schema": "riffle-parallelacc-packet-type-bound-probe-v1",
        "evidence_label": "NUMERICAL_COEFFICIENT_UPPER_BOUND",
        "parameters": {
            "packet_width": PACKET_WIDTH,
            "packet_positions": PACKET_POSITIONS,
            "binary_output_length": 4 * PACKET_POSITIONS,
            "outer_shell": "zero-parity weight-22 occupation-one",
            "outer_shell_words": MINIMUM_SHELL_WORDS,
            "local_packet_types": len(LOCAL_TYPES),
        },
        "method": (
            "For positive x and 0<z<1, coefficient domination bounds "
            "A(h,<=D) by z^(-D)e0^T M(x,z)^N 1 / x^h. The result is "
            "divided by the exact histogram sequence count T_N(h), averaged "
            "over all 136 shell types, and optimized numerically."
        ),
        "distance_rows": rows,
        "validation": {
            "scaled_matrix_power_direct_comparisons": audit_transfer_power(),
            "local_type_probability_mass": probability_mass,
            "all_checks": "PASS",
        },
        "scope": (
            "This probe tests the loss from multivariate coefficient domination "
            "on one complete outer shell. It is not a complete-code distance "
            "bound and the floating-point optimizer is not a proof certificate."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output": str(args.output), "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
