#!/usr/bin/env python3
"""Expected spectrum gate for a scaled nested SplitState constituent."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/"
    "receipts/nested_ensemble_spectrum.json"
)


def xor_fixed_weight_transition(bits: int, column_weight: int) -> np.ndarray:
    transition = np.zeros((bits + 1, bits + 1), dtype=np.float64)
    denominator = math.comb(bits, column_weight)
    for current_weight in range(bits + 1):
        for intersection in range(
            max(0, column_weight - (bits - current_weight)),
            min(column_weight, current_weight) + 1,
        ):
            following = current_weight + column_weight - 2 * intersection
            transition[current_weight, following] += (
                math.comb(current_weight, intersection)
                * math.comb(bits - current_weight, column_weight - intersection)
                / denominator
            )
    return transition


def iterated_distributions(
    transition: np.ndarray, maximum_steps: int
) -> np.ndarray:
    result = np.zeros((maximum_steps + 1, transition.shape[0]), dtype=np.float64)
    result[0, 0] = 1.0
    for step in range(1, maximum_steps + 1):
        result[step] = result[step - 1] @ transition
    return result


def accumulator_weight_transition(length: int) -> np.ndarray:
    current = np.zeros((length + 1, 2, length + 1), dtype=np.float64)
    current[0, 0, 0] = 1.0
    maximum_used = 0
    maximum_weight = 0
    for _ in range(length):
        following = np.zeros_like(current)
        used = slice(0, maximum_used + 1)
        weights = slice(0, maximum_weight + 1)
        following[used, 0, weights] += current[used, 0, weights]
        following[0 : maximum_used + 1, 1, 1 : maximum_weight + 2] += current[
            used, 1, weights
        ]
        following[1 : maximum_used + 2, 1, 1 : maximum_weight + 2] += current[
            used, 0, weights
        ]
        following[1 : maximum_used + 2, 0, 0 : maximum_weight + 1] += current[
            used, 1, weights
        ]
        current = following
        maximum_used += 1
        maximum_weight += 1
    counts = np.sum(current, axis=1)
    for input_weight in range(length + 1):
        counts[input_weight] /= math.comb(length, input_weight)
    return counts


def distinct_probability(population: int, draws: int) -> float:
    probability = 1.0
    for used in range(draws):
        probability *= (population - used) / population
    return probability


def expected_spectrum(
    *,
    state_bits: int,
    p_column_weight: int,
    parity_column_weight: int,
    accumulator_depth: int,
) -> np.ndarray:
    auxiliary_bits = 2 * state_bits
    length = 4 * state_bits
    injection = iterated_distributions(
        xor_fixed_weight_transition(auxiliary_bits, p_column_weight), state_bits
    )
    accumulator = accumulator_weight_transition(auxiliary_bits)
    auxiliary = injection
    for _ in range(accumulator_depth):
        auxiliary = auxiliary @ accumulator

    mixed_coordinates = state_bits + auxiliary_bits
    parity = iterated_distributions(
        xor_fixed_weight_transition(state_bits, parity_column_weight),
        mixed_coordinates,
    )
    spectrum = np.zeros(length + 1, dtype=np.float64)
    for information_weight in range(state_bits + 1):
        information_count = math.comb(state_bits, information_weight)
        for auxiliary_weight, auxiliary_probability in enumerate(
            auxiliary[information_weight]
        ):
            if auxiliary_probability == 0.0:
                continue
            base = information_weight + auxiliary_weight
            distribution = parity[base]
            spectrum[base : base + len(distribution)] += (
                information_count * auxiliary_probability * distribution
            )
    return spectrum


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    spectrum = expected_spectrum(
        state_bits=args.state_bits,
        p_column_weight=args.p_column_weight,
        parity_column_weight=args.parity_column_weight,
        accumulator_depth=args.accumulator_depth,
    )
    mixed_coordinates = 3 * args.state_bits
    population = math.comb(args.state_bits, args.parity_column_weight)
    distinct = distinct_probability(population, mixed_coordinates)
    cutoffs = sorted(set((*args.cutoffs, args.target_distance - 1)))
    rows = []
    for cutoff in cutoffs:
        expected_bad = float(np.sum(spectrum[1 : cutoff + 1]))
        rows.append(
            {
                "cutoff": cutoff,
                "expected_nonzero_words_at_or_below_cutoff": expected_bad,
                "distinct_B_probability": distinct,
                "joint_good_probability_lower_bound": max(0.0, distinct - expected_bad),
                "existence_gate_passes": expected_bad < distinct,
            }
        )
    mass_error = abs(float(np.sum(spectrum)) - math.ldexp(1.0, args.state_bits))
    return {
        "schema": "riffle-splitstate-scaled-nested-ensemble-v1",
        "candidate": (
            "Riffle BCHPerm-TransposeBitShuffle-SplitState "
            f"t={4 * args.state_bits} s={args.state_bits}"
        ),
        "parameters": {
            "state_bits": args.state_bits,
            "length": 4 * args.state_bits,
            "dimension": args.state_bits,
            "auxiliary_bits": 2 * args.state_bits,
            "mixed_coordinates": mixed_coordinates,
            "p_column_weight": args.p_column_weight,
            "parity_column_weight": args.parity_column_weight,
            "accumulator_depth": args.accumulator_depth,
            "target_distance": args.target_distance,
            "distinct_B_probability": distinct,
        },
        "checks": {
            "spectrum_mass_absolute_error": mass_error,
            "expected_zero_weight_words": float(spectrum[0]),
        },
        "distance_gates": rows,
        "expected_weight_spectrum": {
            str(weight): float(value)
            for weight, value in enumerate(spectrum)
            if value != 0.0
        },
        "scope": (
            "Exact expected spectrum for the with-replacement P and B ensemble "
            "and independent uniform interleavers. The joint lower bound uses "
            "a union bound between distinct B columns and absence of low-weight "
            "A codewords. It proves ensemble existence, not a fixed instance."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-bits", type=int, default=32)
    parser.add_argument("--p-column-weight", type=int, default=7)
    parser.add_argument("--parity-column-weight", type=int, default=3)
    parser.add_argument("--accumulator-depth", type=int, default=2)
    parser.add_argument("--target-distance", type=int, default=20)
    parser.add_argument("--cutoffs", type=int, nargs="+", default=(15, 19, 23, 27))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"distinct_B_probability,{payload['parameters']['distinct_B_probability']:.12g}"
    )
    for row in payload["distance_gates"]:
        print(
            f"cutoff,{row['cutoff']},expected_bad,"
            f"{row['expected_nonzero_words_at_or_below_cutoff']:.12g},"
            f"joint_lower,{row['joint_good_probability_lower_bound']:.12g},"
            f"passes,{int(row['existence_gate_passes'])}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
