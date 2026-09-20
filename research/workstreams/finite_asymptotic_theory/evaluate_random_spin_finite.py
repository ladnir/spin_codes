#!/usr/bin/env python3
"""Evaluate the exact finite Random SPIN first-moment formula.

The mathematical model is the proof-model ensemble in
RANDOM_SPIN_PROOF_AUDIT.md. The main calculation uses binary64 arithmetic.
It is a diagnostic, not an outward-rounded certificate.

The small-instance validation is exact. It compares an integer transfer
recurrence with exhaustive enumeration of all inputs and convolution setups.
All requested configurations run sequentially.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from dataclasses import dataclass

import numpy as np


LOG2 = math.log(2.0)


@dataclass(frozen=True)
class Parameters:
    length: int
    block: int
    memory: int
    distance: int

    def validate(self) -> None:
        if self.length <= 0:
            raise ValueError("length must be positive")
        if self.block <= 0 or self.block % 2:
            raise ValueError("block must be positive and even")
        if self.length % self.block:
            raise ValueError("block must divide length")
        if self.memory <= 0:
            raise ValueError("memory must be positive")
        if not 0 <= self.distance <= self.length:
            raise ValueError("distance must lie in [0,length]")


def log2_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
    ) / LOG2


def logaddexp2(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return -math.inf
    maximum = float(finite.max())
    return maximum + math.log2(float(np.exp2(finite - maximum).sum()))


def local_outer_weight_probability(block: int) -> np.ndarray:
    """Return the local weight law for a uniform input and random injection."""

    dimension = block // 2
    probability = np.zeros(block + 1, dtype=np.float64)
    probability[0] = math.ldexp(1.0, -dimension)
    nonzero_input = 1.0 - probability[0]
    denominator_log2 = block + math.log1p(-math.ldexp(1.0, -block)) / LOG2
    for weight in range(1, block + 1):
        probability[weight] = math.exp2(
            math.log2(nonzero_input)
            + log2_choose(block, weight)
            - denominator_log2
        )
    probability /= probability.sum()
    return probability


def global_outer_weight_probability(length: int, block: int) -> np.ndarray:
    """Return the weight law of a uniform message and random block outer."""

    distribution = np.asarray([1.0], dtype=np.float64)
    local = local_outer_weight_probability(block)
    for _ in range(length // block):
        distribution = np.convolve(distribution, local)
    distribution /= distribution.sum()
    return distribution


def inner_joint_lower_tail(parameters: Parameters) -> np.ndarray:
    """Return Pr[H=h and W<=D] for a fair input and random inner.

    Here H is the input weight and W is the output weight. The input bits are
    independent and fair. The convolution setup follows equation (I) in the
    proof audit.
    """

    n = parameters.length
    m = parameters.memory
    distance = parameters.distance
    current = np.zeros((m + 1, n + 1, distance + 1), dtype=np.float64)
    current[m, 0, 0] = 1.0

    for position in range(n):
        following = np.zeros_like(current)
        h_count = position + 1
        w_count = min(position, distance) + 1
        active = current[:m, :h_count, :w_count]
        active_sum = active.sum(axis=0)

        # A nonzero state makes the output fair. The input is also fair.
        following[1 : m + 1, :h_count, :w_count] += 0.25 * active
        following[1 : m + 1, 1 : h_count + 1, :w_count] += 0.25 * active
        if distance:
            source_w = min(position, distance - 1) + 1
            following[0, :h_count, 1 : source_w + 1] += (
                0.25 * active_sum[:, :source_w]
            )
            following[0, 1 : h_count + 1, 1 : source_w + 1] += (
                0.25 * active_sum[:, :source_w]
            )

        # At the zero state, the output equals the input.
        zero = current[m, :h_count, :w_count]
        following[m, :h_count, :w_count] += 0.5 * zero
        if distance:
            source_w = min(position, distance - 1) + 1
            following[0, 1 : h_count + 1, 1 : source_w + 1] += (
                0.5 * zero[:, :source_w]
            )
        current = following

    return current.sum(axis=(0, 2))


def exact_transfer_counts(length: int, memory: int) -> dict[tuple[int, int], int]:
    """Count input and setup pairs with an integer transfer recurrence."""

    current: dict[tuple[int, int, int], int] = {(memory, 0, 0): 1}
    active_multiplicity = 1 << (memory - 1)
    zero_multiplicity = 1 << memory
    for _ in range(length):
        following: defaultdict[tuple[int, int, int], int] = defaultdict(int)
        for (zero_run, input_weight, output_weight), count in current.items():
            if zero_run < memory:
                for input_bit in (0, 1):
                    following[
                        (
                            min(memory, zero_run + 1),
                            input_weight + input_bit,
                            output_weight,
                        )
                    ] += count * active_multiplicity
                    following[
                        (0, input_weight + input_bit, output_weight + 1)
                    ] += count * active_multiplicity
            else:
                following[
                    (memory, input_weight, output_weight)
                ] += count * zero_multiplicity
                following[
                    (0, input_weight + 1, output_weight + 1)
                ] += count * zero_multiplicity
        current = dict(following)

    totals: defaultdict[tuple[int, int], int] = defaultdict(int)
    for (_zero_run, input_weight, output_weight), count in current.items():
        totals[(input_weight, output_weight)] += count
    return dict(totals)


def exhaustive_counts(length: int, memory: int) -> dict[tuple[int, int], int]:
    """Count input and setup pairs by direct exhaustive simulation."""

    totals: defaultdict[tuple[int, int], int] = defaultdict(int)
    state_mask = (1 << memory) - 1
    setup_count = 1 << (length * memory)
    for setup in range(setup_count):
        rows = [
            (setup >> (position * memory)) & state_mask
            for position in range(length)
        ]
        for input_word in range(1 << length):
            state = 0
            output_weight = 0
            for position, row in enumerate(rows):
                input_bit = (input_word >> position) & 1
                output_bit = input_bit ^ ((row & state).bit_count() & 1)
                output_weight += output_bit
                state = ((state << 1) | output_bit) & state_mask
            totals[(input_word.bit_count(), output_weight)] += 1
    return dict(totals)


def validate_exact_transfer(
    maximum_length: int = 5, memory: int = 2
) -> dict[str, object]:
    comparisons = 0
    for length in range(1, maximum_length + 1):
        transfer = exact_transfer_counts(length, memory)
        exhaustive = exhaustive_counts(length, memory)
        keys = set(transfer) | set(exhaustive)
        comparisons += len(keys)
        for key in keys:
            if transfer.get(key, 0) != exhaustive.get(key, 0):
                raise AssertionError(
                    f"exact transfer mismatch at n={length}, key={key}"
                )
    return {
        "maximum_length": maximum_length,
        "memory": memory,
        "coefficient_comparisons": comparisons,
        "status": "EXACT_INTEGER_MATCH",
    }


def evaluate(parameters: Parameters, top: int) -> dict[str, object]:
    parameters.validate()
    outer = global_outer_weight_probability(parameters.length, parameters.block)
    joint = inner_joint_lower_tail(parameters)
    h_values = np.arange(parameters.length + 1)
    log2_binomial_probability = np.asarray(
        [
            log2_choose(parameters.length, int(weight)) - parameters.length
            for weight in h_values
        ]
    )
    log2_outer_expected = np.full(parameters.length + 1, -math.inf)
    positive_outer = outer > 0.0
    log2_outer_expected[positive_outer] = (
        parameters.length / 2.0 + np.log2(outer[positive_outer])
    )
    log2_joint = np.full(parameters.length + 1, -math.inf)
    positive_joint = joint > 0.0
    log2_joint[positive_joint] = np.log2(joint[positive_joint])
    log2_inner_probability = log2_joint - log2_binomial_probability
    log2_terms = log2_outer_expected + log2_inner_probability
    log2_terms[0] = -math.inf
    total = logaddexp2(log2_terms)

    order = np.argsort(log2_terms)[::-1]
    dominant = []
    for index in order:
        weight = int(index)
        if len(dominant) >= top or not math.isfinite(float(log2_terms[index])):
            break
        dominant.append(
            {
                "input_weight": weight,
                "relative_input_weight": weight / parameters.length,
                "outer_expected_count_log2": float(log2_outer_expected[index]),
                "inner_lower_tail_probability_log2": float(
                    log2_inner_probability[index]
                ),
                "first_moment_term_log2": float(log2_terms[index]),
            }
        )

    finite_inner = log2_inner_probability[
        np.isfinite(log2_inner_probability) & (h_values > 0)
    ]
    return {
        "status": (
            "binary64 diagnostic; exact finite formulas; not outward rounded"
        ),
        "parameters": {
            "length": parameters.length,
            "dimension": parameters.length // 2,
            "block": parameters.block,
            "outer_blocks": parameters.length // parameters.block,
            "memory": parameters.memory,
            "distance_threshold_inclusive": parameters.distance,
            "relative_distance": parameters.distance / parameters.length,
        },
        "log2_expected_bad_nonzero_codewords": total,
        "markov_failure_probability_upper_log2": min(0.0, total),
        "dominant_input_weights": dominant,
        "mass_checks": {
            "outer_weight_probability_sum": float(outer.sum()),
            "inner_joint_lower_tail_sum": float(joint.sum()),
            "maximum_conditional_inner_probability": float(
                np.exp2(finite_inner).max()
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=256)
    parser.add_argument("--block", type=int, default=16)
    parser.add_argument("--memory", type=int, default=16)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--distance", type=int)
    parser.add_argument("--top", type=int, default=12)
    parser.add_argument("--skip-validation", action="store_true")
    args = parser.parse_args()

    if not 0.0 <= args.relative_distance <= 1.0:
        raise ValueError("relative distance must lie in [0,1]")
    distance = (
        args.distance
        if args.distance is not None
        else math.floor(args.relative_distance * args.length)
    )
    parameters = Parameters(args.length, args.block, args.memory, distance)
    payload = {
        "schema": "random-spin-finite-evaluation-v1",
        "validation": (
            {"status": "SKIPPED"}
            if args.skip_validation
            else validate_exact_transfer()
        ),
        "evaluation": evaluate(parameters, args.top),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
