#!/usr/bin/env python3
"""Generate an authenticated RM(2,m)-subcode pair for RM2Sub calibration.

For t=2^m, choose s independent evaluation vectors from the constant,
linear, and homogeneous-quadratic monomials on F_2^m.  The resulting map A
has length t and dimension s.  Its generator matrix also defines B, so BA=0
for m>=5.  The script exhaustively enumerates A and derives ker(B)=A^perp by
an exact integer MacWilliams transform.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from collections import Counter
from pathlib import Path


def rank(words: list[int]) -> int:
    basis: dict[int, int] = {}
    for original in words:
        word = original
        while word:
            pivot = word.bit_length() - 1
            if pivot in basis:
                word ^= basis[pivot]
            else:
                basis[pivot] = word
                break
    return len(basis)


def sample_independent_masks(
    rng: random.Random, universe_bits: int, count: int
) -> list[int]:
    words: list[int] = []
    while len(words) < count:
        word = rng.randrange(1, 1 << universe_bits)
        if rank(words + [word]) > len(words):
            words.append(word)
    return words


def evaluate_quadratic(mask: int, point: int, monomials: list[tuple[int, int]]) -> int:
    value = 0
    for index, (left, right) in enumerate(monomials):
        if (mask >> index) & 1:
            value ^= ((point >> left) & 1) & ((point >> right) & 1)
    return value


def coordinate_columns(
    variables: int, quadratics: list[int], monomials: list[tuple[int, int]]
) -> list[int]:
    columns = []
    for point in range(1 << variables):
        column = 1 | (point << 1)
        for output, mask in enumerate(quadratics):
            column |= evaluate_quadratic(mask, point, monomials) << (
                1 + variables + output
            )
        columns.append(column)
    return columns


def kernel_weight_four_count(columns: list[int]) -> int:
    pairs = Counter(
        columns[left] ^ columns[right]
        for left in range(len(columns))
        for right in range(left + 1, len(columns))
    )
    collisions = sum(count * (count - 1) // 2 for count in pairs.values())
    if collisions % 3:
        raise ArithmeticError("pair collisions do not form weight-four words")
    return collisions // 3


def generator_words(columns: list[int], state_bits: int) -> list[int]:
    return [
        sum(((column >> bit) & 1) << coordinate for coordinate, column in enumerate(columns))
        for bit in range(state_bits)
    ]


def enumerate_spectrum(generators: list[int], state_bits: int, length: int) -> list[int]:
    histogram = [0] * (length + 1)
    word = 0
    histogram[0] = 1
    for index in range(1, 1 << state_bits):
        word ^= generators[(index & -index).bit_length() - 1]
        histogram[word.bit_count()] += 1
    return histogram


def krawtchouk(length: int, degree: int, point: int) -> int:
    return sum(
        (-1 if intersection & 1 else 1)
        * math.comb(point, intersection)
        * math.comb(length - point, degree - intersection)
        for intersection in range(
            max(0, degree - (length - point)), min(degree, point) + 1
        )
    )


def macwilliams_kernel(dual: list[int], redundancy: int) -> list[int]:
    length = len(dual) - 1
    denominator = 1 << redundancy
    support = [(weight, count) for weight, count in enumerate(dual) if count]
    kernel = []
    for weight in range(length + 1):
        numerator = sum(
            count * krawtchouk(length, weight, dual_weight)
            for dual_weight, count in support
        )
        if numerator % denominator:
            raise ArithmeticError("nonintegral MacWilliams coefficient")
        kernel.append(numerator // denominator)
    if any(count < 0 for count in kernel):
        raise ArithmeticError("negative MacWilliams coefficient")
    if sum(kernel) != 1 << (length - redundancy):
        raise ArithmeticError("kernel mass mismatch")
    return kernel


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step-bits", type=int, required=True)
    parser.add_argument("--state-bits", type=int, required=True)
    parser.add_argument("--trials", type=int, default=1024)
    parser.add_argument("--finalists", type=int, default=16)
    parser.add_argument("--seed", type=int, default=0xCA110000)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.step_bits < 32 or args.step_bits & (args.step_bits - 1):
        parser.error("step size must be a power of two at least 32")
    variables = args.step_bits.bit_length() - 1
    monomials = [
        (left, right)
        for left in range(variables)
        for right in range(left + 1, variables)
    ]
    quadratic_count = args.state_bits - 1 - variables
    if variables < 5:
        parser.error("the degree-four orthogonality proof requires at least five variables")
    if not 0 <= quadratic_count <= len(monomials):
        parser.error("state dimension is outside the available RM(2,m) basis")

    trials = []
    for offset in range(args.trials):
        seed = args.seed + offset
        quadratics = sample_independent_masks(
            random.Random(seed), len(monomials), quadratic_count
        )
        columns = coordinate_columns(variables, quadratics, monomials)
        a4 = kernel_weight_four_count(columns)
        trials.append((a4, seed, quadratics, columns))
    audited = []
    for a4, seed, quadratics, columns in sorted(trials)[: args.finalists]:
        generators = generator_words(columns, args.state_bits)
        spectrum = enumerate_spectrum(generators, args.state_bits, args.step_bits)
        distance = next(weight for weight, count in enumerate(spectrum[1:], 1) if count)
        audited.append((a4, -distance, seed, quadratics, columns, generators, spectrum))
    selected = min(
        audited,
        key=lambda row: (
            row[0],
            row[1],
            tuple(row[6][weight] for weight in range(1, args.step_bits // 2 + 1)),
            row[2],
        ),
    )
    a4, negative_distance, seed, quadratics, columns, generators, spectrum = selected
    distance = -negative_distance
    if rank(generators) != args.state_bits:
        raise ArithmeticError("selected A generator is not full rank")
    if any((left & right).bit_count() & 1 for left in generators for right in generators):
        raise ArithmeticError("BA=0 audit failed")
    kernel = macwilliams_kernel(spectrum, args.state_bits)
    kernel_distance = next(weight for weight, count in enumerate(kernel[1:], 1) if count)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"t{args.step_bits}_s{args.state_bits}"
    selection_path = args.output_dir / f"{stem}_selection.json"
    a_path = args.output_dir / f"{stem}_a_spectrum.json"
    b_path = args.output_dir / f"{stem}_b_kernel_spectrum.json"
    selection = {
        "schema": "rm2sub-calibration-constituent-v1",
        "status": "EXACTLY_AUDITED_SELECTED_SAMPLE",
        "parameters": {
            "step_bits": args.step_bits,
            "state_bits": args.state_bits,
            "variables": variables,
            "trials": args.trials,
            "finalists_exhaustively_enumerated": args.finalists,
            "seed_start": args.seed,
        },
        "selected": {
            "seed": seed,
            "quadratic_masks_hex": [f"{mask:x}" for mask in quadratics],
            "B_columns_hex": [
                f"{column:0{(args.state_bits + 3) // 4}x}" for column in columns
            ],
            "A_generator_words_hex": [
                f"{word:0{args.step_bits // 4}x}" for word in generators
            ],
            "minimum_A_distance": distance,
            "minimum_kernel_distance": kernel_distance,
            "weight_four_kernel_words": a4,
            "BA_zero": True,
            "distinct_nonzero_columns": len(set(columns)) == args.step_bits and 0 not in columns,
        },
        "limitations": [
            "The search selects one deterministic calibration constituent; it does not prove optimality.",
            "The exact audits apply to the selected constituent after the finite random search.",
        ],
    }
    a_payload = {
        "schema": "rm2sub-calibration-a-spectrum-v1",
        "candidate": stem,
        "checks": {
            "spectrum_mass": sum(spectrum),
            "minimum_distance": distance,
            "BA_zero": True,
        },
        "spectrum": [
            {"weight": weight, "count": count}
            for weight, count in enumerate(spectrum)
            if count
        ],
    }
    b_payload = {
        "schema": "riffle-splitstate-fixed-b-kernel-spectrum-v1",
        "candidate": stem,
        "parameters": {
            "length": args.step_bits,
            "redundancy": args.state_bits,
            "kernel_dimension": args.step_bits - args.state_bits,
        },
        "checks": {
            "dual_mass": sum(spectrum),
            "kernel_mass": sum(kernel),
            "minimum_kernel_distance": kernel_distance,
        },
        "by_total_weight": [
            {
                "total_weight": weight,
                "kernel_words": count,
                "shell_size": math.comb(args.step_bits, weight),
                "fixed_nonactivation_probability": count / math.comb(args.step_bits, weight),
                "uniform_support_average_distinct_upper_bound": count / math.comb(args.step_bits, weight),
            }
            for weight, count in enumerate(kernel)
        ],
        "scope": "Exact integer MacWilliams transform of the selected A spectrum.",
    }
    for path, payload in (
        (selection_path, selection),
        (a_path, a_payload),
        (b_path, b_payload),
    ):
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"wrote,{path}")
    print(f"minimum_A_distance,{distance}")
    print(f"minimum_kernel_distance,{kernel_distance}")
    print(f"weight_four_kernel_words,{a4}")


if __name__ == "__main__":
    main()
