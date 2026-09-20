#!/usr/bin/env python3
"""Select a fixed low-dimensional RM(2,7) subcode for SplitState."""

from __future__ import annotations

import argparse
import json
import random
from collections import Counter
from pathlib import Path


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/smaller_state/s16_rm2sub_selection.json"
)
DEFAULT_HISTOGRAM = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/smaller_state/s16_rm2sub_a_histogram.csv"
)


MONOMIALS = tuple((left, right) for left in range(7) for right in range(left + 1, 7))


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


def sample_quadratics(rng: random.Random, count: int) -> list[int]:
    while True:
        words: list[int] = []
        while len(words) < count:
            word = rng.randrange(1, 1 << len(MONOMIALS))
            if rank(words + [word]) > len(words):
                words.append(word)
        if rank(words) == count:
            return words


def evaluate_quadratic(mask: int, point: int) -> int:
    value = 0
    for index, (left, right) in enumerate(MONOMIALS):
        if (mask >> index) & 1:
            value ^= ((point >> left) & 1) & ((point >> right) & 1)
    return value


def coordinate_columns(quadratics: list[int]) -> list[int]:
    columns = []
    for point in range(128):
        column = 1 | (point << 1)
        for output, mask in enumerate(quadratics):
            column |= evaluate_quadratic(mask, point) << (8 + output)
        columns.append(column)
    return columns


def weight_four_kernel_count(columns: list[int]) -> int:
    pairs = Counter(
        columns[left] ^ columns[right]
        for left in range(len(columns))
        for right in range(left + 1, len(columns))
    )
    collisions = sum(count * (count - 1) // 2 for count in pairs.values())
    if collisions % 3:
        raise ArithmeticError("pair collisions did not form weight-four words")
    return collisions // 3


def generator_words(columns: list[int], state_bits: int) -> list[int]:
    generators = []
    for information in range(state_bits):
        word = 0
        for coordinate, column in enumerate(columns):
            word |= ((column >> information) & 1) << coordinate
        generators.append(word)
    return generators


def enumerate_spectrum(generators: list[int], state_bits: int) -> list[int]:
    histogram = [0] * 129
    current = 0
    histogram[0] = 1
    for index in range(1, 1 << state_bits):
        current ^= generators[(index & -index).bit_length() - 1]
        histogram[current.bit_count()] += 1
    return histogram


def evaluate(args: argparse.Namespace) -> tuple[dict[str, object], list[int]]:
    trials = []
    quadratic_count = args.state_bits - 8
    for offset in range(args.trials):
        seed = args.seed + offset
        quadratics = sample_quadratics(random.Random(seed), quadratic_count)
        columns = coordinate_columns(quadratics)
        a4 = weight_four_kernel_count(columns)
        trials.append((a4, seed, quadratics, columns))
    finalists = sorted(trials)[: args.finalists]
    audited = []
    for a4, seed, quadratics, columns in finalists:
        generators = generator_words(columns, args.state_bits)
        spectrum = enumerate_spectrum(generators, args.state_bits)
        distance = next(weight for weight in range(1, 129) if spectrum[weight])
        audited.append((a4, -distance, seed, quadratics, columns, generators, spectrum))
    if args.selection_order == "spectrum-first":
        selected = min(
            audited,
            key=lambda row: (
                row[1],
                tuple(row[6][weight] for weight in range(1, 65)),
                row[0],
                row[2],
            ),
        )
    elif args.selection_order == "distance-first":
        selected = min(audited, key=lambda row: (row[1], row[0], row[2]))
    else:
        selected = min(audited)
    a4, negative_distance, seed, quadratics, columns, generators, spectrum = selected
    distance = -negative_distance
    payload = {
        "schema": "riffle-splitstate-rm2sub-selection-v1",
        "candidate": (
            "Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul-RM2Sub "
            f"t=128 s={args.state_bits}"
        ),
        "parameters": {
            "seed_start": args.seed,
            "trials": args.trials,
            "finalists_exactly_enumerated": args.finalists,
            "selection_order": args.selection_order,
            "state_bits": args.state_bits,
            "step_bits": 128,
            "basis": (
                "constant, seven linear forms, and "
                f"{quadratic_count} selected homogeneous quadratic forms"
            ),
        },
        "selected": {
            "seed": seed,
            "quadratic_masks_hex": [f"{mask:06x}" for mask in quadratics],
            "B_columns_hex": [
                f"{column:0{(args.state_bits + 3) // 4}x}" for column in columns
            ],
            "A_generator_words_hex": [f"{word:032x}" for word in generators],
            "weight_four_kernel_words": a4,
            "weight_four_nonactivation_probability": a4 / 10668000,
            "minimum_A_distance": distance,
            "ba_zero": True,
            "coordinate_audit": {
                "zero_coordinate_forms": 0,
                "distinct_coordinate_forms": 128,
                "triple_dependencies": 0,
            },
            "coordinate_profile": {
                "zero_forms": 0,
                "distinct_forms": 128,
                "triple_dependencies": 0,
            },
        },
        "proof": {
            "BA_zero": "Products of two selected basis functions have degree at most four, so their sums over F_2^7 vanish.",
            "minimum_A_distance_lower_bound": "The selected code is a subcode of RM(2,7), whose minimum distance is 32.",
            "coordinate_profile": "The constant and all seven linear forms make columns distinct; their constant coordinate is one, so no three columns sum to zero.",
        },
        "top_trials": [
            {"seed": row[1], "weight_four_kernel_words": row[0]}
            for row in sorted(trials)[:20]
        ],
        "audited_finalists": [
            {
                "seed": row[2],
                "weight_four_kernel_words": row[0],
                "minimum_A_distance": -row[1],
                "weight_32_A_words": row[6][32],
            }
            for row in sorted(audited, key=lambda row: (row[1], row[0], row[2]))
        ],
    }
    return payload, spectrum


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=0xA2160000)
    parser.add_argument("--trials", type=int, default=4096)
    parser.add_argument("--finalists", type=int, default=16)
    parser.add_argument("--state-bits", type=int, default=16)
    parser.add_argument(
        "--selection-order",
        choices=("kernel-first", "distance-first", "spectrum-first"),
        default="kernel-first",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--histogram", type=Path, default=DEFAULT_HISTOGRAM)
    args = parser.parse_args()
    if not 9 <= args.state_bits <= 29:
        parser.error("--state-bits must be between 9 and 29")
    payload, spectrum = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.histogram.write_text(
        "weight,count\n"
        + "".join(
            f"{weight},{count}\n"
            for weight, count in enumerate(spectrum)
            if count
        ),
        encoding="utf-8",
    )
    selected = payload["selected"]
    print(f"seed,{selected['seed']}")
    print(f"minimum_A_distance,{selected['minimum_A_distance']}")
    print(f"weight_four_kernel_words,{selected['weight_four_kernel_words']}")
    print(f"wrote,{args.output}")
    print(f"wrote,{args.histogram}")


if __name__ == "__main__":
    main()
