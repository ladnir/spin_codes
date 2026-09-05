#!/usr/bin/env python3
"""Select a fixed sparse 32-by-128 SplitState compressor."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from collections import Counter
from pathlib import Path


DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_t128_s32/"
    "receipts/fixed_b32_selection.json"
)


def weight_three_population() -> list[int]:
    return [
        (1 << left) | (1 << middle) | (1 << right)
        for left, middle, right in itertools.combinations(range(32), 3)
    ]


def weight_four_kernel_count(columns: list[int]) -> int:
    pair_syndromes = Counter(
        columns[left] ^ columns[right]
        for left in range(len(columns))
        for right in range(left + 1, len(columns))
    )
    collision_count = sum(
        count * (count - 1) // 2 for count in pair_syndromes.values()
    )
    if collision_count % 3:
        raise ArithmeticError("pair collisions did not form weight-four words")
    return collision_count // 3


def row_weights(columns: list[int]) -> list[int]:
    return [sum((column >> row) & 1 for column in columns) for row in range(32)]


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    population = weight_three_population()
    identity = [1 << row for row in range(32)]
    best = None
    trial_rows = []
    for trial in range(args.trials):
        seed = args.seed + trial
        rng = random.Random(seed)
        mixed = rng.sample(population, 96)
        columns = mixed + identity
        a4 = weight_four_kernel_count(columns)
        weights = row_weights(columns)
        score = (a4, max(weights) - min(weights), max(weights), seed)
        row = {
            "seed": seed,
            "weight_four_kernel_words": a4,
            "minimum_row_weight": min(weights),
            "maximum_row_weight": max(weights),
            "row_weight_spread": max(weights) - min(weights),
        }
        trial_rows.append(row)
        if best is None or score < best[0]:
            best = (score, row, mixed, columns, weights)
    assert best is not None
    _, best_row, mixed, columns, weights = best
    return {
        "schema": "riffle-splitstate-fixed-b32-selection-v1",
        "candidate": "Riffle BCHPerm-TransposeBitShuffle-SplitState-PreAddMul t=128 s=32",
        "parameters": {
            "seed_start": args.seed,
            "trials": args.trials,
            "syndrome_bits": 32,
            "mixed_columns": 96,
            "mixed_column_weight": 3,
            "parity_columns": 32,
        },
        "selected": {
            **best_row,
            "rank": 32,
            "minimum_kernel_distance_lower_bound": 4,
            "mixed_columns_hex": [f"{column:08x}" for column in mixed],
            "columns_hex": [f"{column:08x}" for column in columns],
            "row_weights": weights,
            "weight_four_nonactivation_probability": (
                best_row["weight_four_kernel_words"] / math.comb(128, 4)
            ),
        },
        "top_trials": sorted(
            trial_rows,
            key=lambda row: (
                row["weight_four_kernel_words"],
                row["row_weight_spread"],
                row["maximum_row_weight"],
                row["seed"],
            ),
        )[:20],
        "scope": (
            "Deterministic seed search. Distinct odd-weight columns prove that "
            "kernel weights one, two, and three are absent. The complete "
            "kernel spectrum requires the separate dual enumeration."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=0xB32A0000)
    parser.add_argument("--trials", type=int, default=4096)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    selected = payload["selected"]
    print(f"seed,{selected['seed']}")
    print(f"weight_four_kernel_words,{selected['weight_four_kernel_words']}")
    print(
        "weight_four_nonactivation_probability,"
        f"{selected['weight_four_nonactivation_probability']:.12g}"
    )
    print(
        f"row_weight_range,{selected['minimum_row_weight']},"
        f"{selected['maximum_row_weight']}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()

