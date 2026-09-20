#!/usr/bin/env python3
"""Exact zero-state activation table for the LDPCSplitState ensemble.

Write an epoch input as X=(x_m,x_p), where x_m contains the mixed
coordinates and x_p contains one systematic parity block. The compressor is

    B(X) = H x_m + x_p,

where every column of H is sampled uniformly from fixed-weight vectors in
the syndrome space. For m=wt(x_m) and v=wt(x_p), this script computes the exact
with-replacement probability that B(X)=0.

The script also records a valid upper bound after conditioning all mixed
columns of H to be distinct. It does not assume independent activation
events across epochs.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/zero_state_activation_table.json"
)


def xor_fixed_weight_transition(bits: int, column_weight: int) -> np.ndarray:
    transition = np.zeros((bits + 1, bits + 1), dtype=np.float64)
    denominator = math.comb(bits, column_weight)
    for current_weight in range(bits + 1):
        minimum_intersection = max(0, column_weight - (bits - current_weight))
        maximum_intersection = min(column_weight, current_weight)
        for intersection in range(minimum_intersection, maximum_intersection + 1):
            following_weight = current_weight + column_weight - 2 * intersection
            transition[current_weight, following_weight] += (
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


def distinct_probability(population: int, draws: int) -> float:
    probability = 1.0
    for used in range(draws):
        probability *= (population - used) / population
    return probability


def evaluate(
    *,
    syndrome_bits: int,
    mixed_coordinates: int,
    parity_coordinates: int,
    column_weight: int,
) -> dict[str, object]:
    if parity_coordinates != syndrome_bits:
        raise ValueError("the systematic parity block must match the syndrome width")
    total_coordinates = mixed_coordinates + parity_coordinates
    population = math.comb(syndrome_bits, column_weight)

    transition = xor_fixed_weight_transition(syndrome_bits, column_weight)
    syndrome_weights = iterated_distributions(transition, mixed_coordinates)
    selected_distinct = [
        distinct_probability(population, draws)
        for draws in range(mixed_coordinates + 1)
    ]

    exact = np.zeros((mixed_coordinates + 1, parity_coordinates + 1))
    conditioned_upper = np.zeros_like(exact)
    for mixed_weight in range(mixed_coordinates + 1):
        for parity_weight in range(parity_coordinates + 1):
            probability = syndrome_weights[mixed_weight, parity_weight] / math.comb(
                parity_coordinates, parity_weight
            )
            exact[mixed_weight, parity_weight] = probability
            conditioned_upper[mixed_weight, parity_weight] = min(
                1.0, probability / selected_distinct[mixed_weight]
            )

    # Conditioning on distinct columns gives exact zeros through weight three.
    for mixed_weight in range(mixed_coordinates + 1):
        for parity_weight in range(parity_coordinates + 1):
            total_weight = mixed_weight + parity_weight
            if total_weight in (1, 2, 3):
                conditioned_upper[mixed_weight, parity_weight] = 0.0

    by_total_weight = []
    for total_weight in range(1, total_coordinates + 1):
        rows = []
        average_exact = 0.0
        average_upper = 0.0
        split_denominator = math.comb(total_coordinates, total_weight)
        for mixed_weight in range(
            max(0, total_weight - parity_coordinates),
            min(mixed_coordinates, total_weight) + 1,
        ):
            parity_weight = total_weight - mixed_weight
            split_count = (
                math.comb(mixed_coordinates, mixed_weight)
                * math.comb(parity_coordinates, parity_weight)
            )
            split_probability = split_count / split_denominator
            p_exact = exact[mixed_weight, parity_weight]
            p_upper = conditioned_upper[mixed_weight, parity_weight]
            average_exact += split_probability * p_exact
            average_upper += split_probability * p_upper
            rows.append((mixed_weight, parity_weight, p_exact, p_upper))
        maximum_exact = max(rows, key=lambda row: row[2])
        maximum_upper = max(rows, key=lambda row: row[3])
        by_total_weight.append(
            {
                "total_weight": total_weight,
                "maximum_with_replacement_probability": maximum_exact[2],
                "maximum_with_replacement_split": {
                    "mixed_weight": maximum_exact[0],
                    "parity_weight": maximum_exact[1],
                },
                "maximum_distinct_conditioned_upper_bound": maximum_upper[3],
                "maximum_distinct_conditioned_split": {
                    "mixed_weight": maximum_upper[0],
                    "parity_weight": maximum_upper[1],
                },
                "uniform_256_support_average_with_replacement": average_exact,
                "uniform_256_support_average_distinct_upper_bound": average_upper,
                "uniform_support_average_with_replacement": average_exact,
                "uniform_support_average_distinct_upper_bound": average_upper,
            }
        )

    selected_rows = []
    for mixed_weight in range(mixed_coordinates + 1):
        for parity_weight in range(parity_coordinates + 1):
            probability = exact[mixed_weight, parity_weight]
            upper = conditioned_upper[mixed_weight, parity_weight]
            if (
                mixed_weight + parity_weight <= 8
                or (mixed_weight, parity_weight)
                in ((16, 16), (32, 32), (64, 0), (64, 64), (128, 64))
            ):
                selected_rows.append(
                    {
                        "mixed_weight": mixed_weight,
                        "parity_weight": parity_weight,
                        "total_weight": mixed_weight + parity_weight,
                        "with_replacement_probability": probability,
                        "distinct_conditioned_upper_bound": upper,
                    }
                )

    return {
        "schema": "riffle-ldpcsplitstate-zero-activation-v1",
        "candidate": (
            "Riffle BitShuffle-SplitState "
            f"t={total_coordinates} s={syndrome_bits}"
        ),
        "compressor": "B(x_m,x_p)=H x_m + x_p",
        "parameters": {
            "syndrome_bits": syndrome_bits,
            "mixed_coordinates": mixed_coordinates,
            "parity_coordinates": parity_coordinates,
            "H_column_weight": column_weight,
            "H_column_population": population,
            "all_mixed_columns_distinct_probability": selected_distinct[-1],
        },
        "structural_claims": {
            "odd_total_weight": (
                "B(X) is nonzero deterministically because every column of B "
                "has odd weight"
            ),
            "distinct_columns": (
                "Conditioned on distinct H columns, ker(B) contains no "
                "nonzero words of weights one, two, or three"
            ),
            "first_canonical_failure": (
                "At split (m,v)=(1,3), nonactivation requires the one H "
                "column to equal the three active identity columns and has "
                f"probability 1/C({syndrome_bits},3)"
            ),
        },
        "selected_split_rows": selected_rows,
        "by_total_weight": by_total_weight,
        "scope": (
            "The split table is exact for independently sampled H columns. "
            "The distinct-column values are upper bounds obtained by "
            "conditioning the selected columns to be distinct. The uniform "
            f"support average assumes a uniform permutation of all {total_coordinates} epoch "
            "positions and is only a reference for the packet-placement law. "
            "No product of activation probabilities across epochs is claimed."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--syndrome-bits", type=int, default=64)
    parser.add_argument("--mixed-coordinates", type=int)
    parser.add_argument("--parity-coordinates", type=int)
    parser.add_argument("--column-weight", type=int, default=3)
    args = parser.parse_args()
    mixed_coordinates = (
        args.mixed_coordinates
        if args.mixed_coordinates is not None
        else 3 * args.syndrome_bits
    )
    parity_coordinates = (
        args.parity_coordinates
        if args.parity_coordinates is not None
        else args.syndrome_bits
    )
    payload = evaluate(
        syndrome_bits=args.syndrome_bits,
        mixed_coordinates=mixed_coordinates,
        parity_coordinates=parity_coordinates,
        column_weight=args.column_weight,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "distinct_probability,"
        f"{payload['parameters']['all_mixed_columns_distinct_probability']:.9f}"
    )
    for total_weight in range(1, min(8, len(payload["by_total_weight"])) + 1):
        row = payload["by_total_weight"][total_weight - 1]
        print(
            f"total_weight,{total_weight},"
            f"max_distinct_upper,{row['maximum_distinct_conditioned_upper_bound']:.9e},"
            "uniform_support_average_upper,"
            f"{row['uniform_256_support_average_distinct_upper_bound']:.9e},"
            "worst_split,"
            f"{row['maximum_distinct_conditioned_split']['mixed_weight']},"
            f"{row['maximum_distinct_conditioned_split']['parity_weight']}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
