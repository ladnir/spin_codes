#!/usr/bin/env python3
"""Verify irreducible-sector energies and the sector-norm variance bound."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from verify_pair_kernel_small import analytic_moments
from verify_subset_covariance_blocks_small import (
    covariance_blocks,
    covariance_kernel,
    shell_coefficient,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "sector_projection_energy_K4_B8_r3_w6_exact.json"


def solve_fraction(matrix: list[list[Fraction]], target: int) -> list[Fraction]:
    size = len(matrix)
    augmented = [
        row.copy() + [Fraction(int(index == target))]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            row for row in range(column, size) if augmented[row][column] != 0
        )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            scale = augmented[row][column]
            if scale:
                augmented[row] = [
                    left - scale * right
                    for left, right in zip(
                        augmented[row], augmented[column], strict=True
                    )
                ]
    return [augmented[row][-1] for row in range(size)]


def level_eigenvalue(
    length: int, level: int, sector: int, intersection: int
) -> int:
    outside = length - 2 * sector
    lifted_level = level - sector
    value = 0
    for same_pairs in range(sector + 1):
        outside_intersection = intersection - same_pairs
        if not 0 <= outside_intersection <= lifted_level:
            continue
        remainder = lifted_level - outside_intersection
        if not 0 <= remainder <= outside - lifted_level:
            continue
        value += (
            (-1) ** (sector - same_pairs)
            * math.comb(sector, same_pairs)
            * math.comb(lifted_level, outside_intersection)
            * math.comb(outside - lifted_level, remainder)
        )
    return value


def sector_energies(
    output_bits: int, shell_weight: int
) -> tuple[list[int], list[Fraction], list[dict[str, object]]]:
    coefficients = {
        subset: shell_coefficient(output_bits, shell_weight, subset)
        for subset in range(1 << output_bits)
    }
    maximum_sector = output_bits // 2
    energies = [Fraction(0) for _ in range(maximum_sector + 1)]
    level_totals = []
    projector_rows = []
    for level in range(output_bits + 1):
        subsets = [
            subset
            for subset in range(1 << output_bits)
            if subset.bit_count() == level
        ]
        minimum_intersection = max(0, 2 * level - output_bits)
        intersections = list(range(minimum_intersection, level + 1))
        sectors = list(range(min(level, output_bits - level) + 1))
        eigenmatrix = [
            [
                Fraction(level_eigenvalue(output_bits, level, sector, intersection))
                for intersection in intersections
            ]
            for sector in sectors
        ]
        correlations = []
        for intersection in intersections:
            correlations.append(
                sum(
                    coefficients[first] * coefficients[second]
                    for first in subsets
                    for second in subsets
                    if (first & second).bit_count() == intersection
                )
            )
        total = sum(coefficients[subset] ** 2 for subset in subsets)
        level_totals.append(total)
        reconstructed = Fraction(0)
        for local_sector, sector in enumerate(sectors):
            projector = solve_fraction(eigenmatrix, local_sector)
            energy = sum(
                coefficient * correlation
                for coefficient, correlation in zip(
                    projector, correlations, strict=True
                )
            )
            if energy < 0:
                raise AssertionError("sector projection has negative energy")
            energies[sector] += energy
            reconstructed += energy
            projector_rows.append(
                {
                    "level": level,
                    "sector": sector,
                    "energy": f"{energy.numerator}/{energy.denominator}",
                }
            )
        if reconstructed != total:
            raise AssertionError("sector energies do not reconstruct level energy")
    return level_totals, energies, projector_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=4)
    parser.add_argument("--output-bits", type=int, default=8)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--shell-weight", type=int, default=6)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output_bits > 10:
        parser.error("the exact subset verifier is limited to output_bits <= 10")

    _returns, kernel = covariance_kernel(
        args.message_bits, args.output_bits, args.right_degree
    )
    blocks = covariance_blocks(args.output_bits, kernel)
    level_totals, energies, projector_rows = sector_energies(
        args.output_bits, args.shell_weight
    )
    parseval = (1 << args.output_bits) * math.comb(
        args.output_bits, args.shell_weight
    )
    if sum(level_totals) != parseval or sum(energies) != parseval:
        raise AssertionError("sector energies violate Parseval")

    sector_bound = 0.0
    block_rows = []
    for block, energy in zip(blocks, energies, strict=True):
        maximum_eigenvalue = float(np.linalg.eigvalsh(block["matrix"])[-1])
        sector_bound += maximum_eigenvalue * float(energy)
        block_rows.append(
            {
                "sector": block["sector"],
                "maximum_covariance_eigenvalue": maximum_eigenvalue,
                "projected_coefficient_energy": (
                    f"{energy.numerator}/{energy.denominator}"
                ),
            }
        )

    mean, factorial = analytic_moments(
        args.message_bits, args.output_bits, args.right_degree
    )
    exact_mean = mean[args.shell_weight]
    exact_variance = (
        exact_mean
        + factorial[args.shell_weight]
        - exact_mean * exact_mean
    )
    scale = 2.0 ** (2 * (args.message_bits - args.output_bits))
    ratio_bound = scale * sector_bound / float(exact_mean)
    exact_ratio = float(exact_variance / exact_mean)

    payload = {
        "schema": "pure-ea-sector-projection-energy-small-v1",
        "status": "EXACT_PROJECTOR_ENERGIES_WITH_BINARY64_EIGENVALUE_BOUND",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "shell_weight": args.shell_weight,
        },
        "parseval_energy": parseval,
        "exact_variance_to_mean": exact_ratio,
        "sector_norm_variance_to_mean_upper": ratio_bound,
        "loss_over_exact": ratio_bound / exact_ratio,
        "sectors": block_rows,
        "level_sector_energies": projector_rows,
        "verified_identities": [
            "sector energies reconstruct every row-subset level energy exactly",
            "total sector energy equals the exact Parseval energy",
        ],
        "scope": [
            "All projector coefficients and energies use exact rational arithmetic.",
            "The covariance eigenvalue bound uses nondirected binary64 arithmetic.",
            "The result validates only the displayed small instance.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "exact_variance_to_mean": exact_ratio,
                "sector_norm_variance_to_mean_upper": ratio_bound,
                "loss_over_exact": ratio_bound / exact_ratio,
                "parseval_energy": parseval,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
