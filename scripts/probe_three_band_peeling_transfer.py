#!/usr/bin/env python3
"""Boundary-budget probe for a peeled fixed-band collision cluster.

Suppose a current tile contains ``r`` active blocks.  Pair capacity one makes
their tile labels in each other band distinct.  Conditional on the union
weight already present in such a boundary tile, the independent local
coordinate permutation gives an exact hypergeometric moment for the number of
new groups introduced by the peeled block.

Generalized Holder splits the one coupled ``r``-way OR kernel in the peeled
band into one rooted diagonal factor per block.  The resulting one-block
factor is bounded by the existing BCH cell-cap transport relaxation.  Across
the ``r`` blocks, the sum of all ``2r`` boundary union weights is at most the
global support cutoff 106.

For speed, boundary weights are rounded *up* to a small grid.  Incremental
moments increase with the existing union weight, so this is safe.  A rounded
grid value is charged the smallest integer weight that rounds to it; hence the
DP keeps the original budget 106 without excluding any real boundary vector.
This script is a floating diagnostic; successful rows should later use exact
moments and an outward factor/DP certificate.
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from probe_bch_three_band_cell_cap_motifs import (
    BAND_SIZES,
    cell_caps,
    factor_table,
    greedy_factor_bound,
)
from probe_bch_three_band_transport_lp import (
    exact_rows,
    remove_trivial_codewords,
    triples,
)


DEFAULT_SUPPORT_CUTOFF = 106


def incremental_table(columns: int, existing: int, pole: float) -> np.ndarray:
    """Moment of newly occupied coordinates for a fixed existing union."""

    result = []
    for weight in range(columns + 1):
        denominator = math.comb(columns, weight)
        value = sum(
            math.comb(existing, intersection)
            * math.comb(columns - existing, weight - intersection)
            / denominator
            * pole ** (weight - intersection)
            for intersection in range(
                max(0, existing + weight - columns), min(existing, weight) + 1
            )
        )
        result.append(value)
    return np.asarray(result)


def upward_grid(maximum: int, step: int) -> list[int]:
    values = list(range(0, maximum + 1, step))
    if values[-1] != maximum:
        values.append(maximum)
    return values


def preimage_lower_bounds(grid: list[int]) -> dict[int, int]:
    """Smallest nonnegative integer rounded upward to each grid value."""

    result = {0: 0}
    for previous, current in zip(grid, grid[1:]):
        result[current] = previous + 1
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r-min", type=int, default=3)
    parser.add_argument("--r-max", type=int, default=7)
    parser.add_argument("--grid-step", type=int, default=4)
    parser.add_argument("--tilt", type=float, default=0.05)
    parser.add_argument("--support-cutoff", type=int, default=DEFAULT_SUPPORT_CUTOFF)
    args = parser.parse_args()
    if not (
        1 <= args.r_min <= args.r_max <= 64
        and 1 <= args.grid_step <= 8
        and 0 < args.tilt < 1
        and 0 <= args.support_cutoff <= 128 * 64
    ):
        raise SystemExit("peeling transfer: invalid arguments")

    states = triples()
    rows = remove_trivial_codewords(states, exact_rows(states))
    caps, diagonal_mass, by_total = cell_caps(states, rows)
    grids = [upward_grid(size, args.grid_step) for size in BAND_SIZES]
    lower_bounds = [preimage_lower_bounds(grid) for grid in grids]
    increments = {
        (band, existing): incremental_table(BAND_SIZES[band], existing, args.tilt)
        for band in range(3)
        for existing in grids[band]
    }

    print("three-band peeled-cluster boundary transfer probe")
    print(
        f"sizes={args.r_min}..{args.r_max} tilt={args.tilt:.12f} "
        f"grid_step={args.grid_step} support_cutoff={args.support_cutoff}"
    )
    for cluster_size in range(args.r_min, args.r_max + 1):
        rooted = {
            band: factor_table(
                BAND_SIZES[band], cluster_size, args.tilt
            )
            for band in range(3)
        }
        # Each option is one peeled block with its two rounded boundary
        # weights.  Take the worst choice of peeled band for that cost.
        best_by_cost: dict[int, float] = {}
        witness_by_cost: dict[int, tuple[int, int, int]] = {}
        for low_band in range(3):
            external = [band for band in range(3) if band != low_band]
            for first_weight in grids[external[0]]:
                for second_weight in grids[external[1]]:
                    cost = (
                        lower_bounds[external[0]][first_weight]
                        + lower_bounds[external[1]][second_weight]
                    )
                    if cost > args.support_cutoff:
                        continue
                    factors = tuple(
                        rooted[band]
                        if band == low_band
                        else increments[
                            (
                                band,
                                first_weight
                                if band == external[0]
                                else second_weight,
                            )
                        ]
                        for band in range(3)
                    )
                    value = math.log2(
                        greedy_factor_bound(
                            states,
                            caps,
                            diagonal_mass,
                            by_total,
                            factors,
                        )
                    )
                    if value > best_by_cost.get(cost, -math.inf):
                        best_by_cost[cost] = value
                        witness_by_cost[cost] = (
                            low_band,
                            first_weight,
                            second_weight,
                        )

        dp = [-math.inf] * (args.support_cutoff + 1)
        paths: list[list[tuple[int, int, int]] | None] = [None] * (
            args.support_cutoff + 1
        )
        dp[0] = 0.0
        paths[0] = []
        for _block in range(cluster_size):
            following = [-math.inf] * (args.support_cutoff + 1)
            following_paths: list[list[tuple[int, int, int]] | None] = [None] * (
                args.support_cutoff + 1
            )
            for used, prefix in enumerate(dp):
                if prefix == -math.inf:
                    continue
                for cost, value in best_by_cost.items():
                    if used + cost > args.support_cutoff:
                        continue
                    candidate = prefix + value
                    if candidate > following[used + cost]:
                        following[used + cost] = candidate
                        following_paths[used + cost] = (
                            paths[used] + [witness_by_cost[cost]]  # type: ignore[operator]
                        )
            dp = following
            paths = following_paths
        best = max(dp)
        used = max(range(len(dp)), key=dp.__getitem__)
        print(
            f"cluster_size={cluster_size} lower_preimage_budget={args.support_cutoff} "
            f"transfer_log2={best:.12f} used_budget={used} "
            f"witness={paths[used]}"
        )
    print("status=SAFE_MONOTONE_GRID_RELAXATION_FLOATING_FACTORS")


if __name__ == "__main__":
    main()
