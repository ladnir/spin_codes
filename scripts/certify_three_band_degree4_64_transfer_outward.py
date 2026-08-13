#!/usr/bin/env python3
"""Directed-outward conditional transfer sweep for degrees 4 through 64.

This is the coarse-grid companion to
``certify_three_band_degree3_transfer_outward.py``.  Boundary weights are
rounded upward to the grid 0,4,8,..., with the physical endpoint appended.
Each grid value is charged the least integer weight that rounds to it, so the
DP budget 106 contains every real boundary vector.  Exact moments, verified
upward float conversion, product guards, and dot-sum guards are identical to
the degree-three certificate.

For every degree r the resulting transfer is multiplied by the complete
connected attachment count ``323 * 378 * C(63,r-1)`` and the composable
pointwise puncture adjustment.  The script proves that all rows are no larger
than the correspondingly adjusted degree-three child envelope.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np

from certify_three_band_degree3_transfer_outward import (
    DOT_GUARD,
    POLE,
    PRODUCT_GUARD,
    ROOT_BITS,
    SUPPORT_CUTOFF,
    incremental_exact,
    upward_float,
)
from outward_log2 import log2_fraction
from probe_bch_three_band_cell_cap_motifs import cell_caps
from probe_bch_three_band_transport_lp import (
    BAND_SIZES,
    exact_rows,
    remove_trivial_codewords,
    triples,
)
from three_band_exact_moments import factor_table_upper
from certify_three_band_exact_length import worst_puncture_block_adjustment


GRID_STEP = 4
MAXIMUM_COMPONENT_BLOCKS = 323
PUNCTURE_BLOCK_ADJUSTMENT = worst_puncture_block_adjustment(POLE)
DEGREE_THREE_CHILD = (
    Fraction(11, 16 * (1 << 30))
    * MAXIMUM_COMPONENT_BLOCKS
    * 378
    * math.comb(63, 2)
    * PUNCTURE_BLOCK_ADJUSTMENT**3
)


def upward_grid(maximum: int) -> list[int]:
    values = list(range(0, maximum + 1, GRID_STEP))
    if values[-1] != maximum:
        values.append(maximum)
    return values


def lower_preimage(grid: list[int]) -> dict[int, int]:
    result = {0: 0}
    for previous, current in zip(grid, grid[1:]):
        result[current] = previous + 1
    return result


def main() -> None:
    state_list = triples()
    states = np.asarray(state_list, dtype=np.int16)
    rows = remove_trivial_codewords(state_list, exact_rows(state_list))
    caps, diagonal_mass, by_total = cell_caps(state_list, rows)
    active_by_total = {
        weight: np.asarray(
            [index for index in by_total[weight] if caps[index]], dtype=np.int32
        )
        for weight, mass in diagonal_mass.items()
        if mass
    }
    grids = [upward_grid(size) for size in BAND_SIZES]
    lower = [lower_preimage(grid) for grid in grids]
    incremental = {
        (band, existing): np.asarray(
            [
                upward_float(incremental_exact(BAND_SIZES[band], existing, weight))
                for weight in range(BAND_SIZES[band] + 1)
            ]
        )
        for band in range(3)
        for existing in grids[band]
    }

    def factor_upper(factors: tuple[np.ndarray, np.ndarray, np.ndarray]) -> Fraction:
        scores = factors[0][states[:, 0]] * factors[1][states[:, 1]]
        scores *= factors[2][states[:, 2]]
        scores = np.nextafter(scores * (1.0 + PRODUCT_GUARD), math.inf)
        diagonal_values = []
        for weight, indices in active_by_total.items():
            order = indices[np.argsort(scores[indices])[::-1]]
            remaining = diagonal_mass[weight]
            selected_indices = []
            selected_caps = []
            for index in order:
                take = min(remaining, caps[int(index)])
                if take:
                    selected_indices.append(int(index))
                    selected_caps.append(upward_float(Fraction(take)))
                    remaining -= take
                if not remaining:
                    break
            if remaining:
                raise AssertionError(f"degree sweep misses diagonal {weight}")
            products = (
                np.asarray(selected_caps)
                * scores[np.asarray(selected_indices, dtype=np.int32)]
            )
            value = float(np.sum(products, dtype=np.float64))
            value = math.nextafter(value * (1.0 + DOT_GUARD), math.inf)
            diagonal_values.append(value)
        return Fraction.from_float(
            math.nextafter(math.fsum(diagonal_values), math.inf)
        )

    worst_child = Fraction(0)
    worst_degree = 0
    for degree in range(4, 65):
        rooted = {
            band: np.asarray(
                [
                    upward_float(value)
                    for value in factor_table_upper(
                        BAND_SIZES[band], degree, POLE, root_bits=ROOT_BITS
                    )
                ]
            )
            for band in range(3)
        }
        best_by_cost: dict[int, Fraction] = {}
        for low_band in range(3):
            external = [band for band in range(3) if band != low_band]
            for first_weight in grids[external[0]]:
                for second_weight in grids[external[1]]:
                    cost = (
                        lower[external[0]][first_weight]
                        + lower[external[1]][second_weight]
                    )
                    if cost > SUPPORT_CUTOFF:
                        continue
                    factors = tuple(
                        rooted[band]
                        if band == low_band
                        else incremental[
                            (
                                band,
                                first_weight
                                if band == external[0]
                                else second_weight,
                            )
                        ]
                        for band in range(3)
                    )
                    value = factor_upper(factors)  # type: ignore[arg-type]
                    if value > best_by_cost.get(cost, Fraction(0)):
                        best_by_cost[cost] = value

        dp = [Fraction(0)] * (SUPPORT_CUTOFF + 1)
        dp[0] = Fraction(1)
        for _block in range(degree):
            following = [Fraction(0)] * (SUPPORT_CUTOFF + 1)
            for used, prefix in enumerate(dp):
                if not prefix:
                    continue
                for cost, value in best_by_cost.items():
                    if used + cost <= SUPPORT_CUTOFF:
                        following[used + cost] = max(
                            following[used + cost], prefix * value
                        )
            dp = following
        transfer = max(dp)
        child = (
            transfer
            * MAXIMUM_COMPONENT_BLOCKS
            * 378
            * math.comb(63, degree - 1)
            * PUNCTURE_BLOCK_ADJUSTMENT**degree
        )
        if child > DEGREE_THREE_CHILD:
            raise SystemExit(
                f"degree sweep: degree {degree} exceeds degree-three child"
            )
        if child > worst_child:
            worst_child = child
            worst_degree = degree
        print(
            f"degree={degree} transfer_log2_upper="
            f"{log2_fraction(transfer).hi} child_log2_upper="
            f"{log2_fraction(child).hi}"
        )

    print("outward three-band degree-4..64 transfer certificate")
    print(f"grid_step={GRID_STEP} worst_degree={worst_degree}")
    print(f"worst_child_log2_upper={log2_fraction(worst_child).hi}")
    print(f"degree_three_child_log2={log2_fraction(DEGREE_THREE_CHILD).hi}")
    print("every_degree_4_64_child_le_degree_three_child=PASS")
    print("status=EXACT_MOMENTS_DIRECTED_OUTWARD_BINARY64")


if __name__ == "__main__":
    main()
