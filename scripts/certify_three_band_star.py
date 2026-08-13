#!/usr/bin/env python3
"""Exact rational certificate for every fixed-band one-tile star.

For a star of r data blocks sharing one tile, the other two bands are
collision-free by pair capacity one.  The committed exact projection rows and
cell caps give an unnormalized inside-weight measure ``mu[a]`` whose outside
coordinates are already charged at a rational Chernoff pole q.

For r independently permuted inside supports, the OR kernel has the exact
positive expansion

    sum_t C(n,t) q^(n-t) (1-q)^t
        (sum_a mu[a] C(n-t,a)/C(n,a))^r.

This avoids both floating arithmetic and the recursive union DP.  The final
row includes the support-106 pole conversion, the exact-length graph factor,
and a composable pointwise puncture adjustment.  All PASS comparisons are
exact rational inequalities; printed logarithms are diagnostics only.
"""

from __future__ import annotations

import math
from collections import Counter
from fractions import Fraction

from certify_three_band_exact_length import (
    graph_factor_upper,
    worst_puncture_block_adjustment,
)
from outward_log2 import log2_fraction
from probe_fixed_three_band_star import (
    BAND_SIZES,
    BLOCKS_PER_TILE,
    TILES,
    load_tables,
)
from punctured_ebch_outer import full_spectrum


TARGET_POLE = Fraction(181, 1000)
SUPPORT_CUTOFF = 106


def exact_capped_inside_measure(
    band: int,
    inside_size: int,
    pole: Fraction,
    projection: dict[int, dict[int, int]],
    outside_slices: Counter,
    inside_boundary: Counter,
    spectrum: list[int],
) -> tuple[Fraction, ...]:
    outside_size = 128 - inside_size
    result = [Fraction(0)] * (inside_size + 1)
    for inside_weight in range(inside_size + 1):
        column_mass = (
            math.comb(inside_size, inside_weight) * (1 << (64 - inside_size))
            - (inside_weight == 0)
            - (inside_weight == inside_size)
        )
        if inside_weight in (0, 1, inside_size - 1, inside_size):
            recovered = 0
            for (row_band, row_inside, outside_weight), raw_count in (
                inside_boundary.items()
            ):
                if row_band != band or row_inside != inside_weight:
                    continue
                count = raw_count
                if (inside_weight, outside_weight) in (
                    (0, 0),
                    (inside_size, outside_size),
                ):
                    count -= 1
                recovered += count
                result[inside_weight] += count * pole**outside_weight
            if recovered != column_mass:
                raise SystemExit("star: incomplete exact inside boundary row")
            continue

        recovered = 0
        for outside_weight in range(1, 9):
            count = outside_slices.get(
                (band, inside_weight, outside_weight), 0
            )
            recovered += count
            result[inside_weight] += count * pole**outside_weight
        remaining = column_mass - recovered
        for outside_weight in range(9, outside_size + 1):
            if not remaining:
                break
            row_cap = projection[band].get(outside_weight, 0)
            if outside_weight == outside_size:
                row_cap -= 1
            total_weight = inside_weight + outside_weight
            diagonal_cap = spectrum[total_weight]
            if total_weight == 128:
                diagonal_cap -= 1
            take = min(remaining, row_cap, diagonal_cap)
            result[inside_weight] += take * pole**outside_weight
            remaining -= take
        if remaining:
            raise SystemExit("star: transportation caps miss a column")
    return tuple(result)


def exact_star_moment(
    columns: int, blocks: int, pole: Fraction, measure: tuple[Fraction, ...]
) -> Fraction:
    total = Fraction(0)
    for avoided in range(columns + 1):
        one_block = Fraction(0)
        available = columns - avoided
        for weight, mass in enumerate(measure):
            if mass and weight <= available:
                one_block += mass * Fraction(
                    math.comb(available, weight), math.comb(columns, weight)
                )
        total += (
            math.comb(columns, avoided)
            * pole ** (columns - avoided)
            * (1 - pole) ** avoided
            * one_block**blocks
        )
    return total


def scheduled_pole(band: int, blocks: int) -> Fraction:
    if blocks <= 5:
        return TARGET_POLE
    if blocks == 6:
        return Fraction(7, 50)
    if blocks == 7:
        return Fraction(2, 25)
    if band == 0:
        if blocks == 8:
            return Fraction(3, 100)
        if blocks == 9:
            return Fraction(1, 50)
        return Fraction(1, 100)
    if band == 1:
        if blocks <= 10:
            return Fraction(3, 100)
        if blocks == 11:
            return Fraction(1, 50)
        return Fraction(1, 100)
    if blocks == 8:
        return Fraction(3, 100)
    if blocks == 9:
        return Fraction(3, 200)
    return Fraction(1, 100)


def main() -> None:
    projection, outside_slices, inside_boundary = load_tables()
    spectrum = full_spectrum()
    measures: dict[tuple[int, Fraction], tuple[Fraction, ...]] = {}
    terms: list[Fraction] = []
    worst_by_band: list[tuple[Fraction, int, Fraction]] = []
    for band, columns in enumerate(BAND_SIZES):
        band_terms: list[tuple[Fraction, int, Fraction]] = []
        for blocks in range(2, BLOCKS_PER_TILE + 1):
            pole = scheduled_pole(band, blocks)
            key = (band, pole)
            if key not in measures:
                measures[key] = exact_capped_inside_measure(
                    band,
                    columns,
                    pole,
                    projection,
                    outside_slices,
                    inside_boundary,
                    spectrum,
                )
            moment = exact_star_moment(
                columns, blocks, pole, measures[key]
            )
            value = (
                TILES
                * math.comb(BLOCKS_PER_TILE, blocks)
                * moment
                * (TARGET_POLE / pole) ** SUPPORT_CUTOFF
                * graph_factor_upper(pole)
                * worst_puncture_block_adjustment(pole) ** blocks
            )
            terms.append(value)
            band_terms.append((value, blocks, pole))
        worst_by_band.append(max(band_terms))

    total = sum(terms, Fraction(0))
    tail_pole = Fraction(1, 20)
    tail_conversion = (TARGET_POLE / tail_pole) ** SUPPORT_CUTOFF
    tail_terms: list[Fraction] = []
    for band, columns in enumerate(BAND_SIZES):
        key = (band, tail_pole)
        if key not in measures:
            measures[key] = exact_capped_inside_measure(
                band,
                columns,
                tail_pole,
                projection,
                outside_slices,
                inside_boundary,
                spectrum,
            )
        for blocks in range(6, BLOCKS_PER_TILE + 1):
            tail_terms.append(
                TILES
                * math.comb(BLOCKS_PER_TILE, blocks)
                * exact_star_moment(
                    columns, blocks, tail_pole, measures[key]
                )
                * worst_puncture_block_adjustment(tail_pole) ** blocks
                * tail_conversion
                * graph_factor_upper(tail_pole)
            )
    common_tail = sum(tail_terms, Fraction(0))
    print("exact three-band one-tile-star certificate")
    for band, (value, blocks, pole) in enumerate(worst_by_band):
        print(
            f"band={band} worst_log2_upper={log2_fraction(value).hi} "
            f"blocks={blocks} pole={pole}"
        )
    print(f"all_bands_s_ge_2_log2_upper={log2_fraction(total).hi}")
    if total > Fraction(1, 1 << 56):
        raise SystemExit("star: all-band total exceeds 2^-56")
    print("all_bands_s_ge_2_le_2^-56=PASS")
    print(
        "common_pole_1_20_s_ge_6_log2_upper="
        f"{log2_fraction(common_tail).hi}"
    )
    if common_tail > Fraction(1, 1 << 53):
        raise SystemExit("star: common-pole size>=6 tail exceeds 2^-53")
    print("common_pole_1_20_s_ge_6_le_2^-53=PASS")
    print("status=EXACT_RATIONAL_STAR_CERTIFICATE")


if __name__ == "__main__":
    main()
