#!/usr/bin/env python3
"""Robust low-support ledger for fixed-band one-tile star families.

All selected data blocks share one tile in one of the three fixed coordinate
bands.  Pairwise cell capacity one makes their coordinates in the other two
bands collision-free.  The unknown inside/outside split enumerator is bounded
using:

* the exact complete outside-projection spectra;
* exact outside weights through eight;
* exact inside weights zero, one, and complements;
* the ordinary EBCH total spectrum; and
* exact binomial inside marginals from full projection rank.

For every inside weight, a greedy transportation cap maximizes its outside
MGF.  Independent within-band coordinate permutations then give an exact
hypergeometric union transition for that capped measure.  A second Chernoff
tilt restricts the result to total group support at most 106.

This is a long-double design diagnostic; the final certificate still needs
outward-rounded arithmetic.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from punctured_ebch_outer import full_spectrum


ROOT = Path(__file__).resolve().parent
PROJECTION_SPECTRA = ROOT / "ebch128_fixed_band_projection_spectra.csv"
OUTSIDE_SLICES = ROOT / "ebch128_fixed_band_outside_slices.csv"
INSIDE_BOUNDARY = ROOT / "ebch128_fixed_band_inside_boundary.csv"
BAND_SIZES = (42, 43, 43)
BLOCKS_PER_TILE = 64
TILES = 256
GRAPH_CODIMENSION = 24


def cluster_moments(
    columns: int, s_max: int, pole: float, weights: np.ndarray
) -> np.ndarray:
    """Union moments in the actual 42/43-column fixed band."""

    distribution = np.zeros(columns + 1, dtype=np.longdouble)
    distribution[0] = 1
    moments = np.zeros(s_max + 1, dtype=np.longdouble)
    moments[0] = 1
    powers = np.asarray(
        [np.longdouble(pole) ** union for union in range(columns + 1)]
    )
    for blocks in range(1, s_max + 1):
        following = np.zeros_like(distribution)
        for union, prior in enumerate(distribution):
            if not prior:
                continue
            for weight, mass in enumerate(weights):
                if not mass:
                    continue
                denominator = math.comb(columns, weight)
                for intersection in range(
                    max(0, union + weight - columns), min(union, weight) + 1
                ):
                    following[union + weight - intersection] += (
                        prior
                        * mass
                        * math.comb(union, intersection)
                        * math.comb(columns - union, weight - intersection)
                        / denominator
                    )
        distribution = following
        moments[blocks] = np.sum(distribution * powers)
    return moments


def load_tables():
    projection: dict[int, dict[int, int]] = defaultdict(dict)
    outside = Counter()
    inside = Counter()
    with PROJECTION_SPECTRA.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["code"] == "primal":
                projection[int(row["band"])][int(row["weight"])] = int(row["count"])
    with OUTSIDE_SLICES.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            outside[
                int(row["band"]),
                int(row["inside_weight"]),
                int(row["outside_weight"]),
            ] = int(row["count"])
    with INSIDE_BOUNDARY.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            inside[
                int(row["band"]),
                int(row["inside_weight"]),
                int(row["outside_weight"]),
            ] += int(row["count"])
    return projection, outside, inside


def capped_inside_measure(
    band: int,
    inside_size: int,
    pole: float,
    projection: dict[int, dict[int, int]],
    outside_slices: Counter,
    inside_boundary: Counter,
    spectrum: list[int],
) -> np.ndarray:
    outside_size = 128 - inside_size
    result = np.zeros(inside_size + 1, dtype=np.longdouble)
    for inside_weight in range(inside_size + 1):
        column_mass = (
            math.comb(inside_size, inside_weight) * (1 << (64 - inside_size))
            - (inside_weight == 0)
            - (inside_weight == inside_size)
        )
        if inside_weight in (0, 1, inside_size - 1, inside_size):
            recovered = 0
            for (row_band, row_inside, outside_weight), raw_count in inside_boundary.items():
                if row_band != band or row_inside != inside_weight:
                    continue
                count = raw_count
                if (inside_weight, outside_weight) in (
                    (0, 0),
                    (inside_size, outside_size),
                ):
                    count -= 1
                recovered += count
                result[inside_weight] += (
                    np.longdouble(count) * np.longdouble(pole) ** outside_weight
                )
            if recovered != column_mass:
                raise SystemExit("fixed star: incomplete exact inside boundary row")
            continue

        recovered = 0
        for outside_weight in range(1, 9):
            count = outside_slices.get(
                (band, inside_weight, outside_weight), 0
            )
            recovered += count
            result[inside_weight] += (
                np.longdouble(count) * np.longdouble(pole) ** outside_weight
            )
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
            result[inside_weight] += (
                np.longdouble(take) * np.longdouble(pole) ** outside_weight
            )
            remaining -= take
        if remaining:
            raise SystemExit("fixed star: transportation caps did not cover a column")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group-pole", type=float, default=0.181)
    parser.add_argument("--support-cutoff", type=int, default=106)
    parser.add_argument(
        "--tilts",
        default=".181,.16,.14,.12,.1,.08,.06,.05,.04,.03,.025,.02,.015,.01",
    )
    args = parser.parse_args()
    tilts = [float(value) for value in args.tilts.split(",")]
    if not 0 < min(tilts) <= max(tilts) <= args.group_pole < 1:
        raise SystemExit("fixed star: invalid poles")

    projection, outside_slices, inside_boundary = load_tables()
    spectrum = full_spectrum()
    candidates = []
    for band, inside_size in enumerate(BAND_SIZES):
        band_candidates = []
        for tilt in tilts:
            measure = capped_inside_measure(
                band,
                inside_size,
                tilt,
                projection,
                outside_slices,
                inside_boundary,
                spectrum,
            )
            mass = float(np.sum(measure))
            clusters = cluster_moments(
                inside_size, BLOCKS_PER_TILE, tilt, measure / mass
            )
            band_candidates.append((tilt, mass, clusters))
        candidates.append(band_candidates)

    all_rows = []
    print("fixed three-band star low-support ledger")
    print(f"group_pole={args.group_pole:.12f}")
    print(f"support_cutoff={args.support_cutoff}")
    for band, band_candidates in enumerate(candidates):
        band_rows = []
        for blocks in range(1, BLOCKS_PER_TILE + 1):
            rows = []
            for tilt, mass, clusters in band_candidates:
                value = (
                    math.log2(TILES)
                    + math.log2(math.comb(BLOCKS_PER_TILE, blocks))
                    - GRAPH_CODIMENSION
                    + blocks * math.log2(mass)
                    + math.log2(float(clusters[blocks]))
                    + args.support_cutoff * math.log2(args.group_pole / tilt)
                )
                rows.append((value, tilt))
            best = min(rows)
            all_rows.append((best[0], band, blocks, best[1]))
            band_rows.append((blocks, best[1]))
        worst_multi = max(row for row in all_rows if row[1] == band and row[2] >= 2)
        print(
            f"band={band} worst_s_ge_2_log2={worst_multi[0]:.12f} "
            f"at_s={worst_multi[2]} tilt={worst_multi[3]:.12f}"
        )
        runs = []
        run_start, run_pole = band_rows[1]
        previous = run_start
        for blocks, pole in band_rows[2:]:
            if pole != run_pole:
                runs.append((run_start, previous, run_pole))
                run_start, run_pole = blocks, pole
            previous = blocks
        runs.append((run_start, previous, run_pole))
        print(f"band={band} s_ge_2_pole_runs={runs}")
    multi_total_rows = [row[0] for row in all_rows if row[2] >= 2]
    maximum = max(multi_total_rows)
    multi_total = maximum + math.log2(
        sum(2 ** (value - maximum) for value in multi_total_rows)
    )
    print(f"all_bands_s_ge_2_log2={multi_total:.12f}")
    print("one-block rows are deliberately replaced by the exact PA1 ledger")
    print("status=DIAGNOSTIC_NEEDS_OUTWARD_ARITHMETIC")


if __name__ == "__main__":
    main()
