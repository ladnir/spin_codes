#!/usr/bin/env python3
"""Transportation-LP probe for fixed three-band BCH correlations.

The unknown exact table N[a,b,c] is constrained by every exact projection
fact currently available:

* the ordinary [128,64,22] weight spectrum (diagonal sums);
* the exact binomial marginal of each full-rank fixed-band projection;
* the exact complement-projection spectra;
* all exact complement-weight rows through eight; and
* all exact inside-weight 0/1 rows and their complements.

For a pair collision in a band, the within-band random coordinate permutation
gives a hypergeometric union kernel K(a,b).  This kernel is positive
semidefinite, hence K(a,b) <= sqrt(K(a,a) K(b,b)).  A four-block collision
motif consequently decouples into four one-codeword factors.  This script
maximizes each such factor over all nonnegative N consistent with the exact
tables.  It is a floating solver probe.  A useful result must later be
exported as a dual solution and verified with outward/exact arithmetic.

The optional ``highspy`` package is deliberately not a proof dependency; it
is only used to discover a compact dual certificate.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    import highspy
except ImportError:
    highspy = None

from certify_bch_band_projections import BANDS
from punctured_ebch_outer import full_spectrum


ROOT = Path(__file__).resolve().parent
PROJECTION_SPECTRA = ROOT / "ebch128_fixed_band_projection_spectra.csv"
OUTSIDE_SLICES = ROOT / "ebch128_fixed_band_outside_slices.csv"
INSIDE_BOUNDARY = ROOT / "ebch128_fixed_band_inside_boundary.csv"
BAND_SIZES = tuple(end - begin for begin, end in BANDS)
CODE_SIZE = 1 << 64
OBJECTIVE_SCALE_LOG2 = 80


def load_exact_tables():
    outside_spectra: dict[int, dict[int, int]] = defaultdict(dict)
    outside_slices: dict[tuple[int, int, int], int] = {}
    inside_boundary: dict[tuple[int, int, int], int] = {}
    with PROJECTION_SPECTRA.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["code"] == "primal":
                outside_spectra[int(row["band"])][int(row["weight"])] = int(
                    row["count"]
                )
    with OUTSIDE_SLICES.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (
                int(row["band"]),
                int(row["inside_weight"]),
                int(row["outside_weight"]),
            )
            outside_slices[key] = int(row["count"])
    with INSIDE_BOUNDARY.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (
                int(row["band"]),
                int(row["inside_weight"]),
                int(row["outside_weight"]),
            )
            inside_boundary[key] = inside_boundary.get(key, 0) + int(row["count"])
    return outside_spectra, outside_slices, inside_boundary


def triples() -> list[tuple[int, int, int]]:
    return [
        (first, second, third)
        for first in range(BAND_SIZES[0] + 1)
        for second in range(BAND_SIZES[1] + 1)
        for third in range(BAND_SIZES[2] + 1)
    ]


def hypergeometric_pair_diagonal(columns: int, weight: int, pole: float) -> float:
    denominator = math.comb(columns, weight)
    return sum(
        math.comb(weight, intersection)
        * math.comb(columns - weight, weight - intersection)
        / denominator
        * pole ** (2 * weight - intersection)
        for intersection in range(max(0, 2 * weight - columns), weight + 1)
    )


def factor_tables(pole: float) -> tuple[tuple[np.ndarray, ...], ...]:
    singleton = tuple(
        np.array([pole**weight for weight in range(size + 1)])
        for size in BAND_SIZES
    )
    paired = tuple(
        np.array(
            [
                math.sqrt(hypergeometric_pair_diagonal(size, weight, pole))
                for weight in range(size + 1)
            ]
        )
        for size in BAND_SIZES
    )
    # Dense-band matching in band zero; one extra pair in each other band.
    # The four rows are the four vertex types up to relabeling.
    return (
        (paired[0], paired[1], paired[2]),
        (paired[0], singleton[1], singleton[2]),
        (paired[0], paired[1], singleton[2]),
        (paired[0], singleton[1], paired[2]),
    )


def add_equality(
    rows: list[tuple[str, int, list[int]]], name: str, rhs: int, indices: list[int]
) -> None:
    if not indices and rhs:
        raise SystemExit(f"transport LP: empty nonzero row {name}")
    rows.append((name, rhs, indices))


def exact_rows(
    states: list[tuple[int, int, int]],
) -> list[tuple[str, int, list[int]]]:
    outside_spectra, outside_slices, inside_boundary = load_exact_tables()
    by_total: dict[int, list[int]] = defaultdict(list)
    by_inside: list[dict[int, list[int]]] = [defaultdict(list) for _ in range(3)]
    by_inside_outside: list[dict[tuple[int, int], list[int]]] = [
        defaultdict(list) for _ in range(3)
    ]
    for index, state in enumerate(states):
        by_total[sum(state)].append(index)
        for band in range(3):
            inside = state[band]
            outside = sum(state) - inside
            by_inside[band][inside].append(index)
            by_inside_outside[band][inside, outside].append(index)

    rows: list[tuple[str, int, list[int]]] = []
    spectrum = full_spectrum()
    for total in range(129):
        add_equality(rows, f"total_{total}", spectrum[total], by_total[total])

    for band, inside_size in enumerate(BAND_SIZES):
        outside_size = 128 - inside_size
        for inside in range(inside_size + 1):
            rhs = math.comb(inside_size, inside) << (64 - inside_size)
            add_equality(
                rows,
                f"band_{band}_inside_{inside}",
                rhs,
                by_inside[band][inside],
            )
        for outside in range(outside_size + 1):
            indices = [
                index
                for (inside, row_outside), cell in by_inside_outside[band].items()
                if row_outside == outside
                for index in cell
            ]
            add_equality(
                rows,
                f"band_{band}_outside_{outside}",
                outside_spectra[band].get(outside, 0),
                indices,
            )

        boundary_inside = (0, 1, inside_size - 1, inside_size)
        for inside in range(inside_size + 1):
            outside_range = (
                range(outside_size + 1) if inside in boundary_inside else range(9)
            )
            for outside in outside_range:
                source = (
                    inside_boundary
                    if inside in boundary_inside
                    else outside_slices
                )
                add_equality(
                    rows,
                    f"band_{band}_cell_{inside}_{outside}",
                    source.get((band, inside, outside), 0),
                    by_inside_outside[band].get((inside, outside), []),
                )
    return rows


def remove_trivial_codewords(
    states: list[tuple[int, int, int]], rows: list[tuple[str, int, list[int]]]
) -> list[tuple[str, int, list[int]]]:
    endpoints = {
        states.index((0, 0, 0)),
        states.index(tuple(BAND_SIZES)),
    }
    result = []
    for name, rhs, indices in rows:
        removed = sum(index in endpoints for index in indices)
        result.append(
            (name, rhs - removed, [index for index in indices if index not in endpoints])
        )
    return result


def build_model(
    states: list[tuple[int, int, int]],
    rows: list[tuple[str, int, list[int]]],
    solver_output: bool = False,
) -> tuple[highspy.Highs, np.ndarray]:
    if highspy is None:
        raise SystemExit(
            "transport LP probe requires highspy; install it in a temporary "
            "tool environment and add that directory to PYTHONPATH"
        )
    # Write N_i = cap_i z_i, 0 <= z_i <= 1, where cap_i is the smallest
    # exact aggregate containing cell i.  Every normalized row coefficient
    # cap_i/rhs then lies in [0,1], avoiding the 1..2^64 dynamic range of the
    # raw-count formulation.
    caps = np.full(len(states), CODE_SIZE, dtype=object)
    covered = np.zeros(len(states), dtype=bool)
    for _name, rhs, row_indices in rows:
        for index in row_indices:
            covered[index] = True
            if rhs < caps[index]:
                caps[index] = rhs
    # ``remove_trivial_codewords`` removes the zero and all-ones cells from
    # every equality while deliberately retaining stable state indices.  Such
    # uncovered columns are fixed to zero here; otherwise they are free
    # variables, and the zero cell in particular makes every positive-objective
    # transport model unbounded.
    caps[~covered] = 0
    caps_float = np.asarray(caps, dtype=float)

    model = highspy.Highs()
    model.setOptionValue("output_flag", solver_output)
    model.setOptionValue("presolve", "on")
    columns = len(states)
    model.addCols(
        columns,
        np.zeros(columns),
        np.zeros(columns),
        np.ones(columns),
        0,
        np.zeros(columns + 1, dtype=np.int32),
        np.zeros(0, dtype=np.int32),
        np.zeros(0),
    )

    starts = np.zeros(len(rows) + 1, dtype=np.int32)
    indices: list[int] = []
    values: list[float] = []
    positive_rows = [row for row in rows if row[1] > 0]
    starts = np.zeros(len(positive_rows) + 1, dtype=np.int32)
    lower = np.ones(len(positive_rows))
    upper = np.ones(len(positive_rows))
    for row_index, (_name, rhs, row_indices) in enumerate(positive_rows):
        starts[row_index] = len(indices)
        for index in row_indices:
            coefficient = caps_float[index] / rhs
            if coefficient:
                indices.append(index)
                values.append(coefficient)
    starts[len(positive_rows)] = len(indices)
    model.addRows(
        len(positive_rows),
        lower,
        upper,
        len(indices),
        starts,
        np.asarray(indices, dtype=np.int32),
        np.asarray(values),
    )
    model.changeObjectiveSense(highspy.ObjSense.kMaximize)
    return model, caps_float


def main() -> None:
    if highspy is None:
        raise SystemExit(
            "transport LP probe requires highspy; install it in a temporary "
            "tool environment and add that directory to PYTHONPATH"
        )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group-pole", type=float, default=0.181)
    parser.add_argument("--solver-output", action="store_true")
    args = parser.parse_args()
    if not 0 < args.group_pole < 1:
        raise SystemExit("transport LP: invalid pole")

    states = triples()
    print(f"building_exact_rows states={len(states)}", flush=True)
    rows = remove_trivial_codewords(states, exact_rows(states))
    print(f"building_model equalities={len(rows)}", flush=True)
    model, caps = build_model(states, rows, args.solver_output)
    print("fixed three-band BCH transportation LP")
    print(
        f"states={len(states)} equalities={len(rows)} "
        f"group_pole={args.group_pole:.12f}"
    )
    for factor_index, factors in enumerate(factor_tables(args.group_pole)):
        costs = np.array(
            [
                factors[0][state[0]]
                * factors[1][state[1]]
                * factors[2][state[2]]
                for state in states
            ]
        )
        costs *= caps / CODE_SIZE * (1 << OBJECTIVE_SCALE_LOG2)
        model.changeColsCost(
            len(states), np.arange(len(states), dtype=np.int32), costs
        )
        print(f"solving_vertex_type={factor_index}", flush=True)
        model.run()
        status = model.getModelStatus()
        if status != highspy.HighsModelStatus.kOptimal:
            raise SystemExit(f"transport LP: solver status {status}")
        scaled_normalized = model.getObjectiveValue()
        normalized_log2 = math.log2(scaled_normalized) - OBJECTIVE_SCALE_LOG2
        raw_log2 = normalized_log2 + 64
        print(
            f"vertex_type={factor_index} "
            f"nontrivial_raw_log2={raw_log2:.12f} "
            f"normalized_log2={normalized_log2:.12f}"
        )
    print("status=FLOATING_LP_PROBE_NEEDS_DUAL_CERTIFICATE")


if __name__ == "__main__":
    main()
