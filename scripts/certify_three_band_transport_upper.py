#!/usr/bin/env python3
"""Exact-dual certificate for a stable three-band transport relaxation.

For exact codeword counts ``N_i`` write ``N_i = cap_i z_i``.  Every exact
aggregate row implies

    sum_i (cap_i / rhs) z_i <= 1,       0 <= z_i <= 1.

Terms below a selected threshold may be omitted: all coefficients are
nonnegative, so this only weakens the inequality.  The resulting LP is well
conditioned enough to discover a dual with HiGHS.  The solver output is *not*
trusted.  Row multipliers are rounded upward to dyadic rationals, and exact
fraction arithmetic constructs the nonnegative bound slacks

    u_i = max(0, c_i - sum_r A[r,i] y_r).

Then ``sum(y)+sum(u)`` is an independently checked exact dual upper bound.
The objective coefficients use exact BCH diagonal moments and outward dyadic
roots, so the complete reported factor is safe in the theorem direction.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

try:
    import highspy
except ImportError as error:
    raise SystemExit(
        "transport upper certificate requires highspy in a temporary tool environment"
    ) from error

from outward_log2 import log2_fraction
from probe_bch_three_band_cell_cap_motifs import cell_caps
from probe_bch_three_band_transport_lp import (
    BAND_SIZES,
    CODE_SIZE,
    exact_rows,
    remove_trivial_codewords,
    triples,
)
from three_band_exact_moments import factor_table_upper


def parse_degrees(text: str) -> tuple[int, int, int]:
    values = tuple(int(field) for field in text.split(","))
    if len(values) != 3 or any(not 1 <= value <= 64 for value in values):
        raise argparse.ArgumentTypeError("degrees must be d0,d1,d2 in 1..64")
    return values  # type: ignore[return-value]


def dyadic_ceiling(value: float, bits: int) -> Fraction:
    if not math.isfinite(value):
        raise ValueError("nonfinite dual multiplier")
    if value <= 0:
        return Fraction(0)
    scale = 1 << bits
    return Fraction(math.ceil(value * scale), scale)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degrees", type=parse_degrees, required=True)
    parser.add_argument("--pole", type=Fraction, default=Fraction(1, 20))
    parser.add_argument("--coefficient-threshold", type=Fraction, default=Fraction(1, 10**9))
    parser.add_argument("--root-bits", type=int, default=160)
    parser.add_argument("--dual-bits", type=int, default=160)
    parser.add_argument(
        "--solver-scale-bits",
        type=int,
        default=-1,
        help="binary objective scale for discovery, or -1 to choose it automatically",
    )
    parser.add_argument("--solver-output", action="store_true")
    parser.add_argument("--manifest", type=Path, default=None)
    args = parser.parse_args()
    if not (
        0 < args.pole < 1
        and 0 < args.coefficient_threshold < 1
        and args.root_bits > 0
        and args.dual_bits > 0
        and -1 <= args.solver_scale_bits <= 100
    ):
        raise SystemExit("transport upper certificate: invalid arguments")

    states = triples()
    rows = remove_trivial_codewords(states, exact_rows(states))
    caps, _diagonal_mass, _by_total = cell_caps(states, rows)
    covered = {index for _name, _rhs, indices in rows for index in indices}
    active = [index for index, cap in enumerate(caps) if cap > 0 and index in covered]
    active_set = set(active)

    factor_tables = tuple(
        factor_table_upper(
            BAND_SIZES[band],
            args.degrees[band],
            args.pole,
            root_bits=args.root_bits,
        )
        for band in range(3)
    )
    exact_costs = [Fraction(0)] * len(states)
    for index in active:
        first, second, third = states[index]
        exact_costs[index] = (
            factor_tables[0][first]
            * factor_tables[1][second]
            * factor_tables[2][third]
            * Fraction(caps[index], CODE_SIZE)
        )

    row_entries: list[list[tuple[int, Fraction]]] = []
    row_names: list[str] = []
    for name, rhs, indices in rows:
        if rhs <= 0:
            continue
        entries = []
        for index in indices:
            if index not in active_set:
                continue
            coefficient = Fraction(caps[index], rhs)
            if coefficient >= args.coefficient_threshold:
                entries.append((index, coefficient))
        if entries:
            row_names.append(name)
            row_entries.append(entries)

    model = highspy.Highs()
    model.setOptionValue("output_flag", args.solver_output)
    model.setOptionValue("solver", "ipm")
    model.setOptionValue("run_crossover", "on")
    maximum_cost = max(exact_costs)
    if args.solver_scale_bits < 0:
        # Put the largest discovery coefficient near 2^10 while respecting
        # HiGHS' 1e20 infinity threshold.  This avoids both underflow for
        # strongly decaying profiles and overflow for asymmetric dense ones.
        maximum_log2 = math.log2(float(maximum_cost))
        solver_scale_bits = max(0, min(100, math.floor(10 - maximum_log2)))
    else:
        solver_scale_bits = args.solver_scale_bits
    scale = 1 << solver_scale_bits
    costs_float = np.array([float(cost * scale) for cost in exact_costs])
    columns = len(states)
    model.addCols(
        columns,
        costs_float,
        np.zeros(columns),
        np.ones(columns),
        0,
        np.zeros(columns + 1, dtype=np.int32),
        np.zeros(0, dtype=np.int32),
        np.zeros(0),
    )
    starts = [0]
    indices: list[int] = []
    values: list[float] = []
    for entries in row_entries:
        indices.extend(index for index, _coefficient in entries)
        values.extend(float(coefficient) for _index, coefficient in entries)
        starts.append(len(indices))
    model.addRows(
        len(row_entries),
        np.full(len(row_entries), -highspy.kHighsInf),
        np.ones(len(row_entries)),
        len(indices),
        np.asarray(starts, dtype=np.int32),
        np.asarray(indices, dtype=np.int32),
        np.asarray(values),
    )
    model.changeObjectiveSense(highspy.ObjSense.kMaximize)
    model.run()
    if model.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        raise SystemExit(f"transport upper certificate: {model.getModelStatus()}")

    solution = model.getSolution()
    # For a maximization LP HiGHS reports nonnegative multipliers on upper
    # rows.  Divide out the discovery objective scale before rounding.
    discovered = [dual / scale for dual in solution.row_dual]
    if min(discovered, default=0.0) < -1e-8:
        raise SystemExit("transport upper certificate: unexpected negative row dual")
    dual = [dyadic_ceiling(max(0.0, value), args.dual_bits) for value in discovered]

    coverage = [Fraction(0)] * len(states)
    for multiplier, entries in zip(dual, row_entries, strict=True):
        if not multiplier:
            continue
        for index, coefficient in entries:
            coverage[index] += multiplier * coefficient
    slacks = [Fraction(0)] * len(states)
    for index in active:
        slacks[index] = max(Fraction(0), exact_costs[index] - coverage[index])

    dual_bound_normalized = sum(dual, Fraction(0)) + sum(slacks, Fraction(0))
    for index in active:
        if coverage[index] + slacks[index] < exact_costs[index]:
            raise AssertionError(f"dual coverage failed at cell {index}")
    raw_factor_bound = dual_bound_normalized * CODE_SIZE
    raw_log = log2_fraction(raw_factor_bound)
    nonzero_rows = sum(bool(value) for value in dual)
    nonzero_slacks = sum(bool(value) for value in slacks)

    payload = {
        "degrees": args.degrees,
        "pole": f"{args.pole.numerator}/{args.pole.denominator}",
        "coefficient_threshold": (
            f"{args.coefficient_threshold.numerator}/"
            f"{args.coefficient_threshold.denominator}"
        ),
        "root_bits": args.root_bits,
        "dual_bits": args.dual_bits,
        "solver_scale_bits": solver_scale_bits,
        "states": len(states),
        "active_variables": len(active),
        "retained_rows": len(row_entries),
        "retained_nonzeros": len(indices),
        "nonzero_dual_rows": nonzero_rows,
        "nonzero_bound_slacks": nonzero_slacks,
        "raw_factor_log2_lower": str(raw_log.lo),
        "raw_factor_log2_upper": str(raw_log.hi),
        "status": "EXACT_DYADIC_DUAL_VERIFIED",
    }
    if args.manifest is not None:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    print("three-band aggregate-upper transport certificate")
    for key, value in payload.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
