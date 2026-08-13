#!/usr/bin/env python3
"""Export an exact rational transportation LP for one BCH degree type.

The older floating probe normalized every cell by its smallest aggregate.
That creates matrix coefficients below ``1e-18``, which floating LP solvers
silently discard.  This exporter instead keeps the exact integer codeword
counts ``N[a,b,c]`` as variables and every projection equality with unit
coefficients.

The rooted cluster moments are rounded upward to exact dyadic rationals by
``three_band_exact_moments``.  Their product is rounded upward once more to an
integer objective grid.  Consequently the exported continuous LP is a safe
relaxation with an entirely rational input suitable for exact SoPlex/SCIP.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from pathlib import Path

from outward_log2 import log2_fraction
from probe_bch_three_band_cell_cap_motifs import cell_caps
from probe_bch_three_band_transport_lp import (
    BAND_SIZES,
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


def ceil_fraction(value: Fraction) -> int:
    return -(-value.numerator // value.denominator)


def write_terms(handle, terms: list[str], *, width: int = 8) -> None:
    for begin in range(0, len(terms), width):
        handle.write(" ")
        handle.write(" ".join(terms[begin : begin + width]))
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degrees", type=parse_degrees, required=True)
    parser.add_argument("--pole", type=Fraction, default=Fraction(1, 20))
    parser.add_argument("--root-bits", type=int, default=160)
    parser.add_argument("--objective-bits", type=int, default=256)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not 0 < args.pole < 1 or args.root_bits <= 0 or args.objective_bits <= 0:
        raise SystemExit("exact transport exporter: invalid arguments")

    states = triples()
    rows = remove_trivial_codewords(states, exact_rows(states))
    caps, _diagonal_mass, _by_total = cell_caps(states, rows)
    active = {index for index, cap in enumerate(caps) if cap > 0}
    # The endpoints were removed from every equality.  cell_caps initializes
    # unseen cells to 2^64, so remove every index that is absent from all rows.
    covered = {index for _name, _rhs, indices in rows for index in indices}
    active &= covered

    tables = tuple(
        factor_table_upper(
            BAND_SIZES[band],
            args.degrees[band],
            args.pole,
            root_bits=args.root_bits,
        )
        for band in range(3)
    )
    objective_scale = 1 << args.objective_bits
    costs: dict[int, int] = {}
    maximum_rounding_error = Fraction(0)
    for index in active:
        first, second, third = states[index]
        score = tables[0][first] * tables[1][second] * tables[2][third]
        cost = ceil_fraction(score * objective_scale)
        costs[index] = cost
        maximum_rounding_error += Fraction(caps[index], objective_scale)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("Maximize\n")
        # Keep coefficients below SCIP's floating infinity threshold while
        # preserving the exact dyadic value.  In exact mode the LP reader
        # parses these rational tokens without converting them to binary64.
        objective = [
            f"+ {cost}/{objective_scale} x{index}" for index, cost in costs.items()
        ]
        write_terms(handle, objective, width=2)
        handle.write("Subject To\n")
        written_rows = 0
        for name, rhs, indices in rows:
            if rhs <= 0:
                continue
            present = [index for index in indices if index in active]
            if not present:
                raise SystemExit(f"positive row {name} has no active cells")
            handle.write(f" {name}:\n")
            write_terms(handle, [f"+ x{index}" for index in present], width=16)
            handle.write(f" = {rhs}\n")
            written_rows += 1
        handle.write("Bounds\n")
        for index in sorted(active):
            handle.write(f" 0 <= x{index} <= {caps[index]}\n")
        handle.write("End\n")

    print("exact three-band BCH transportation LP export")
    print(f"degrees={args.degrees} pole={args.pole}")
    print(
        f"states={len(states)} active_variables={len(active)} "
        f"positive_equalities={written_rows}"
    )
    print(
        f"root_bits={args.root_bits} objective_bits={args.objective_bits} "
        f"objective_rounding_error_log2_upper="
        f"{log2_fraction(maximum_rounding_error).hi}"
    )
    print(f"out={args.out}")
    print("status=EXACT_RATIONAL_SAFE_RELAXATION")


if __name__ == "__main__":
    main()
