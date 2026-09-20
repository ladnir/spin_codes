#!/usr/bin/env python3
"""Merge Arb receipts into one complete sector-zero diagonal receipt."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import arb, ctx


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "sector_zero_hoeffding_radial_complete.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--precision", type=int, default=2048)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.precision < 128:
        parser.error("precision must be at least 128 bits")
    ctx.prec = args.precision
    ctx.threads = 1

    common = None
    selected: dict[int, tuple[arb, Path]] = {}
    skipped_nonfinite = []
    for path in args.inputs:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != "pure-ea-sector-zero-hoeffding-radial-arb-v1":
            raise ValueError(f"unexpected schema in {path}")
        parameters = payload["parameters"]
        identity = {
            name: parameters[name]
            for name in (
                "message_bits",
                "output_bits",
                "right_degree",
                "level",
            )
        }
        if common is None:
            common = identity
        elif identity != common:
            raise ValueError(f"parameter mismatch in {path}")
        for row in payload["orders"]:
            order = int(row["order"])
            value = arb(row["scaled_diagonal_contribution"])
            if not value.is_finite():
                skipped_nonfinite.append(
                    {"order": order, "source": str(path), "value": str(value)}
                )
                continue
            if order in selected:
                previous, previous_path = selected[order]
                if not previous.overlaps(value):
                    raise ArithmeticError(
                        f"disjoint duplicate order {order}: "
                        f"{previous_path} gives {previous}, {path} gives {value}"
                    )
                if value.rad() < previous.rad():
                    selected[order] = (value, path)
            else:
                selected[order] = (value, path)

    assert common is not None
    required = list(range(1, common["output_bits"] + 1))
    missing = [order for order in required if order not in selected]
    if missing:
        raise ValueError(f"missing finite order contributions: {missing}")

    total = sum((selected[order][0] for order in required), arb(0))
    if not total.is_finite():
        raise ArithmeticError(f"non-finite complete diagonal: {total}")
    chosen_rows = [
        {
            "order": order,
            "source": str(selected[order][1]),
            "scaled_diagonal_contribution": str(selected[order][0]),
        }
        for order in required
    ]
    output = {
        "schema": "pure-ea-sector-zero-hoeffding-radial-merged-v1",
        "status": "OUTWARD_COMPLETE_SECTOR_ZERO_DIAGONAL",
        "parameters": {**common, "precision_bits": args.precision},
        "scaled_sector_zero_diagonal": str(total),
        "orders": chosen_rows,
        "source_receipts": [str(path) for path in args.inputs],
        "skipped_nonfinite_rows": skipped_nonfinite,
        "scope": [
            "Every Hoeffding order from one through output_bits has one finite Arb enclosure.",
            "Their Arb sum encloses the complete sector-zero diagonal at the reported output level.",
            "This receipt does not bound other levels, sectors one and two, complete accumulator shells, or rank failure.",
        ],
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(
        f"receipt,{args.output},status,{output['status']},scaled,{total}",
        flush=True,
    )


if __name__ == "__main__":
    main()
