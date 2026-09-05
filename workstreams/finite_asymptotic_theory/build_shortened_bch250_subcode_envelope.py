#!/usr/bin/env python3
"""Restrict the audited BCH250 packing envelope to a lower-dimensional subcode."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_PARENT = WORKSTREAM / "shortened_bch250_125_constant_weight_envelope.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", type=Path, default=DEFAULT_PARENT)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.dimension <= 125:
        parser.error("dimension must lie in [1,125]")

    parent = json.loads(args.parent.read_text(encoding="utf-8"))
    total_nonzero = (1 << args.dimension) - 1
    spectrum = []
    for row in parent["spectrum"]:
        updated = dict(row)
        if int(updated["weight"]) == 0:
            updated["multiplicity_upper"] = 1
        else:
            updated["multiplicity_upper"] = min(
                int(updated["multiplicity_upper"]), total_nonzero
            )
        spectrum.append(updated)

    payload = {
        "schema": "shortened-bch-subcode-constant-weight-envelope-v1",
        "candidate": f"fixed [{parent['parameters']['outer_bits']},{args.dimension},>=38] subcode of audited BCH250",
        "parameters": {
            **parent["parameters"],
            "outer_dimension": args.dimension,
            "parent_dimension": int(parent["parameters"]["outer_dimension"]),
            "exact_nonzero_mass": total_nonzero,
        },
        "spectrum": spectrum,
        "scope": (
            "Exact pointwise envelope for every fixed subcode of the audited "
            "BCH250 code with the stated dimension. The bounds combine the "
            "parent packing caps with the subcode cardinality. They do not "
            "claim an exact subcode spectrum."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
