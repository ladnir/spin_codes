#!/usr/bin/env python3
"""Export a native PB-XOR shell-counting instance for ApproxMCPB.

The Boolean variables are the coordinates of one extended primitive BCH word.
An exact pseudo-Boolean equality fixes its weight.  Native XOR constraints
enforce membership in the code.  For affine-invariant codes, fixing two
coordinates to one reduces the count by the exact 2-design incidence factor.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from approx_count_low_shells import systematic_generator_rows
from prepare_wambach_bases import CODE_SPECS, generator_polynomial_generic


ROOT = Path(__file__).resolve().parents[1]
EXACT_MINIMUM_SHELLS = {
    "b8": (4, 14),
    "b32": (8, 620),
    "b128": (22, 243840),
}


def extended_standard_rows(code_name: str) -> list[int]:
    spec = CODE_SPECS[code_name]
    m = spec["m"]
    punctured_length = (1 << m) - 1
    generator = generator_polynomial_generic(m, spec["delta"], spec["modulus"])
    rows = []
    for shift in range(spec["dimension"]):
        punctured = generator << shift
        rows.append(punctured | ((punctured.bit_count() & 1) << punctured_length))
    return systematic_generator_rows(rows)


def export_instance(code_name: str, weight: int, output: Path) -> dict:
    rows = extended_standard_rows(code_name)
    dimension = len(rows)
    length = 1 << CODE_SPECS[code_name]["m"]
    pivots = [(row & -row).bit_length() - 1 for row in rows]
    if len(set(pivots)) != dimension:
        raise ArithmeticError("systematic pivots are not unique")
    pivot_set = set(pivots)

    xor_constraints: list[list[int]] = []
    for coordinate in range(length):
        if coordinate in pivot_set:
            continue
        variables = [coordinate + 1]
        variables.extend(
            pivots[row_index] + 1
            for row_index, row in enumerate(rows)
            if (row >> coordinate) & 1
        )
        xor_constraints.append(variables)

    lines = [
        f"* #variable= {length} #constraint= {len(xor_constraints) + 3}",
        "* Native PB-XOR extended-BCH shell instance",
        "* ind " + " ".join(str(coordinate + 1) for coordinate in pivots) + " 0",
    ]
    lines.extend(
        "* xor " + " ".join(f"x{variable}" for variable in variables) + " 0"
        for variables in xor_constraints
    )
    lines.append(
        " ".join(f"+1 x{index}" for index in range(1, length + 1))
        + f" = {weight};"
    )
    lines.extend(("+1 x1 = 1;", "+1 x2 = 1;"))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="ascii")

    digest = hashlib.sha256(
        b"".join(row.to_bytes((length + 7) // 8, "little") for row in rows)
    ).hexdigest()
    metadata = {
        "classification": "native PB-XOR projected model-counting instance",
        "code": code_name,
        "parameters": [length, dimension],
        "designed_distance": CODE_SPECS[code_name]["delta"],
        "target_weight": weight,
        "fixed_coordinates": [0, 1],
        "native_xor_constraints": len(xor_constraints),
        "pseudo_boolean_constraints": 3,
        "variables": length,
        "projection_variables": [coordinate + 1 for coordinate in pivots],
        "projection_size": dimension,
        "generator_rows_sha256": digest,
        "opb_path": str(output.resolve()),
    }
    if code_name in EXACT_MINIMUM_SHELLS:
        exact_weight, exact_shell = EXACT_MINIMUM_SHELLS[code_name]
        if weight == exact_weight:
            exact_fixed_pair = exact_shell * weight * (weight - 1) // (
                length * (length - 1)
            )
            metadata["exact_calibration"] = {
                "full_shell": exact_shell,
                "fixed_pair_count": exact_fixed_pair,
            }
    output.with_suffix(".json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--code", choices=tuple(CODE_SPECS), required=True)
    parser.add_argument("--weight", type=int, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or (
        ROOT / "generated" / f"approxmcpb_{args.code}_w{args.weight}_pair.opb"
    )
    print(json.dumps(export_instance(args.code, args.weight, output), indent=2))


if __name__ == "__main__":
    main()
