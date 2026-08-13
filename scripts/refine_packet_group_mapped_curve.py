#!/usr/bin/env python3
"""Apply the same shaped-density refinement to mapped hard-family rows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from packet_group_profile_bound import refine_inner_shaped_density, split_cap_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--groups", default="1,2,4,8,16,32,64")
    parser.add_argument("--coordinate-iterations", type=int, default=3)
    parser.add_argument("--passes", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    selected = {int(value) for value in args.groups.split(",")}
    source = []
    for path in args.input:
        source.extend(json.loads(path.read_text(encoding="utf-8")))
    rows = []
    split_caps = split_cap_table()
    for row in sorted(source, key=lambda item: item["group_bits"]):
        group_bits = int(row["group_bits"])
        if group_bits not in selected:
            continue
        refined = refine_inner_shaped_density(
            group_bits,
            row["profile"],
            split_caps,
            initial_pole=float(row["pole"]),
            initial_exponent=float(row["density_exponent"]),
            coordinate_iterations=args.coordinate_iterations,
            passes=args.passes,
            witness_iterations=40,
        )
        combined = float(row["outer_profile_log2"]) + refined[0]
        result = {
            **row,
            "unrefined_combined_log2": row["combined_log2"],
            "unrefined_margin_log2": row["diagnostic_margin_log2"],
            "inner_probability_log2": refined[0],
            "combined_log2": combined,
            "diagnostic_margin_log2": float(row["union_target_log2"]) - combined,
            "pole": refined[1],
            "fugacities": refined[2].tolist(),
            "inner_details": refined[3],
            "shape_parameters": refined[4].tolist(),
            "refinement_evaluations": refined[5],
        }
        rows.append(result)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
        print(
            f"g={group_bits} margin={result['diagnostic_margin_log2']:.6f} "
            f"combined={combined:.6f} pole={refined[1]:.9f} "
            f"parameters={','.join(f'{value:.6f}' for value in refined[4][:3])}",
            flush=True,
        )
    print("status=DIAGNOSTIC_MAPPED_HARD_FAMILY_SHAPED_REFINEMENT")


if __name__ == "__main__":
    main()
