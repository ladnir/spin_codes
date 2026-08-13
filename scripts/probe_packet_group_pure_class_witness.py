#!/usr/bin/env python3
"""Build a sharp or class-mass-floor witness for one pure packet class."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from packet_group_outer_profile import atom_count, optimize_spectrum_outer
from packet_group_profile_bound import inner_probability, split_cap_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--class-weight", type=int, required=True)
    parser.add_argument("--class-mass-floor", type=float, default=0.1)
    parser.add_argument("--poles", default="0.1,0.2,0.3,0.4,0.5")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    g = args.group
    if not 0 <= args.class_weight <= g or args.class_mass_floor <= 0:
        raise SystemExit("pure class witness: invalid class or floor")
    profile = [0] * (g + 1)
    profile[args.class_weight] = atom_count(g)
    fugacities = np.asarray(
        [args.class_mass_floor / math.comb(g, weight) for weight in range(g + 1)]
    )
    fugacities[args.class_weight] = 1.0 / math.comb(g, args.class_weight)
    split_caps = split_cap_table()
    candidates = []
    for pole in (float(value) for value in args.poles.split(",")):
        probability, details = inner_probability(
            g, profile, pole, fugacities, split_caps, iterations=40
        )
        candidates.append((probability, pole, details))
    inner, pole, inner_details = min(candidates, key=lambda row: row[0])
    outer_log, outer_result, outer_components = optimize_spectrum_outer(g, profile)
    outer = outer_log / math.log(2.0)
    report = {
        "status": "DIAGNOSTIC_PURE_CLASS_MASS_FLOOR_SHARED_DRIVE",
        "group_bits": g,
        "class_weight": args.class_weight,
        "class_mass_floor": args.class_mass_floor,
        "profile": profile,
        "pole": pole,
        "fugacities": fugacities.tolist(),
        "inner_probability_log2": inner,
        "inner_details": inner_details,
        "spectrum_outer_log2": outer,
        "outer_point": outer_result.x.tolist(),
        "outer_components": [outer_components[1], outer_components[2]],
        "combined_before_graph_log2": outer + inner,
    }
    rendered = json.dumps(report, indent=2)
    print(rendered)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
