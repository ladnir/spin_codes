#!/usr/bin/env python3
"""Optimize every nonzero pure atom-weight profile for one or more g."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from packet_group_drive_stratified import (
    block_histograms,
    point_caps,
    profile_classes,
    profile_count,
)
from packet_group_outer_profile import (
    D,
    N,
    atom_count,
    pure_profile_total_weight_outer,
)
from packet_group_profile_bound import split_cap_table
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", default="64")
    parser.add_argument("--poles", default="0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    groups = [int(value) for value in args.groups.split(",")]
    poles = [float(value) for value in args.poles.split(",")]
    split_caps = split_cap_table()
    reports = []
    for group_bits in groups:
        atoms = atom_count(group_bits)
        classes = profile_classes(group_bits)
        target = -40.0 - math.log2(profile_count(group_bits, N))
        rows = []
        for class_weight in range(1, group_bits + 1):
            fugacities = np.zeros(group_bits + 1)
            fugacities[class_weight] = 1.0
            histograms = block_histograms(group_bits, fugacities)
            caps = point_caps(group_bits, fugacities)
            normalization = atoms * math.log2(classes[class_weight])
            candidates = []
            for pole in poles:
                kernel = SharedDriveStratifiedKernel(
                    histograms, caps, split_caps, pole
                )
                eigenvalue, domination, values, worst_state = witness(kernel, 40)
                mgf = (
                    math.log2(domination)
                    + INNER_BLOCKS * math.log2(eigenvalue)
                    + math.log2(float(values[0]))
                )
                probability = mgf - normalization - D * math.log2(pole)
                candidates.append(
                    (probability, pole, mgf, math.log2(eigenvalue), worst_state)
                )
            inner = min(candidates, key=lambda row: row[0])
            outer, outer_result, interval = pure_profile_total_weight_outer(
                group_bits, class_weight
            )
            combined = outer + inner[0]
            row = {
                "class_weight": class_weight,
                "profile_count": atoms,
                "pole": inner[1],
                "inner_probability_log2": inner[0],
                "inner_mgf_log2": inner[2],
                "lambda_log2": inner[3],
                "worst_state": inner[4],
                "outer_log2": outer,
                "outer_log_pole": float(outer_result.x),
                "pre_replacement_weight_interval": list(interval),
                "combined_log2": combined,
                "target_log2": target,
                "margin_log2": target - combined,
            }
            rows.append(row)
            print(
                f"g={group_bits} class={class_weight} margin={row['margin_log2']:.6f} "
                f"pole={inner[1]:.6f}",
                flush=True,
            )
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(reports + [{"group_bits": group_bits, "rows": rows}], indent=2)
                + "\n",
                encoding="utf-8",
            )
        reports.append({"group_bits": group_bits, "rows": rows})
    args.output.write_text(json.dumps(reports, indent=2) + "\n", encoding="utf-8")
    print("status=DIAGNOSTIC_PURE_CLASS_SHARED_DRIVE_CURVE")


if __name__ == "__main__":
    main()
