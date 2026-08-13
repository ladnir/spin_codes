#!/usr/bin/env python3
"""Locate the sparse-to-transfer handoff on a binary packet-value edge."""

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
from packet_group_outer_profile import D, N, atom_count, normalization_log2, total_weight_outer
from packet_group_profile_bound import outer_probability, split_cap_table
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS


def parse_counts(text: str) -> list[int]:
    result = []
    for part in text.split(","):
        if ":" in part:
            begin, end, step = (int(value) for value in part.split(":"))
            result.extend(range(begin, end + 1, step))
        else:
            result.append(int(part))
    return sorted(set(result))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, default=2)
    parser.add_argument("--base-class", type=int, default=0)
    parser.add_argument("--active-class", type=int, required=True)
    parser.add_argument(
        "--counts",
        default="11,16,21,32,48,64,96,128,192,256,384,512,768,1024,1536,2048,3072,4096,6144,8192,12288,16384,24576,32768,49152,65536,98304,131072,172032,262144,524288,1048576",
    )
    parser.add_argument("--ratio-offsets", default="-8,-6,-4,-2,0,2,4,6,8")
    parser.add_argument("--poles", default="0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8")
    parser.add_argument("--screen-iterations", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    g = args.group
    if not 0 <= args.base_class <= g or not 0 <= args.active_class <= g:
        raise SystemExit("binary edge: invalid class")
    if args.base_class == args.active_class:
        raise SystemExit("binary edge: base and active classes must differ")
    atoms = atom_count(g)
    classes = profile_classes(g)
    counts = [value for value in parse_counts(args.counts) if 0 < value <= atoms]
    offsets = [float(value) for value in args.ratio_offsets.split(",")]
    poles = [float(value) for value in args.poles.split(",")]
    split_caps = split_cap_table()
    target = -40.0 - math.log2(profile_count(g, N))
    rows = []
    for active in counts:
        profile = [0] * (g + 1)
        profile[args.base_class] = atoms - active
        profile[args.active_class] = active
        normalization = float(normalization_log2(g, np.asarray(profile))[0])
        empirical = (
            (active / classes[args.active_class])
            / max(1.0 / classes[args.base_class], (atoms - active) / classes[args.base_class])
        )
        candidates = []
        for offset in offsets:
            ratio = max(2.0**-128, empirical * 2.0**offset)
            fugacities = np.zeros(g + 1)
            fugacities[args.base_class] = 1.0
            fugacities[args.active_class] = ratio
            histograms = block_histograms(g, fugacities)
            caps = point_caps(g, fugacities)
            charge = active * math.log2(ratio)
            for pole in poles:
                eigenvalue, domination, values, worst_state = witness(
                    SharedDriveStratifiedKernel(
                        histograms, caps, split_caps, pole
                    ),
                    args.screen_iterations,
                )
                mgf = (
                    math.log2(domination)
                    + INNER_BLOCKS * math.log2(eigenvalue)
                    + math.log2(float(values[0]))
                )
                probability = mgf - charge - normalization - D * math.log2(pole)
                candidates.append(
                    (probability, ratio, pole, worst_state)
                )
        screened = min(candidates, key=lambda row: row[0])
        ratio, pole = screened[1], screened[2]
        fugacities = np.zeros(g + 1)
        fugacities[args.base_class] = 1.0
        fugacities[args.active_class] = ratio
        histograms = block_histograms(g, fugacities)
        caps = point_caps(g, fugacities)
        eigenvalue, domination, values, worst_state = witness(
            SharedDriveStratifiedKernel(histograms, caps, split_caps, pole), 40
        )
        mgf = (
            math.log2(domination)
            + INNER_BLOCKS * math.log2(eigenvalue)
            + math.log2(float(values[0]))
        )
        inner = (
            mgf
            - active * math.log2(ratio)
            - normalization
            - D * math.log2(pole)
        )
        physical_weight = (
            args.base_class * (atoms - active) + args.active_class * active
        )
        total_outer, outer_result, interval = total_weight_outer(physical_weight)
        profile_outer, profile_outer_details = outer_probability(
            g, profile, fast=True
        )
        outer = min(total_outer, profile_outer)
        outer_type = "total_weight" if total_outer <= profile_outer else "linear_bl"
        combined = outer + min(0.0, inner)
        row = {
            "active_atoms": active,
            "profile": profile,
            "physical_weight": physical_weight,
            "ratio": ratio,
            "pole": pole,
            "inner_probability_log2": inner,
            "outer_log2": outer,
            "outer_type": outer_type,
            "total_weight_outer_log2": total_outer,
            "linear_bl_outer_log2": profile_outer,
            "linear_bl_outer_details": profile_outer_details,
            "outer_log_pole": float(outer_result.x),
            "outer_interval": list(interval),
            "combined_log2": combined,
            "target_log2": target,
            "margin_bits": target - combined,
            "worst_state": int(worst_state),
        }
        rows.append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
        print(
            f"edge={args.base_class}->{args.active_class} active={active} inner={inner:.6f} "
            f"outer={outer:.6f} margin={row['margin_bits']:.6f}",
            flush=True,
        )
    print("status=DIAGNOSTIC_BINARY64_GENERALIZED_BINARY_EDGE_CURVE")


if __name__ == "__main__":
    main()
