#!/usr/bin/env python3
"""First optimized cross-g curve on the canonical packet-8 hard distribution."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from packet_group_drive_stratified import SUPPORTED_GROUPS, profile_count
from packet_group_outer_profile import N, atom_count
from packet_group_profile_bound import (
    outer_probability,
    split_cap_table,
    tune_inner_density,
)
from probe_packet8_hard_face_drive_witness import DEFAULT_ANCHORS, DEFAULT_NAME, load_anchor


def largest_remainder(probabilities: np.ndarray, total: int) -> list[int]:
    raw = probabilities * total
    result = np.floor(raw).astype(np.int64)
    order = np.argsort(raw - result)[::-1]
    result[order[: total - int(np.sum(result))]] += 1
    return result.tolist()


def mapped_distribution(base: np.ndarray, group_bits: int) -> np.ndarray:
    if group_bits == 8:
        return base.copy()
    if group_bits < 8:
        result = np.zeros(group_bits + 1)
        for old_weight, mass in enumerate(base):
            for new_weight in range(group_bits + 1):
                if new_weight <= old_weight and old_weight - new_weight <= 8 - group_bits:
                    result[new_weight] += mass * (
                        math.comb(group_bits, new_weight)
                        * math.comb(8 - group_bits, old_weight - new_weight)
                        / math.comb(8, old_weight)
                    )
        return result / np.sum(result)
    copies = group_bits // 8
    result = np.array([1.0])
    for _ in range(copies):
        result = np.convolve(result, base)
    return result / np.sum(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchors", type=Path, default=DEFAULT_ANCHORS)
    parser.add_argument("--name", default=DEFAULT_NAME)
    parser.add_argument("--groups", default="1,2,4,8,16,32,64")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    groups = [int(value) for value in args.groups.split(",")]
    if any(group not in SUPPORTED_GROUPS for group in groups):
        raise SystemExit("mapped hard curve: unsupported group")
    anchor = load_anchor(args.anchors, args.name)
    base = np.asarray(anchor.profile, dtype=np.float64) / sum(anchor.profile)
    split_caps = split_cap_table()
    rows = []
    for group_bits in groups:
        distribution = mapped_distribution(base, group_bits)
        profile = largest_remainder(distribution, atom_count(group_bits))
        inner = tune_inner_density(
            group_bits,
            profile,
            split_caps,
            poles=(0.05, 0.1, 0.2, 0.3),
            exponents=(0.5, 0.75, 1.0),
            iterations=40,
        )
        outer, outer_details = outer_probability(group_bits, profile)
        combined = outer + inner[0]
        union_target = -40.0 - math.log2(profile_count(group_bits, N))
        row = {
            "group_bits": group_bits,
            "profile": profile,
            "profile_proportions": distribution.tolist(),
            "pole": inner[1],
            "density_exponent": inner[2],
            "fugacities": inner[3].tolist(),
            "inner_probability_log2": inner[0],
            "outer_profile_log2": outer,
            "combined_log2": combined,
            "union_target_log2": union_target,
            "diagnostic_margin_log2": union_target - combined,
            "inner_details": inner[4],
            "outer_details": outer_details,
        }
        rows.append(row)
        print(
            f"g={group_bits} combined={combined:.6f} target={union_target:.6f} "
            f"margin={union_target-combined:.6f} pole={inner[1]} "
            f"density_exponent={inner[2]}",
            flush=True,
        )
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    print("status=DIAGNOSTIC_MAPPED_HARD_FAMILY_OPTIMIZED_DENSITY_GRID")


if __name__ == "__main__":
    main()
