#!/usr/bin/env python3
"""Replace every cached robust atlas inner constant by shared-drive constants."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from probe_packet8_adaptive_simplex_atlas import load_witness_cache, split_cap_table
from probe_packet8_drive_residual_anchor import mgf
from probe_packet8_drive_stratified_caps import block_histograms, exact_point_caps
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_profile_simplex_landscape import ANCHORS, load_anchor_file
from probe_packet8_weight_profile_state_transfer import RobustProfileKernel


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output-cache", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--iterations", type=int, default=64)
    args = parser.parse_args()

    anchors = ANCHORS + load_anchor_file(args.atlas)
    witnesses = load_witness_cache(
        args.atlas.with_suffix(".witnesses.npz"), anchors
    )
    if witnesses is None:
        raise SystemExit("shared atlas upgrade: matching robust cache missing")
    by_name = {row.name: row for row in witnesses}
    completed = {}
    if args.resume and args.checkpoint.exists():
        rows = json.loads(args.checkpoint.read_text(encoding="utf-8"))
        completed = {str(row["name"]): row for row in rows}

    split_caps = split_cap_table()
    selected = anchors if args.limit is None else anchors[: args.limit]
    for index, anchor in enumerate(selected):
        if anchor.name in completed:
            continue
        fugacities = np.asarray(anchor.fugacities, dtype=np.float64)
        histograms = block_histograms(fugacities)
        robust = RobustProfileKernel(
            np.sum(histograms, axis=1),
            split_caps,
            anchor.pole,
            float(np.max(fugacities)) ** 8,
        )
        shared = SharedDriveStratifiedKernel(
            histograms,
            exact_point_caps(fugacities),
            split_caps,
            anchor.pole,
        )
        robust_mgf = mgf(robust, args.iterations)
        shared_mgf = mgf(shared, args.iterations)
        improvement = robust_mgf - shared_mgf
        completed[anchor.name] = {
            "name": anchor.name,
            "pole": anchor.pole,
            "robust_mgf_log2": robust_mgf,
            "shared_mgf_log2": shared_mgf,
            "improvement_log2": improvement,
        }
        args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
        args.checkpoint.write_text(
            json.dumps(list(completed.values()), indent=2) + "\n",
            encoding="utf-8",
        )
        print(
            f"anchor={index + 1}/{len(selected)} name={anchor.name} "
            f"improvement={improvement:.9f}",
            flush=True,
        )

    names = []
    constants = []
    charges = []
    upgraded = 0
    for fixed in witnesses:
        base_name = (
            fixed.name[: -len("_spectrum")]
            if fixed.name.endswith("_spectrum")
            else fixed.name
        )
        improvement = max(
            0.0,
            float(completed.get(base_name, {}).get("improvement_log2", 0.0)),
        )
        upgraded += improvement != 0.0
        names.append(fixed.name)
        constants.append(float(fixed.constant_log2) - improvement)
        charges.append(np.asarray(fixed.linear_charge))
    args.output_cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_cache,
        names=np.asarray(names),
        constants=np.asarray(constants),
        charges=np.vstack(charges),
        upgraded_fixed_witnesses=np.asarray(upgraded),
        upgraded_anchors=np.asarray(len(completed)),
    )
    print(
        f"completed_anchors={len(completed)} "
        f"upgraded_fixed_witnesses={upgraded}",
        flush=True,
    )


if __name__ == "__main__":
    main()
