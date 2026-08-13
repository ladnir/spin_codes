#!/usr/bin/env python3
"""Compose the shared-drive inner transfer with one fixed outer witness.

The old and tightened inner operators use the same pole and profile
fugacities.  Therefore their difference is exactly the difference of their
Perron MGF bounds; the fixed outer Cauchy bound, profile charge, and orbit
normalization cancel.  This script first rebuilds the old fixed witness, then
subtracts that inner improvement.  All arithmetic is binary64 and the result
is a proof-design diagnostic until an outward certificate is supplied.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import EXACT_SLICES, SPECTRUM
from certificate_spectra import load_ebch128_spectrum
from probe_packet8_drive_stratified_caps import (
    block_histograms as drive_histograms,
    exact_point_caps,
)
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_profile_simplex_landscape import (
    M,
    PROFILE_COUNT_LOG2,
    Anchor,
    build_fixed_witness,
    normalization_log2,
)
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_weight_transfer import build_split_caps, load_exact
from probe_packet8_weight_profile_state_transfer import (
    D,
    INNER_BLOCKS,
    RobustProfileKernel,
)


DEFAULT_ANCHORS = Path("out/packet8_adaptive_atlas_wide_probe.json")
DEFAULT_NAME = "face_100_1e0_1ee_494_sharp"


def load_anchor(path: Path, name: str) -> Anchor:
    with path.open(encoding="utf-8") as handle:
        rows = json.load(handle)
    for row in rows:
        if row["name"] == name:
            return Anchor(
                str(row["name"]),
                tuple(int(value) for value in row["profile"]),
                float(row["pole"]),
                tuple(float(value) for value in row["fugacities"]),
                float(row["outer_log_bound"]),
            )
    raise SystemExit(f"hard-face drive witness: anchor {name!r} not found in {path}")


def inner_mgf_log2(kernel, iterations: int) -> tuple[float, dict[str, float | int]]:
    eigenvalue, domination, values, worst_state = witness(kernel, iterations)
    value = (
        math.log2(domination)
        + INNER_BLOCKS * math.log2(eigenvalue)
        + math.log2(float(values[0]))
    )
    return value, {
        "lambda_log2": math.log2(eigenvalue),
        "domination_log2": math.log2(domination),
        "initial_witness_log2": math.log2(float(values[0])),
        "worst_state": int(worst_state),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchors", type=Path, default=DEFAULT_ANCHORS)
    parser.add_argument("--name", default=DEFAULT_NAME)
    parser.add_argument("--iterations", type=int, default=64)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    anchor = load_anchor(args.anchors, args.name)
    spectrum = {
        weight: count
        for weight, count in load_ebch128_spectrum(SPECTRUM)
        if weight and count
    }
    split_caps = build_split_caps(spectrum, load_exact(EXACT_SLICES))

    # Rebuild the complete old fixed-witness value, including the identical
    # outer bound and exact profile normalization used by the atlas.
    old_fixed = build_fixed_witness(anchor, split_caps)
    profile = np.asarray(anchor.profile, dtype=np.float64)
    normalization = float(normalization_log2(profile[None, :])[0])
    old_combined = (
        old_fixed.constant_log2
        - float(profile @ old_fixed.linear_charge)
        - normalization
    )

    fugacities = np.asarray(anchor.fugacities, dtype=np.float64)
    histograms = drive_histograms(fugacities)
    old_kernel = RobustProfileKernel(
        np.sum(histograms, axis=1),
        split_caps,
        anchor.pole,
        float(np.max(fugacities)) ** 8,
    )
    new_kernel = SharedDriveStratifiedKernel(
        histograms,
        exact_point_caps(fugacities),
        split_caps,
        anchor.pole,
    )
    old_inner, old_details = inner_mgf_log2(old_kernel, args.iterations)
    new_inner, new_details = inner_mgf_log2(new_kernel, args.iterations)
    improvement = old_inner - new_inner
    tightened_combined = old_combined - improvement
    target = -40.0 - PROFILE_COUNT_LOG2

    report = {
        "status": "DIAGNOSTIC_FIXED_OUTER_SHARED_DRIVE_BINARY64",
        "anchor": anchor.name,
        "profile": list(anchor.profile),
        "profile_sum": int(sum(anchor.profile)),
        "expected_profile_sum": M,
        "pole": anchor.pole,
        "fugacities": list(anchor.fugacities),
        "iterations": args.iterations,
        "old_combined_log2": old_combined,
        "old_inner_mgf_log2": old_inner,
        "tightened_inner_mgf_log2": new_inner,
        "inner_improvement_log2": improvement,
        "inner_improvement_per_group": improvement / INNER_BLOCKS,
        "tightened_combined_log2": tightened_combined,
        "per_profile_target_log2": target,
        "certification_margin_log2": target - tightened_combined,
        "old_inner": old_details,
        "tightened_inner": new_details,
        "distance_cutoff": D,
        "inner_groups": INNER_BLOCKS,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
