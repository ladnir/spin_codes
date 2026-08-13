#!/usr/bin/env python3
"""Retune one residual profile, then replace its robust inner by shared-drive."""

from __future__ import annotations

import argparse
import json
import math
from types import SimpleNamespace

import numpy as np

from probe_packet8_adaptive_simplex_atlas import make_anchor, split_cap_table
from probe_packet8_drive_stratified_caps import block_histograms, exact_point_caps
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_profile_simplex_landscape import build_fixed_witness, normalization_log2
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_weight_profile_scalar import parse_profile
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS, RobustProfileKernel


def mgf(kernel, iterations: int) -> float:
    eigenvalue, domination, values, _worst = witness(kernel, iterations)
    return (
        math.log2(domination)
        + INNER_BLOCKS * math.log2(eigenvalue)
        + math.log2(float(values[0]))
    )


def retune_shared_report(
    profile: list[int],
    name: str,
    *,
    poles: list[float],
    coordinate_iterations: int,
) -> dict:
    split_caps = split_cap_table()
    tuning = SimpleNamespace(
        poles=poles,
        tuning_passes=1,
        coordinate_iterations=coordinate_iterations,
        witness_iterations=32,
        outer_bounds=[3.0, 6.0, 10.0, 40.0],
        anchor_target=-50000.0,
    )
    anchor, robust_tuned_value, _outer, _inner, evaluations = make_anchor(
        profile, name, split_caps, tuning, 1e-30
    )
    fixed = build_fixed_witness(anchor, split_caps)
    vector = np.asarray(profile, dtype=np.float64)
    robust_combined = (
        fixed.constant_log2
        - float(vector @ fixed.linear_charge)
        - float(normalization_log2(vector[None, :])[0])
    )
    fugacities = np.asarray(anchor.fugacities)
    histograms = block_histograms(fugacities)
    robust = RobustProfileKernel(
        np.sum(histograms, axis=1), split_caps, anchor.pole,
        float(np.max(fugacities)) ** 8,
    )
    shared = SharedDriveStratifiedKernel(
        histograms, exact_point_caps(fugacities), split_caps, anchor.pole
    )
    improvement = mgf(robust, 64) - mgf(shared, 64)
    return {
        "status": "DIAGNOSTIC_RETUNED_RESIDUAL_SHARED_DRIVE_BINARY64",
        "profile": profile,
        "anchor": {
            "name": anchor.name,
            "pole": anchor.pole,
            "fugacities": list(anchor.fugacities),
            "outer_log_bound": anchor.outer_log_bound,
        },
        "tuning_evaluations": evaluations,
        "robust_tuned_value_from_make_anchor": robust_tuned_value,
        "robust_fixed_combined_log2": robust_combined,
        "shared_drive_improvement_log2": improvement,
        "shared_drive_combined_log2": robust_combined - improvement,
        "shared_drive_improvement_per_group": improvement / INNER_BLOCKS,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--name", default="residual")
    parser.add_argument("--poles", default="0.1,0.2,0.3")
    parser.add_argument("--coordinate-iterations", type=int, default=4)
    parser.add_argument("--output")
    args = parser.parse_args()
    profile = list(parse_profile(args.profile))
    report = retune_shared_report(
        profile,
        args.name,
        poles=[float(value) for value in args.poles.split(",")],
        coordinate_iterations=args.coordinate_iterations,
    )
    rendered = json.dumps(report, indent=2)
    print(rendered)
    if args.output:
        from pathlib import Path
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
