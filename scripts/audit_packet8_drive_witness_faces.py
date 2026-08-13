#!/usr/bin/env python3
"""Audit the frozen hard-face drive witness on every nested triangle.

This is an exhaustive finite audit of the 198580 nested-prefix triangle
barycenters and vertices.  It combines the existing atlas with the one new
shared-drive witness; it does not tune new parameters.  Binary64 is used for
the mass audit, while the hard anchor itself is separately closed by the
outward composition certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from probe_packet8_adaptive_simplex_atlas import load_witness_cache, split_cap_table
from probe_packet8_hard_face_drive_witness import DEFAULT_ANCHORS, DEFAULT_NAME, load_anchor
from probe_packet8_nested_face_landscape import nested_triples
from probe_packet8_profile_edge_interval_cover import orbit_value
from probe_packet8_profile_ordered_chambers import add_full_bijection, uniform_subset
from probe_packet8_profile_simplex_landscape import (
    ANCHORS,
    Anchor,
    PROFILE_COUNT_LOG2,
    build_fixed_witness,
    integral_profile,
    load_anchor_file,
    normalization_log2,
)
from probe_packet8_three_band_profile_enumerator import D, GROUPS, K
from certify_packet8_hard_face_drive_composition import (
    GRAPH_REPLACEMENTS,
    OUTER_PACKET_VARIABLES,
    outer_s2,
)


DEFAULT_ATLAS = Path("out/packet8_adaptive_atlas_wide_probe.json")
DEFAULT_CERTIFICATE = Path("out/packet8_hard_face_drive_inner_certificate.json")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--inner-certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--name", default=DEFAULT_NAME)
    parser.add_argument("--show-worst", type=int, default=20)
    parser.add_argument(
        "--extra-shared-report", type=Path, action="append", default=[]
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    anchor = load_anchor(args.atlas, args.name)
    certificate = json.loads(args.inner_certificate.read_text(encoding="utf-8"))
    inner_mgf_upper = float(certificate["inner_mgf_log2_interval"][1])
    target = -40.0 - PROFILE_COUNT_LOG2

    anchors = ANCHORS + load_anchor_file(args.atlas)
    witnesses = load_witness_cache(
        args.atlas.with_suffix(".witnesses.npz"), anchors
    )
    if witnesses is None:
        raise SystemExit("drive face audit: matching witness cache missing")
    witnesses = add_full_bijection(witnesses)
    old_constants = np.asarray([row.constant_log2 for row in witnesses])
    old_charges = np.vstack([row.linear_charge for row in witnesses])
    if args.extra_shared_report:
        split_caps = split_cap_table()
        extra_constants = []
        extra_charges = []
        for path in args.extra_shared_report:
            report = json.loads(path.read_text(encoding="utf-8"))
            source = report["anchor"]
            residual_anchor = Anchor(
                str(source["name"]),
                tuple(int(value) for value in report["profile"]),
                float(source["pole"]),
                tuple(float(value) for value in source["fugacities"]),
                float(source["outer_log_bound"]),
            )
            fixed = build_fixed_witness(residual_anchor, split_caps)
            extra_constants.append(
                fixed.constant_log2
                - float(report["shared_drive_improvement_log2"])
            )
            extra_charges.append(fixed.linear_charge)
        old_constants = np.concatenate(
            (old_constants, np.asarray(extra_constants, dtype=np.float64))
        )
        old_charges = np.vstack((old_charges, np.vstack(extra_charges)))

    outer_base = K + (GROUPS / 2) * math.log2(float(outer_s2()))
    adjacent_ratio = max(
        OUTER_PACKET_VARIABLES[new] / OUTER_PACKET_VARIABLES[old]
        for old in range(9)
        for new in range(9)
        if abs(new - old) <= 1
    )
    graph_correction = GRAPH_REPLACEMENTS * math.log2(adjacent_ratio)
    new_constant = (
        outer_base
        + graph_correction
        + inner_mgf_upper
        - D * math.log2(anchor.pole)
    )
    new_charge = np.log2(np.asarray(OUTER_PACKET_VARIABLES, dtype=np.float64))
    new_charge += np.log2(np.asarray(anchor.fugacities, dtype=np.float64))

    triples = nested_triples()
    points = {subset: uniform_subset(subset) for subset in range(2, 1 << 9)}
    profiles = np.vstack(
        [
            (points[left] + points[middle] + points[right]) / 3.0
            for left, middle, right in triples
        ]
    )
    normalizations = normalization_log2(profiles)
    old_best = np.full(len(profiles), np.inf)
    for start in range(0, len(old_constants), 128):
        stop = min(start + 128, len(old_constants))
        old_best = np.minimum(
            old_best,
            np.min(
                old_constants[None, start:stop]
                - profiles @ old_charges[start:stop].T
                - normalizations[:, None],
                axis=1,
            ),
        )
    orbit = np.asarray(
        [orbit_value(0.0, point, np.zeros(9), 0.0) for point in profiles]
    )
    old_best = np.minimum(old_best, orbit)
    new_values = new_constant - profiles @ new_charge - normalizations
    combined = np.minimum(old_best, new_values)

    new_vertex_safe = {
        subset: (
            new_constant
            - float(points[subset] @ new_charge)
            - float(normalization_log2(points[subset][None, :])[0])
            <= target
        )
        for subset in points
    }
    whole_face_new = np.asarray(
        [
            new_vertex_safe[left]
            and new_vertex_safe[middle]
            and new_vertex_safe[right]
            for left, middle, right in triples
        ],
        dtype=bool,
    )

    worst_indices = np.argsort(combined)[-args.show_worst :][::-1]
    worst = []
    for index in worst_indices:
        left, middle, right = triples[int(index)]
        worst.append(
            {
                "face": f"{left:03x}<{middle:03x}<{right:03x}",
                "profile": integral_profile(profiles[int(index)]).astype(int).tolist(),
                "old_best_log2": float(old_best[index]),
                "drive_witness_log2": float(new_values[index]),
                "combined_log2": float(combined[index]),
            }
        )
    report = {
        "status": "DIAGNOSTIC_BINARY64_EXHAUSTIVE_DRIVE_FACE_AUDIT",
        "target_log2": target,
        "triangles": len(triples),
        "existing_atlas_uncovered_barycenters": int(np.sum(old_best > target)),
        "extra_shared_reports": len(args.extra_shared_report),
        "drive_witness_covered_existing_residual": int(
            np.sum((old_best > target) & (new_values <= target))
        ),
        "combined_uncovered_barycenters": int(np.sum(combined > target)),
        "whole_triangles_closed_by_drive_witness_convexity": int(np.sum(whole_face_new)),
        "drive_witness_safe_vertices": int(sum(new_vertex_safe.values())),
        "total_vertices": len(new_vertex_safe),
        "worst_combined_barycenters": worst,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
