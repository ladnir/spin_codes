#!/usr/bin/env python3
"""Greedy shared-drive atlas over every nested packet-profile barycenter."""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path

import numpy as np

from audit_packet8_drive_witness_faces import (
    DEFAULT_ATLAS,
    DEFAULT_CERTIFICATE,
)
from certify_packet8_hard_face_drive_composition import (
    GRAPH_REPLACEMENTS,
    OUTER_PACKET_VARIABLES,
    outer_s2,
)
from probe_packet8_adaptive_simplex_atlas import (
    load_witness_cache,
    split_cap_table,
)
from probe_packet8_drive_residual_anchor import retune_shared_report
from probe_packet8_hard_face_drive_witness import (
    DEFAULT_NAME,
    load_anchor,
)
from probe_packet8_nested_face_landscape import nested_triples
from probe_packet8_profile_edge_interval_cover import orbit_value
from probe_packet8_profile_ordered_chambers import (
    add_full_bijection,
    uniform_subset,
)
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


def fixed_from_report(report: dict, split_caps):
    source = report["anchor"]
    anchor = Anchor(
        str(source["name"]),
        tuple(int(value) for value in report["profile"]),
        float(source["pole"]),
        tuple(float(value) for value in source["fugacities"]),
        float(source["outer_log_bound"]),
    )
    fixed = build_fixed_witness(anchor, split_caps)
    return (
        fixed.constant_log2
        - float(report["shared_drive_improvement_log2"]),
        fixed.linear_charge,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument(
        "--inner-certificate", type=Path, default=DEFAULT_CERTIFICATE
    )
    parser.add_argument("--hard-name", default=DEFAULT_NAME)
    parser.add_argument(
        "--existing-report", type=Path, action="append", default=[]
    )
    parser.add_argument("--upgraded-cache", type=Path)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--coordinate-iterations", type=int, default=4)
    parser.add_argument("--poles", default="0.1,0.2,0.3")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    if args.iterations < 0 or args.coordinate_iterations <= 0:
        raise SystemExit("adaptive shared atlas: invalid iteration limit")
    poles = [float(value) for value in args.poles.split(",")]

    hard_anchor = load_anchor(args.atlas, args.hard_name)
    certificate = json.loads(
        args.inner_certificate.read_text(encoding="utf-8")
    )
    inner_mgf_upper = float(certificate["inner_mgf_log2_interval"][1])
    target = -40.0 - PROFILE_COUNT_LOG2

    anchors = ANCHORS + load_anchor_file(args.atlas)
    witnesses = load_witness_cache(
        args.atlas.with_suffix(".witnesses.npz"), anchors
    )
    if witnesses is None:
        raise SystemExit("adaptive shared atlas: matching base cache missing")
    if args.upgraded_cache:
        with np.load(args.upgraded_cache, allow_pickle=False) as cache:
            constants = [float(value) for value in cache["constants"]]
            charges = [np.asarray(row) for row in cache["charges"]]
        # The upgraded cache contains the paired atlas witnesses but not the
        # separate full-bijection branch.
        full = add_full_bijection([])[0]
        constants.append(float(full.constant_log2))
        charges.append(np.asarray(full.linear_charge))
    else:
        witnesses = add_full_bijection(witnesses)
        constants = [float(row.constant_log2) for row in witnesses]
        charges = [np.asarray(row.linear_charge) for row in witnesses]

    split_caps = split_cap_table()
    reports = []
    for path in args.existing_report:
        report = json.loads(path.read_text(encoding="utf-8"))
        constant, charge = fixed_from_report(report, split_caps)
        constants.append(constant)
        charges.append(charge)
        reports.append({"path": str(path), "report": report})

    outer_base = K + (GROUPS / 2) * math.log2(float(outer_s2()))
    adjacent_ratio = max(
        OUTER_PACKET_VARIABLES[new] / OUTER_PACKET_VARIABLES[old]
        for old in range(9)
        for new in range(9)
        if abs(new - old) <= 1
    )
    hard_constant = (
        outer_base
        + GRAPH_REPLACEMENTS * math.log2(adjacent_ratio)
        + inner_mgf_upper
        - D * math.log2(hard_anchor.pole)
    )
    hard_charge = np.log2(
        np.asarray(OUTER_PACKET_VARIABLES, dtype=np.float64)
    ) + np.log2(np.asarray(hard_anchor.fugacities, dtype=np.float64))
    constants.append(hard_constant)
    charges.append(hard_charge)

    triples = nested_triples()
    points = {subset: uniform_subset(subset) for subset in range(2, 1 << 9)}
    profiles = np.vstack(
        [
            (points[left] + points[middle] + points[right]) / 3.0
            for left, middle, right in triples
        ]
    )
    normalizations = normalization_log2(profiles)
    values = np.full(len(profiles), np.inf)
    for start in range(0, len(constants), 128):
        stop = min(start + 128, len(constants))
        batch_constants = np.asarray(constants[start:stop])
        batch_charges = np.vstack(charges[start:stop])
        values = np.minimum(
            values,
            np.min(
                batch_constants[None, :]
                - profiles @ batch_charges.T
                - normalizations[:, None],
                axis=1,
            ),
        )
    orbit = np.asarray(
        [orbit_value(0.0, point, np.zeros(9), 0.0) for point in profiles]
    )
    values = np.minimum(values, orbit)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"profiles={len(profiles)} initial_uncovered={int(np.sum(values > target))} "
        f"initial_worst={float(np.max(values)):.9f} target={target:.9f}",
        flush=True,
    )
    generated = []
    for iteration in range(args.iterations):
        worst_index = int(np.argmax(values))
        left, middle, right = triples[worst_index]
        face = f"{left:03x}_{middle:03x}_{right:03x}"
        profile = integral_profile(profiles[worst_index]).astype(int).tolist()
        before = float(values[worst_index])
        report = retune_shared_report(
            profile,
            f"adaptive_shared_{iteration:03d}_{face}",
            poles=poles,
            coordinate_iterations=args.coordinate_iterations,
        )
        path = args.output_dir / f"adaptive_shared_{iteration:03d}_{face}.json"
        path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        constant, charge = fixed_from_report(report, split_caps)
        candidate = constant - profiles @ charge - normalizations
        improved = int(np.sum(candidate < values))
        values = np.minimum(values, candidate)
        uncovered = int(np.sum(values > target))
        generated.append(
            {
                "path": str(path),
                "face": face,
                "selected_before_log2": before,
                "anchor_combined_log2": float(
                    report["shared_drive_combined_log2"]
                ),
                "improved_barycenters": improved,
                "uncovered_after": uncovered,
                "worst_after_log2": float(np.max(values)),
            }
        )
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(
            json.dumps(
                {
                    "status": "DIAGNOSTIC_ADAPTIVE_SHARED_FACE_ATLAS",
                    "target_log2": target,
                    "existing_reports": [str(path) for path in args.existing_report],
                    "generated": generated,
                    "uncovered": uncovered,
                    "worst_log2": float(np.max(values)),
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        print(
            f"iteration={iteration} face={face} before={before:.9f} "
            f"anchor={report['shared_drive_combined_log2']:.9f} "
            f"improved={improved} uncovered={uncovered} "
            f"worst={float(np.max(values)):.9f}",
            flush=True,
        )
        if uncovered == 0:
            break
    print(
        f"final_uncovered={int(np.sum(values > target))} "
        f"final_worst={float(np.max(values)):.9f}",
        flush=True,
    )
    uncovered_indices = np.flatnonzero(values > target)
    support_counts = Counter(
        sum(
            1 << weight
            for weight, count in enumerate(profiles[index])
            if count > 1e-9
        )
        for index in uncovered_indices
    )
    worst_indices = uncovered_indices[
        np.argsort(values[uncovered_indices])[-20:][::-1]
    ] if len(uncovered_indices) else np.asarray([], dtype=np.int64)
    final_report = {
        "status": "DIAGNOSTIC_ADAPTIVE_SHARED_FACE_ATLAS",
        "target_log2": target,
        "existing_reports": [str(path) for path in args.existing_report],
        "generated": generated,
        "uncovered": int(len(uncovered_indices)),
        "worst_log2": float(np.max(values)),
        "uncovered_support_masks": [
            {"mask": f"{mask:03x}", "count": count}
            for mask, count in support_counts.most_common()
        ],
        "worst_residuals": [
            {
                "face": "%03x_%03x_%03x" % triples[int(index)],
                "value_log2": float(values[index]),
                "profile": integral_profile(profiles[index]).astype(int).tolist(),
            }
            for index in worst_indices
        ],
    }
    args.manifest.write_text(
        json.dumps(final_report, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
