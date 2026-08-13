#!/usr/bin/env python3
"""Generate g=4 witness-dominance transition anchors from hard mesh cells.

For each selected cell, the script finds the best finite frozen witness at
every vertex.  If an edge has different endpoint winners, their affine parts
are equated along that edge; the shared multinomial normalization cancels.
Interior crossings are rounded deterministically to legal exact-support
integer profiles and emitted in the residual-ledger format consumed by
``probe_packet_group_fixed_atlas_parallel.py``.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import N
from probe_packet_group_g4_stellar_cover import Evaluator, Witness, load_witnesses


GROUP_BITS = 4
MINIMUM_PHYSICAL_WEIGHT = 21

Profile = tuple[Fraction, ...]


def parse_fraction(value: Any) -> Fraction:
    if isinstance(value, bool) or isinstance(value, float):
        raise ValueError("transition anchors require exact non-binary64 profiles")
    if isinstance(value, int):
        return Fraction(value)
    return Fraction(str(value))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def select_stratum(payload: dict[str, Any], mask: int) -> dict[str, Any]:
    if "support_strata" in payload:
        matches = [
            row for row in payload["support_strata"] if int(row["support_mask"]) == mask
        ]
        if len(matches) != 1:
            raise ValueError(f"transition anchors: support mask {mask} is not unique")
        return matches[0]
    return payload


def load_geometry(
    payload: dict[str, Any], mask: int
) -> tuple[list[Profile], list[dict[str, Any]], tuple[int, ...]]:
    stratum = select_stratum(payload, mask)
    support = tuple(
        int(value)
        for value in stratum.get(
            "support", [index for index in range(5) if mask & (1 << index)]
        )
    )
    raw_anchors = stratum.get("anchor_profiles", payload.get("anchor_profiles"))
    if not isinstance(raw_anchors, list) or not raw_anchors:
        raise ValueError("transition anchors: artifact lacks anchor_profiles")
    anchors = [tuple(parse_fraction(value) for value in row) for row in raw_anchors]
    raw_cells = stratum.get("covered_cell_ledger")
    if not raw_cells:
        raw_cells = stratum.get("root_cell_ledger", payload.get("root_cell_ledger"))
    if not isinstance(raw_cells, list) or not raw_cells:
        raise ValueError("transition anchors: artifact lacks a usable cell ledger")
    return anchors, list(raw_cells), support


def cell_contribution(row: dict[str, Any]) -> float:
    for key in (
        "cell_local_contribution_log2_upperish",
        "contribution_log2_upperish",
        "support_contribution_log2",
    ):
        if key in row:
            return float(row[key])
    score = float(row.get("maximum_vertex_log2", row.get("best_max_vertex_value_log2", 0.0)))
    count = int(row.get("integer_profile_count_upper", 1))
    return score + (math.log2(count) if count else -math.inf)


def finite_winner(values: tuple[float, ...]) -> int:
    finite = [index for index, value in enumerate(values) if math.isfinite(value)]
    if not finite:
        raise ValueError("transition anchors: cell vertex has no finite witness")
    return min(finite, key=lambda index: (values[index], index))


def affine_difference(witnesses: list[Witness], left: int, right: int, profile: Profile) -> float:
    a = witnesses[left]
    b = witnesses[right]
    return (a.constant_log2 - b.constant_log2) - math.fsum(
        float(count) * (left_charge - right_charge)
        for count, left_charge, right_charge in zip(profile, a.charge, b.charge)
    )


def affine_value_exact(witness: Witness, profile: Profile) -> Fraction:
    """Evaluate the normalization-free affine part, treating binary64 as dyadic."""
    return Fraction.from_float(witness.constant_log2) - sum(
        (count * Fraction.from_float(coefficient)
         for count, coefficient in zip(profile, witness.charge)),
        Fraction(),
    )


def lower_envelope_transitions(
    witnesses: list[Witness], left: Profile, right: Profile
) -> list[tuple[Fraction, int, int]]:
    """Return every strict interior transition of the full affine lower envelope.

    Lines are sorted by decreasing slope and scanned into a lower hull.  This is
    O(w log w), rather than testing O(w^2) witness pairs on every mesh edge.
    """
    lines_by_slope: dict[Fraction, tuple[Fraction, int]] = {}
    for index, witness in enumerate(witnesses):
        if any(
            (left[class_index] or right[class_index])
            and class_index not in witness.support
            for class_index in range(5)
        ):
            continue
        intercept = affine_value_exact(witness, left)
        slope = affine_value_exact(witness, right) - intercept
        current = lines_by_slope.get(slope)
        if current is None or (intercept, index) < current:
            lines_by_slope[slope] = (intercept, index)
    ordered = sorted(
        ((slope, intercept, index) for slope, (intercept, index) in lines_by_slope.items()),
        key=lambda row: (-row[0], row[1], row[2]),
    )
    # Entries are (slope, intercept, witness index, first t where this line wins).
    hull: list[tuple[Fraction, Fraction, int, Fraction | None]] = []
    for slope, intercept, index in ordered:
        start: Fraction | None = None
        while hull:
            old_slope, old_intercept, _old_index, old_start = hull[-1]
            start = (old_intercept - intercept) / (slope - old_slope)
            if old_start is None or start > old_start:
                break
            hull.pop()
        if not hull:
            start = None
        hull.append((slope, intercept, index, start))
    active: list[tuple[int, Fraction, Fraction]] = []
    for position, (_slope, _intercept, index, start) in enumerate(hull):
        lo = max(Fraction(), start if start is not None else Fraction())
        next_start = hull[position + 1][3] if position + 1 < len(hull) else None
        hi = min(Fraction(1), next_start if next_start is not None else Fraction(1))
        if lo < hi:
            active.append((index, lo, hi))
    transitions = []
    for before, after in zip(active, active[1:]):
        parameter = before[2]
        if parameter == after[1] and 0 < parameter < 1 and before[0] != after[0]:
            transitions.append((parameter, before[0], after[0]))
    return transitions


def round_exact_support(
    point: Profile,
    support: tuple[int, ...],
    total_mass: int,
) -> tuple[int, ...]:
    profile = [0] * 5
    floors = {}
    for index in support:
        floor = point[index].numerator // point[index].denominator
        floors[index] = max(1, floor)
        profile[index] = floors[index]
    difference = total_mass - sum(profile)
    if difference > 0:
        order = sorted(
            support,
            key=lambda index: (point[index] - floors[index], -index),
            reverse=True,
        )
        for offset in range(difference):
            profile[order[offset % len(order)]] += 1
    elif difference < 0:
        order = sorted(
            support,
            key=lambda index: (point[index] - floors[index], index),
        )
        for _ in range(-difference):
            candidate = next(index for index in order if profile[index] > 1)
            profile[candidate] -= 1
    if sum(profile) != total_mass or any((profile[index] > 0) != (index in support) for index in range(5)):
        raise RuntimeError("transition anchor rounding lost exact support or mass")
    physical = sum(index * profile[index] for index in range(5))
    if physical < MINIMUM_PHYSICAL_WEIGHT:
        low = min((index for index in support if profile[index] > 1), default=None)
        high = max(support)
        if low is None or low == high:
            raise RuntimeError("transition anchor cannot satisfy physical-weight floor")
        needed = math.ceil((MINIMUM_PHYSICAL_WEIGHT - physical) / (high - low))
        move = min(needed, profile[low] - 1)
        profile[low] -= move
        profile[high] += move
    if sum(index * profile[index] for index in range(5)) < MINIMUM_PHYSICAL_WEIGHT:
        raise RuntimeError("transition anchor rounding missed physical-weight floor")
    return tuple(profile)


def generate(
    payload: dict[str, Any],
    witnesses: list[Witness],
    mask: int,
    max_cells: int,
    max_profiles: int,
) -> dict[str, Any]:
    normalization_modes = {witness.subtract_normalization for witness in witnesses}
    if len(normalization_modes) != 1:
        raise ValueError(
            "transition anchors: all witnesses must have identical normalization behavior"
        )
    anchors, cells, support = load_geometry(payload, mask)
    total_mass = int(sum(anchors[0]))
    evaluator = Evaluator(witnesses)
    selected = sorted(
        cells,
        key=lambda row: (-cell_contribution(row), str(row.get("cell_id", ""))),
    )[:max_cells]
    unique: dict[tuple[int, ...], dict[str, Any]] = {}
    audit_cells = []
    for row in selected:
        indices = tuple(int(value) for value in row["anchor_indices"])
        profiles = [anchors[index] for index in indices]
        values = [evaluator.evaluate(profile)[1] for profile in profiles]
        winners = [finite_winner(vertex_values) for vertex_values in values]
        contribution = cell_contribution(row)
        crossings = 0
        vanished_after_rounding = 0
        for left_local, right_local in itertools.combinations(range(len(indices)), 2):
            transitions = lower_envelope_transitions(
                witnesses, profiles[left_local], profiles[right_local]
            )
            for exact_parameter, left_winner, right_winner in transitions:
                point = tuple(
                    left + exact_parameter * (right - left)
                    for left, right in zip(profiles[left_local], profiles[right_local])
                )
                profile = round_exact_support(point, support, total_mass)
                _normalization, rounded_values = evaluator.evaluate(tuple(map(Fraction, profile)))
                finite_order = sorted(
                    (value, index)
                    for index, value in enumerate(rounded_values)
                    if math.isfinite(value)
                )
                transition_pair = {left_winner, right_winner}
                # Rounding may move the profile under an unrelated lower-envelope
                # branch.  Such a point no longer diagnoses this transition.
                if len(finite_order) < 2 or {
                    finite_order[0][1], finite_order[1][1]
                } != transition_pair:
                    vanished_after_rounding += 1
                    continue
                key = tuple(profile)
                record = {
                    "name": f"transition_{row['cell_id']}_{left_local}_{right_local}",
                    "profile": list(profile),
                    "value_log2": contribution,
                    "target_log2": -40.0 - math.log2(profile_count(GROUP_BITS, N)),
                    "gap_bits": contribution + 40.0 + math.log2(profile_count(GROUP_BITS, N)),
                    "parent_cell_id": str(row["cell_id"]),
                    "parent_contribution_log2": contribution,
                    "edge_local_vertices": [left_local, right_local],
                    "edge_anchor_indices": [indices[left_local], indices[right_local]],
                    "transition_witness_references": [
                        witnesses[left_winner].reference,
                        witnesses[right_winner].reference,
                    ],
                    "rounded_profile_best_two_references": [
                        witnesses[finite_order[0][1]].reference,
                        witnesses[finite_order[1][1]].reference,
                    ],
                    "rounded_profile_best_two_gap_log2": (
                        finite_order[1][0] - finite_order[0][0]
                    ),
                    "crossing_parameter_binary64": float(exact_parameter),
                    "crossing_parameter_exact_dyadic": str(exact_parameter),
                    "full_lower_envelope_verified_after_rounding": True,
                }
                current = unique.get(key)
                if current is None or contribution > current["parent_contribution_log2"]:
                    unique[key] = record
                crossings += 1
        audit_cells.append(
            {
                "cell_id": str(row["cell_id"]),
                "contribution_log2": contribution,
                "vertex_winner_references": [witnesses[index].reference for index in winners],
                "crossing_edges": crossings,
                "transitions_vanished_after_integer_rounding": vanished_after_rounding,
            }
        )
    residuals = sorted(
        unique.values(),
        key=lambda row: (-row["parent_contribution_log2"], tuple(row["profile"])),
    )[:max_profiles]
    return {
        "status": "DIAGNOSTIC_G4_WITNESS_DOMINANCE_TRANSITION_ANCHORS",
        "group_bits": GROUP_BITS,
        "support_mask": mask,
        "support": list(support),
        "selected_cells": len(selected),
        "unique_transition_profiles_before_limit": len(unique),
        "emitted_profiles": len(residuals),
        "uncovered_integer_residuals": residuals,
        "cell_audit": audit_cells,
        "assumptions": [
            "Only fixed affine witness parts are equated; common normalization cancels exactly.",
            "All witnesses were required to have identical subtract_normalization behavior.",
            "Every mesh edge uses the complete affine lower envelope, not only endpoint winners.",
            "Binary64 witness constants and charges make this a tuning diagnostic, not a certificate.",
            "A rounded point is retained only if its full-envelope best two are the crossing pair.",
            "Every emitted profile has exact support, fixed mass, and physical weight at least 21.",
        ],
    }


def run_self_test() -> None:
    witnesses = [
        Witness("left", 0.0, (0.0, 1.0, 0.0, 0.0, 0.0), frozenset(range(5)), True),
        Witness("right", -12.0, (0.0, -1.0, 0.0, 0.0, 0.0), frozenset(range(5)), True),
    ]
    left = tuple(map(Fraction, (1, 1, 1, 1, 96)))
    right = tuple(map(Fraction, (1, 11, 1, 1, 86)))
    transitions = lower_envelope_transitions(witnesses, left, right)
    if transitions != [(Fraction(1, 2), 1, 0)]:
        raise SystemExit(f"transition anchor self-test: affine envelope mismatch {transitions}")
    parameter = transitions[0][0]
    profile = round_exact_support(
        tuple((a + b) / 2 for a, b in zip(left, right)), tuple(range(5)), 100
    )
    if profile != (1, 6, 1, 1, 91):
        raise SystemExit(f"transition anchor self-test: rounding mismatch {profile}")
    print("self_test_crossing_parameter=1/2")
    print("status=PASS_G4_WITNESS_DOMINANCE_TRANSITION_ANCHOR_SELF_TEST")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path)
    parser.add_argument("--atlas", type=Path, action="append", default=[])
    parser.add_argument("--support-mask", type=lambda value: int(value, 0), default=31)
    parser.add_argument("--max-cells", type=int, default=128)
    parser.add_argument("--max-profiles", type=int, default=512)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if args.mesh is None or not args.atlas or args.output is None:
        raise SystemExit("transition anchors: --mesh, --atlas, and --output are required")
    if not 1 <= args.support_mask <= 31 or args.max_cells <= 0 or args.max_profiles <= 0:
        raise SystemExit("transition anchors: invalid search parameter")
    payload = json.loads(args.mesh.read_text(encoding="utf-8"))
    witnesses, sources, _anchors, duplicate_rows = load_witnesses(args.atlas)
    report = generate(
        payload, witnesses, args.support_mask, args.max_cells, args.max_profiles
    )
    report["mesh"] = {"path": str(args.mesh.resolve()), "sha256": sha256(args.mesh)}
    report["sources"] = sources
    report["duplicate_atlas_anchor_rows"] = duplicate_rows
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"selected_cells={report['selected_cells']}")
    print(f"transition_profiles={report['emitted_profiles']}")
    print("status=DIAGNOSTIC_G4_WITNESS_DOMINANCE_TRANSITION_ANCHORS")


if __name__ == "__main__":
    main()
