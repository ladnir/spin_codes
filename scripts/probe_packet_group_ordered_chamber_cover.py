#!/usr/bin/env python3
"""Diagnostic ordered-chamber cover for generalized packet profiles.

For every permutation of the ``g+1`` class counts, the corresponding ordered
region of the profile simplex is itself a simplex.  Its vertices are uniform
profiles on the nested prefixes of that permutation.  Intersecting with the
physical-weight-21 halfspace retains the valid vertices and adds the cut-edge
intersections.  One fixed affine witness or one fixed convex mixture checked
at those vertices therefore covers the whole real chamber.

Closed chambers overlap on ordering boundaries.  The diagnostic global ledger
retains that harmless duplication and charges every chamber by the full
profile count, a deliberately coarse but simple upper bound.  Binary64 LP and
witness values are discovery evidence only until outward hardened.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import N, atom_count, normalization_log2
from probe_packet_group_support_face_cover import (
    MINIMUM_PHYSICAL_WEIGHT,
    parse_atlas,
    rounded_integer_profile,
    solve_fixed_mixture,
)


def chamber_vertices(
    group_bits: int, order: tuple[int, ...]
) -> list[tuple[str, tuple[Fraction, ...]]]:
    atoms = atom_count(group_bits)
    simplex = []
    for prefix_end in range(group_bits + 1):
        profile = [Fraction(0)] * (group_bits + 1)
        value = Fraction(atoms, prefix_end + 1)
        for index in order[: prefix_end + 1]:
            profile[index] = value
        simplex.append((f"prefix_{prefix_end + 1}", tuple(profile)))

    valid = [
        sum(index * value for index, value in enumerate(profile))
        >= MINIMUM_PHYSICAL_WEIGHT
        for _name, profile in simplex
    ]
    vertices: dict[tuple[Fraction, ...], str] = {}
    for (name, profile), keep in zip(simplex, valid):
        if keep:
            vertices[profile] = name

    for left, ((left_name, left_profile), left_valid) in enumerate(zip(simplex, valid)):
        if left_valid:
            continue
        left_weight = sum(index * value for index, value in enumerate(left_profile))
        for right, ((right_name, right_profile), right_valid) in enumerate(zip(simplex, valid)):
            if not right_valid or left == right:
                continue
            right_weight = sum(index * value for index, value in enumerate(right_profile))
            if right_weight <= left_weight:
                continue
            toward_right = Fraction(
                MINIMUM_PHYSICAL_WEIGHT - left_weight,
                right_weight - left_weight,
            )
            if not Fraction(0) <= toward_right <= Fraction(1):
                continue
            profile = tuple(
                left_value * (Fraction(1) - toward_right) + right_value * toward_right
                for left_value, right_value in zip(left_profile, right_profile)
            )
            vertices[profile] = f"weight21_{left_name}_{right_name}"
    return sorted(((name, row) for row, name in vertices.items()), key=lambda item: item[0])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--atlas", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    g = args.group
    witnesses = parse_atlas(args.atlas, g)
    constants = np.asarray([row["constant_log2"] for row in witnesses])
    charges = np.asarray([row["charge"] for row in witnesses])
    total_profiles = profile_count(g, N)
    chambers = math.factorial(g + 1)
    per_chamber_target = -40.0 - math.log2(total_profiles) - math.log2(chambers)
    rows = []
    residuals = []
    terms = []

    for order in itertools.permutations(range(g + 1)):
        vertices = chamber_vertices(g, order)
        matrix = np.asarray(
            [[float(value) for value in profile] for _name, profile in vertices]
        )
        normalizations = normalization_log2(g, matrix)
        values = constants[:, None] - charges @ matrix.T - normalizations[None, :]
        finite = np.flatnonzero(np.all(np.isfinite(values), axis=1))
        result = solve_fixed_mixture(values[finite])
        if not result.success:
            raise RuntimeError(f"ordered chamber {order} LP failed: {result.message}")
        weights = result.x[: len(finite)]
        mixed = weights @ values[finite]
        worst_index = int(np.argmax(mixed))
        maximum = float(mixed[worst_index])
        contribution = maximum + math.log2(total_profiles)
        terms.append(contribution)
        worst_name, worst_fractional = vertices[worst_index]
        worst_profile = rounded_integer_profile(worst_fractional)
        row = {
            "order": list(order),
            "vertices": len(vertices),
            "maximum_vertex_log2": maximum,
            "per_chamber_target_log2": per_chamber_target,
            "target_margin_bits": per_chamber_target - maximum,
            "chamber_profile_count_upper": total_profiles,
            "chamber_contribution_log2": contribution,
            "worst_vertex": worst_name,
            "worst_integer_profile": worst_profile,
            "mixture": [
                {
                    "witness": witnesses[int(finite[index])]["name"],
                    "weight": float(weight),
                }
                for index, weight in enumerate(weights)
                if weight > 1e-10
            ],
            "passed_allocated_target": bool(maximum <= per_chamber_target),
        }
        rows.append(row)
        if maximum > per_chamber_target:
            residuals.append(
                {
                    "name": "chamber_" + "_".join(map(str, order)) + "_worst",
                    "profile": worst_profile,
                    "value_log2": maximum,
                    "target_log2": per_chamber_target,
                }
            )

    global_union = float(
        logsumexp(np.asarray(terms) * math.log(2.0)) / math.log(2.0)
    )
    worst = max(rows, key=lambda row: row["chamber_contribution_log2"])
    report = {
        "status": "DIAGNOSTIC_BINARY64_ORDERED_CHAMBER_CONVEX_COVER",
        "group_bits": g,
        "chambers": chambers,
        "witnesses": len(witnesses),
        "profile_count": total_profiles,
        "per_chamber_target_log2": per_chamber_target,
        "complete_chamber_cover": not residuals,
        "global_union_log2": global_union,
        "global_security_margin_bits": -global_union,
        "global_slack_beyond_40_bits": -40.0 - global_union,
        "passed_global_40_bits": bool(global_union <= -40.0),
        "worst_chamber_contribution": worst,
        "chamber_rows": rows,
        "uncovered_integer_residuals": residuals,
        "sources": [str(path) for path in args.atlas],
        "assumptions": [
            "binary64 witness evaluation and LP selection are diagnostic",
            "each closed chamber is charged the full global profile count",
            "the clipped ordered chamber contains every feasible profile in that ordering",
        ],
    }
    rendered = json.dumps(report, indent=2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    print(f"chambers={chambers} witnesses={len(witnesses)}")
    print(f"residual_chambers={len(residuals)}")
    print(f"global_union_log2={global_union:.12f}")
    print(f"global_security_margin_bits={-global_union:.12f}")
    print(f"passed_global_40_bits={global_union <= -40.0}")


if __name__ == "__main__":
    main()
