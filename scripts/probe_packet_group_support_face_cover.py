#!/usr/bin/env python3
"""Diagnostic support-stratified convex cover for generalized packet profiles.

Integer profiles are partitioned by their exact nonzero class support.  After
subtracting one from every active class, one support stratum is a simplex
clipped by the physical-weight-21 halfspace.  A fixed affine Cauchy witness,
or one fixed convex mixture of such witnesses, is convex after the Gamma
normalization extension, so checking the explicit continuous-hull vertices
covers every integer profile in that support.

The final diagnostic sums each support maximum against the safe positive-
composition count ``C(M-1,s-1)``.  This count deliberately includes any
profiles removed by the weight-21 clip.  Arithmetic and LP selection are
binary64 discovery tools; selected witnesses and mixture weights still need
outward recomputation before theorem use.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.special import logsumexp

from packet_group_drive_stratified import profile_count
from packet_group_outer_profile import D, K, N, atom_count, normalization_log2


MINIMUM_PHYSICAL_WEIGHT = 21


def support_hull_vertices(
    group_bits: int, support: tuple[int, ...]
) -> list[tuple[str, tuple[Fraction, ...]]]:
    """Vertices of the positive-support simplex after the weight-21 clip."""

    atoms = atom_count(group_bits)
    if not support or support == (0,):
        return []
    base_weight = sum(support)
    residual = atoms - len(support)
    if residual < 0:
        return []

    pure = {}
    valid = {}
    for favored in support:
        profile = [Fraction(0)] * (group_bits + 1)
        for index in support:
            profile[index] = Fraction(1)
        profile[favored] += residual
        row = tuple(profile)
        pure[favored] = row
        valid[favored] = sum(index * value for index, value in enumerate(row)) >= 21

    vertices: dict[tuple[Fraction, ...], str] = {}
    for favored in support:
        if valid[favored]:
            vertices[pure[favored]] = f"pure_residual_{favored}"

    # A simplex clipped by one halfspace has, in addition to retained
    # vertices, exactly the intersections of cut simplex edges.
    for invalid_class in support:
        if valid[invalid_class]:
            continue
        invalid_weight = base_weight + residual * invalid_class
        for valid_class in support:
            if not valid[valid_class]:
                continue
            denominator = residual * (valid_class - invalid_class)
            if denominator <= 0:
                continue
            toward_valid = Fraction(MINIMUM_PHYSICAL_WEIGHT - invalid_weight, denominator)
            if not Fraction(0) <= toward_valid <= Fraction(1):
                continue
            profile = [Fraction(0)] * (group_bits + 1)
            for index in support:
                profile[index] = Fraction(1)
            profile[invalid_class] += residual * (Fraction(1) - toward_valid)
            profile[valid_class] += residual * toward_valid
            row = tuple(profile)
            if sum(index * value for index, value in enumerate(row)) != 21:
                raise RuntimeError("support hull edge intersection missed weight 21")
            vertices[row] = f"weight21_edge_{invalid_class}_{valid_class}"

    return sorted(((name, row) for row, name in vertices.items()), key=lambda item: item[0])


def parse_atlas(paths: list[Path], group_bits: int) -> list[dict]:
    rows = []
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        candidates = payload.get("rows", payload) if isinstance(payload, dict) else payload
        if not isinstance(candidates, list):
            raise ValueError(f"support cover: atlas {path} has no witness rows")
        for index, row in enumerate(candidates):
            if int(row.get("group_bits", group_bits)) != group_bits:
                raise ValueError(f"support cover: group mismatch in {path}:{index}")
            charge = row.get("charge")
            if charge is None or len(charge) != group_bits + 1:
                continue
            if not math.isfinite(float(row["constant_log2"])):
                continue
            rows.append(
                {
                    "name": f"{path.name}:{index}",
                    "display_name": str(row.get("name", f"witness_{index}")),
                    "constant_log2": float(row["constant_log2"]),
                    "charge": [float(value) for value in charge],
                    "source": str(path),
                }
            )
    if not rows:
        raise ValueError("support cover: no affine witnesses loaded")

    # The exact inner bijection/Hamming-ball moment is a valid full-support
    # affine branch with zero charge and the same profile normalization.
    rows.append(
        {
            "name": "full_bijection",
            "display_name": "full_bijection",
            "constant_log2": K + N * math.log2(11.0) - (N - D) * math.log2(10.0),
            "charge": [0.0] * (group_bits + 1),
            "source": "built-in bijection moment",
        }
    )
    return rows


def solve_fixed_mixture(values: np.ndarray):
    """Minimize the maximum vertex value using one fixed convex mixture."""

    components, vertices = values.shape
    objective = np.zeros(components + 1)
    objective[-1] = 1.0
    inequalities = np.hstack((values.T, -np.ones((vertices, 1))))
    return linprog(
        objective,
        A_ub=inequalities,
        b_ub=np.zeros(vertices),
        A_eq=np.asarray([[1.0] * components + [0.0]]),
        b_eq=np.asarray([1.0]),
        bounds=[(0.0, None)] * components + [(None, None)],
        method="highs-ds",
    )


def rounded_integer_profile(profile: tuple[Fraction, ...]) -> list[int]:
    floors = [value.numerator // value.denominator for value in profile]
    remainder = sum(profile) - sum(floors)
    if remainder.denominator != 1:
        raise RuntimeError("profile rounding remainder is not integral")
    order = sorted(
        range(len(profile)),
        key=lambda index: profile[index] - floors[index],
        reverse=True,
    )
    for index in order[: int(remainder)]:
        floors[index] += 1
    if sum(floors) != int(sum(profile)):
        raise RuntimeError("rounded profile has the wrong mass")
    if sum(index * count for index, count in enumerate(floors)) < 21:
        low = next(index for index, count in enumerate(floors) if count > 1)
        high = max(index for index, value in enumerate(profile) if value > 0)
        floors[low] -= 1
        floors[high] += 1
    return floors


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--atlas", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    g = args.group
    atoms = atom_count(g)
    witnesses = parse_atlas(args.atlas, g)
    constants = np.asarray([row["constant_log2"] for row in witnesses])
    charges = np.asarray([row["charge"] for row in witnesses])
    uniform_target = -40.0 - math.log2(profile_count(g, N))
    reports = []
    residuals = []
    aggregate_terms = []

    for mask in range(1, 1 << (g + 1)):
        support = tuple(index for index in range(g + 1) if mask & (1 << index))
        vertices = support_hull_vertices(g, support)
        if not vertices:
            continue
        matrix = np.asarray(
            [[float(value) for value in profile] for _name, profile in vertices]
        )
        normalizations = normalization_log2(g, matrix)
        values = constants[:, None] - charges @ matrix.T - normalizations[None, :]
        finite = np.flatnonzero(np.all(np.isfinite(values), axis=1))
        if not len(finite):
            raise RuntimeError(f"support {support} has no finite witness")
        result = solve_fixed_mixture(values[finite])
        if not result.success:
            raise RuntimeError(f"support {support} mixture LP failed: {result.message}")
        weights = result.x[: len(finite)]
        mixed = weights @ values[finite]
        worst_index = int(np.argmax(mixed))
        maximum = float(mixed[worst_index])
        count_upper = math.comb(atoms - 1, len(support) - 1)
        contribution = maximum + math.log2(count_upper)
        aggregate_terms.append(contribution)
        active = [
            {
                "witness": witnesses[int(finite[index])]["name"],
                "weight": float(weight),
            }
            for index, weight in enumerate(weights)
            if weight > 1e-10
        ]
        worst_name, worst_fractional = vertices[worst_index]
        worst_profile = rounded_integer_profile(worst_fractional)
        row = {
            "support": list(support),
            "vertices": len(vertices),
            "maximum_vertex_log2": maximum,
            "uniform_per_profile_target_log2": uniform_target,
            "uniform_target_margin_bits": uniform_target - maximum,
            "positive_composition_count_upper": count_upper,
            "support_contribution_log2": contribution,
            "worst_vertex": worst_name,
            "worst_integer_profile": worst_profile,
            "mixture": active,
            "passed_uniform_target": bool(maximum <= uniform_target),
        }
        reports.append(row)
        if maximum > uniform_target:
            residuals.append(
                {
                    "name": "support_" + "_".join(map(str, support)) + "_worst",
                    "profile": worst_profile,
                    "value_log2": maximum,
                    "target_log2": uniform_target,
                }
            )

    # scipy.logsumexp uses natural exponentials, while aggregate_terms are
    # base-two logs, so scale by ln(2) in both directions.
    global_union = float(
        logsumexp(np.asarray(aggregate_terms) * math.log(2.0)) / math.log(2.0)
    )
    worst = max(reports, key=lambda row: row["support_contribution_log2"])
    report = {
        "status": "DIAGNOSTIC_BINARY64_SUPPORT_STRATIFIED_CONVEX_COVER",
        "group_bits": g,
        "atom_count": atoms,
        "supports": len(reports),
        "witnesses": len(witnesses),
        "uniform_per_profile_target_log2": uniform_target,
        "all_supports_pass_uniform_target": not residuals,
        "global_union_log2": global_union,
        "global_security_margin_bits": -global_union,
        "global_slack_beyond_40_bits": -40.0 - global_union,
        "passed_global_40_bits": bool(global_union <= -40.0),
        "worst_support_contribution": worst,
        "support_rows": reports,
        "uncovered_integer_residuals": residuals,
        "sources": [str(path) for path in args.atlas],
        "assumptions": [
            "binary64 witness evaluation and LP selection are diagnostic",
            "each support count uses the safe positive-composition upper bound",
            "continuous support hull includes every exact-support integer profile",
        ],
    }
    rendered = json.dumps(report, indent=2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    print(f"supports={len(reports)} witnesses={len(witnesses)}")
    print(f"uniform_residuals={len(residuals)}")
    print(f"global_union_log2={global_union:.12f}")
    print(f"global_security_margin_bits={-global_union:.12f}")
    print(f"passed_global_40_bits={global_union <= -40.0}")


if __name__ == "__main__":
    main()
