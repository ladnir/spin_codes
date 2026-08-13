#!/usr/bin/env python3
"""Polynomial global profile cover from a mixture of full-support witnesses.

For a frozen inner/outer Cauchy witness the per-profile log bound has form

    F(a) = C - <a, q> - log2 Q(a).

Any convex combination of valid log bounds is again an upper bound because
``min_i F_i(a) <= sum_i lambda_i F_i(a)``.  The combined function is convex
in ``a``.  On the packet-profile simplex clipped by physical weight at least
21, it is therefore enough to check the nonzero pure vertices and the
zero-to-pure edge intersections with the weight-21 plane: only ``2g`` points.

Each component uses a class-mass floor ``epsilon / binom(g,j)`` so it has full
support, one favored pure class, the shared-drive inner bound, and a frozen
exact-spectrum total-weight pole.  The graph replacement interval is enlarged
to all offsets -128..128, making the frozen outer bound affine in physical
weight even at the clipped endpoints.

This is a binary64 discovery diagnostic.  Its finite LP and every selected
component must still be recomputed outward for a certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.special import gammaln, logsumexp

from packet_group_drive_stratified import (
    block_histograms,
    point_caps,
    profile_classes,
    profile_count,
)
from packet_group_outer_profile import (
    D,
    K,
    N,
    atom_count,
    full_spectrum,
    normalization_log2,
    total_weight_outer,
)
from packet_group_profile_bound import split_cap_table
from probe_packet8_drive_stratified_transfer import SharedDriveStratifiedKernel
from probe_packet8_repeated_value_state_transfer import witness
from probe_packet8_weight_profile_state_transfer import INNER_BLOCKS


MINIMUM_PHYSICAL_WEIGHT = 21
GRAPH_REPLACEMENTS = 128


def clipped_simplex_vertices(group_bits: int) -> tuple[np.ndarray, list[str]]:
    atoms = atom_count(group_bits)
    rows = []
    names = []
    for weight in range(1, group_bits + 1):
        pure = np.zeros(group_bits + 1, dtype=np.float64)
        pure[weight] = atoms
        rows.append(pure)
        names.append(f"pure_{weight}")

        boundary = np.zeros(group_bits + 1, dtype=np.float64)
        boundary[weight] = MINIMUM_PHYSICAL_WEIGHT / weight
        boundary[0] = atoms - boundary[weight]
        rows.append(boundary)
        names.append(f"weight21_edge_0_{weight}")
    return np.vstack(rows), names


def weight_slab_vertices(
    group_bits: int, low_average: float, high_average: float
) -> tuple[np.ndarray, list[str]]:
    """Vertices of a simplex intersected with one average-weight slab."""

    atoms = atom_count(group_bits)
    rows: dict[tuple[float, ...], tuple[np.ndarray, str]] = {}

    def add_boundary(average: float, label: str) -> None:
        rounded = round(average)
        if abs(average - rounded) < 1e-12 and 0 <= rounded <= group_bits:
            point = np.zeros(group_bits + 1, dtype=np.float64)
            point[int(rounded)] = atoms
            rows[tuple(point)] = (point, f"{label}_pure_{int(rounded)}")
        for left in range(group_bits + 1):
            if left >= average:
                break
            for right in range(left + 1, group_bits + 1):
                if right <= average:
                    continue
                point = np.zeros(group_bits + 1, dtype=np.float64)
                point[left] = atoms * (right - average) / (right - left)
                point[right] = atoms - point[left]
                key = tuple(round(float(value), 12) for value in point)
                rows[key] = (point, f"{label}_edge_{left}_{right}")

    add_boundary(low_average, "low")
    add_boundary(high_average, "high")
    values = list(rows.values())
    return np.vstack([row[0] for row in values]), [row[1] for row in values]


def solve_mixture(values: np.ndarray):
    component_count, vertex_count = values.shape
    objective = np.zeros(component_count + 1, dtype=np.float64)
    objective[-1] = 1.0
    inequalities = np.hstack(
        (values.T, -np.ones((vertex_count, 1), dtype=np.float64))
    )
    return linprog(
        objective,
        A_ub=inequalities,
        b_ub=np.zeros(vertex_count, dtype=np.float64),
        A_eq=np.asarray(
            [[1.0] * component_count + [0.0]], dtype=np.float64
        ),
        b_eq=np.asarray([1.0]),
        bounds=[(0.0, 1.0)] * component_count + [(None, None)],
        method="highs",
    )


def full_offset_outer_constant(log_pole: float) -> float:
    spectrum = full_spectrum()
    weights = np.asarray(
        [weight for weight, count in enumerate(spectrum) if count],
        dtype=np.float64,
    )
    logs = np.asarray(
        [math.log(spectrum[int(weight)]) for weight in weights],
        dtype=np.float64,
    )
    local = logsumexp(logs + weights * log_pole)
    offsets = np.arange(
        -GRAPH_REPLACEMENTS,
        GRAPH_REPLACEMENTS + 1,
        dtype=np.float64,
    )
    replacement = logsumexp(-offsets * log_pole)
    return ((K // 64) * local + replacement) / math.log(2.0)


def load_pure_rows(path: Path, group_bits: int) -> dict[int, dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    reports = data if isinstance(data, list) else [data]
    report = next(row for row in reports if int(row["group_bits"]) == group_bits)
    return {int(row["class_weight"]): row for row in report["rows"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--pure-curve", type=Path)
    parser.add_argument("--class-mass-floor", type=float, default=0.1)
    parser.add_argument("--poles", default="0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6")
    parser.add_argument("--screen-iterations", type=int, default=8)
    parser.add_argument(
        "--heuristic-poles",
        action="store_true",
        help="use the edge-distance pole schedule instead of a per-class grid screen",
    )
    parser.add_argument(
        "--retain-pole-grid",
        action="store_true",
        help="retain every verified inner pole as an LP component",
    )
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument(
        "--component-cache",
        type=Path,
        action="append",
        help="reuse the frozen component witnesses from a previous report",
    )
    parser.add_argument("--slabs-per-class", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.class_mass_floor <= 0.0:
        raise SystemExit("global mixture: class-mass floor must be positive")

    g = args.group
    atoms = atom_count(g)
    classes = np.asarray(profile_classes(g), dtype=np.float64)
    split_caps = split_cap_table()
    poles = [float(value) for value in args.poles.split(",")]
    if args.component_cache:
        components = []
        for cache_path in args.component_cache:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if int(cached["group_bits"]) != g:
                raise SystemExit("global mixture: component cache group mismatch")
            components.extend(cached["components"])
        favored_classes = ()
    else:
        components = []
        favored_classes = range(0, g + 1)
    for favored in favored_classes:
        profile = [0] * (g + 1)
        profile[favored] = atoms
        fugacities = args.class_mass_floor / classes
        # Normalize aggregate class mass, not one representative atom.  Thus
        # the favored class has total fugacity one and every other class has
        # total fugacity ``class_mass_floor``.
        fugacities[favored] = 1.0 / classes[favored]
        histograms = block_histograms(g, fugacities)
        caps = point_caps(g, fugacities)

        def inner_row(pole: float, iterations: int):
            eigenvalue, domination, values, worst_state = witness(
                SharedDriveStratifiedKernel(
                    histograms, caps, split_caps, pole
                ),
                iterations,
            )
            mgf = (
                math.log2(domination)
                + INNER_BLOCKS * math.log2(eigenvalue)
                + math.log2(float(values[0]))
            )
            return (
                mgf - D * math.log2(pole),
                pole,
                {
                    "inner_mgf_log2": mgf,
                    "charge_log2": 0.0,
                    "normalization_log2": atoms * math.log2(classes[favored]),
                    "lambda_log2": math.log2(eigenvalue),
                    "domination_log2": math.log2(domination),
                    "worst_state": int(worst_state),
                },
            )

        if args.retain_pole_grid:
            selected_inner_rows = [
                inner_row(pole, args.iterations) for pole in poles
            ]
        elif args.heuristic_poles:
            edge_distance = min(favored, g - favored)
            if edge_distance <= 1:
                pole = 0.6
            elif edge_distance <= 4:
                pole = 0.4
            elif edge_distance <= 10:
                pole = 0.2
            else:
                pole = 0.1
            selected_inner_rows = [inner_row(pole, args.iterations)]
        else:
            screened = min(
                (inner_row(pole, args.screen_iterations) for pole in poles),
                key=lambda row: row[0],
            )
            selected_inner_rows = [inner_row(screened[1], args.iterations)]
        outer_anchor_weight = (
            MINIMUM_PHYSICAL_WEIGHT if favored == 0 else atoms * favored
        )
        _, outer_result, _ = total_weight_outer(outer_anchor_weight)
        outer_log_pole = float(outer_result.x)
        for _, pole, inner_details in selected_inner_rows:
            constant = (
                full_offset_outer_constant(outer_log_pole)
                + float(inner_details["inner_mgf_log2"])
                - D * math.log2(pole)
            )
            charge = (
                np.arange(g + 1, dtype=np.float64)
                * outer_log_pole
                / math.log(2.0)
                + np.log2(fugacities)
            )
            components.append(
                {
                    "favored_class": favored,
                    "pole": pole,
                    "outer_log_pole": outer_log_pole,
                    "constant_log2": constant,
                    "charge": charge.tolist(),
                    "fugacities": fugacities.tolist(),
                    "inner_details": inner_details,
                }
            )
            print(
                f"component={favored}/{g} pole={pole:.6g} "
                f"outer_log_pole={outer_log_pole:.9g}",
                flush=True,
            )

    vertices, vertex_names = clipped_simplex_vertices(g)
    normalizations = normalization_log2(g, vertices)
    constants = np.asarray([row["constant_log2"] for row in components])
    charges = np.asarray([row["charge"] for row in components])
    values = constants[:, None] - charges @ vertices.T - normalizations[None, :]

    result = solve_mixture(values)
    if not result.success:
        raise SystemExit(f"global mixture LP failed: {result.message}")
    component_count = len(components)
    mixture = result.x[:component_count]
    mixed_values = mixture @ values
    worst = int(np.argmax(mixed_values))
    target = -40.0 - math.log2(profile_count(g, N))
    if args.slabs_per_class <= 0:
        raise SystemExit("global mixture: slabs per class must be positive")
    inner_constants = np.asarray(
        [
            float(row["inner_details"]["inner_mgf_log2"])
            - D * math.log2(float(row["pole"]))
            for row in components
        ]
    )
    inner_log_fugacities = np.asarray(
        [np.log2(np.asarray(row["fugacities"])) for row in components]
    )
    class_weights = np.arange(g + 1, dtype=np.float64)
    slab_reports = []
    slab_count = g * args.slabs_per_class
    for slab in range(slab_count):
        nominal_low = slab / args.slabs_per_class
        low = max(MINIMUM_PHYSICAL_WEIGHT / atoms, nominal_low)
        high = (slab + 1) / args.slabs_per_class
        if high <= low:
            continue
        outer_anchor_weight = int(round(atoms * (low + high) / 2.0))
        outer_anchor_weight = min(N, max(MINIMUM_PHYSICAL_WEIGHT, outer_anchor_weight))
        _, local_outer_result, _ = total_weight_outer(outer_anchor_weight)
        local_outer_log_pole = float(local_outer_result.x)
        local_constants = (
            inner_constants + full_offset_outer_constant(local_outer_log_pole)
        )
        local_charges = (
            inner_log_fugacities
            + class_weights[None, :]
            * local_outer_log_pole
            / math.log(2.0)
        )
        slab_vertices, slab_names = weight_slab_vertices(g, low, high)
        slab_normalizations = normalization_log2(g, slab_vertices)
        slab_values = (
            local_constants[:, None]
            - local_charges @ slab_vertices.T
            - slab_normalizations[None, :]
        )
        slab_result = solve_mixture(slab_values)
        if not slab_result.success:
            raise SystemExit(
                f"global mixture slab {slab} LP failed: {slab_result.message}"
            )
        slab_mixture = slab_result.x[:component_count]
        mixed_slab_values = slab_mixture @ slab_values
        slab_worst = int(np.argmax(mixed_slab_values))
        slab_reports.append(
            {
                "slab": slab,
                "outer_log_pole": local_outer_log_pole,
                "low_average_weight": low,
                "high_average_weight": high,
                "vertices": len(slab_vertices),
                "worst_vertex_log2": float(mixed_slab_values[slab_worst]),
                "margin_log2": float(target - mixed_slab_values[slab_worst]),
                "passed": bool(mixed_slab_values[slab_worst] <= target),
                "worst_vertex": slab_names[slab_worst],
                "mixture": [
                    {
                        "favored_class": int(
                            components[index]["favored_class"]
                        ),
                        "weight": float(weight),
                    }
                    for index, weight in enumerate(slab_mixture)
                    if weight > 1e-10
                ],
            }
        )
    slab_worst = max(slab_reports, key=lambda row: row["worst_vertex_log2"])

    report = {
        "status": "DIAGNOSTIC_BINARY64_GLOBAL_PROFILE_MIXTURE_COVER",
        "group_bits": g,
        "class_mass_floor": args.class_mass_floor,
        "iterations": args.iterations,
        "screen_iterations": args.screen_iterations,
        "pole_grid": poles,
        "heuristic_poles": args.heuristic_poles,
        "retain_pole_grid": args.retain_pole_grid,
        "slabs_per_class": args.slabs_per_class,
        "target_log2": target,
        "worst_vertex_log2": float(mixed_values[worst]),
        "margin_log2": float(target - mixed_values[worst]),
        "passed": bool(mixed_values[worst] <= target),
        "worst_vertex": vertex_names[worst],
        "mixture": [
            {
                "favored_class": int(components[index]["favored_class"]),
                "weight": float(weight),
            }
            for index, weight in enumerate(mixture)
            if weight > 1e-10
        ],
        "vertices": [
            {"name": name, "value_log2": float(value)}
            for name, value in zip(vertex_names, mixed_values)
        ],
        "components": components,
        "lp_message": result.message,
        "slab_cover": {
            "passed": all(row["passed"] for row in slab_reports),
            "worst_margin_log2": float(
                min(row["margin_log2"] for row in slab_reports)
            ),
            "worst_slab": int(slab_worst["slab"]),
            "rows": slab_reports,
        },
    }
    rendered = json.dumps(report, indent=2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    print(
        f"worst={report['worst_vertex_log2']:.9f} target={target:.9f} "
        f"margin={report['margin_log2']:.9f} passed={report['passed']} "
        f"vertex={report['worst_vertex']} active={len(report['mixture'])}",
        flush=True,
    )
    print(
        f"slab_cover_passed={report['slab_cover']['passed']} "
        f"worst_slab={report['slab_cover']['worst_slab']} "
        f"worst_margin={report['slab_cover']['worst_margin_log2']:.9f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
