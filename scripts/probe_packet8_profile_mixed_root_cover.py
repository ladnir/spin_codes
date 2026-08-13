#!/usr/bin/env python3
"""Find one convex mixture of fixed witnesses covering the full profile polytope.

If ``L_i(a)`` is any logged fixed-witness upper bound, then the true logged
quantity is at most every ``L_i(a)`` and hence at most every convex average of
them.  All fixed packet-profile witnesses have the form

    constant_i - <a, charge_i> - log2 Q(a),

so a fixed convex average is convex in ``a``.  It therefore suffices to check
the vertices of the simplex clipped by the minimum outer-weight halfspace.

This script discovers the minimax mixture in binary64.  It is diagnostic, not
an outward-rounded certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from probe_packet8_profile_convex_cell_cover import Cell, cell_vertices
from probe_packet8_profile_ordered_chambers import add_full_bijection
from probe_packet8_profile_simplex_landscape import M, PROFILE_COUNT_LOG2, normalization_log2
from probe_packet8_shared_witness_io import load_shared_witness_arrays


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upgraded-cache", type=Path, required=True)
    parser.add_argument(
        "--extra-shared-report", type=Path, action="append", default=[]
    )
    parser.add_argument("--minimum-outer-weight", type=int, default=21)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    names, constants, charges = load_shared_witness_arrays(
        args.upgraded_cache, args.extra_shared_report
    )
    full = add_full_bijection([])[0]
    names.append(full.name)
    constants = np.append(constants, full.constant_log2)
    charges = np.vstack((charges, full.linear_charge))

    root = Cell((0,) * 9, (M,) * 9)
    vertices = cell_vertices(root, args.minimum_outer_weight)
    normalizations = normalization_log2(vertices)
    values = constants[None, :] - vertices @ charges.T - normalizations[:, None]

    witness_count = len(constants)
    objective = np.zeros(witness_count + 1)
    objective[-1] = 1.0
    inequalities = np.hstack((values, -np.ones((len(vertices), 1))))
    result = linprog(
        objective,
        A_ub=inequalities,
        b_ub=np.zeros(len(vertices)),
        A_eq=np.asarray([[1.0] * witness_count + [0.0]]),
        b_eq=np.asarray([1.0]),
        bounds=[(0.0, None)] * witness_count + [(None, None)],
        method="highs",
    )
    if not result.success:
        raise SystemExit(f"mixed root cover LP failed: {result.message}")

    weights = result.x[:-1]
    mixed_values = values @ weights
    target = -40.0 - PROFILE_COUNT_LOG2
    active = np.flatnonzero(weights > 1e-10)
    report = {
        "status": "DIAGNOSTIC_BINARY64_GLOBAL_CONVEX_MIXTURE",
        "target_log2": target,
        "minimum_outer_weight": args.minimum_outer_weight,
        "vertices": vertices.tolist(),
        "vertex_values_log2": mixed_values.tolist(),
        "maximum_vertex_value_log2": float(np.max(mixed_values)),
        "margin_bits": float(target - np.max(mixed_values)),
        "active_witnesses": [
            {
                "name": names[int(index)],
                "weight": float(weights[index]),
                "constant_log2": float(constants[index]),
                "charge": charges[index].tolist(),
            }
            for index in active
        ],
        "weight_sum": float(np.sum(weights)),
        "lp_objective_log2": float(result.fun),
        "global_profile_cover": bool(np.max(mixed_values) <= target),
    }
    print("packet-8 global mixed-witness root cover")
    print(f"witnesses={witness_count} root_vertices={len(vertices)}")
    print(f"active_witnesses={len(active)}")
    print(f"maximum_vertex_value_log2={np.max(mixed_values):.12f}")
    print(f"target_log2={target:.12f}")
    print(f"margin_bits={target - np.max(mixed_values):.12f}")
    for index in active:
        print(f"weight={weights[index]:.17g} witness={names[int(index)]}")
    print(
        "global_profile_cover="
        + ("PASS" if np.max(mixed_values) <= target else "NO")
    )
    print("status=DIAGNOSTIC_BINARY64_GLOBAL_CONVEX_MIXTURE")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
