#!/usr/bin/env python3
"""Certify the minimum distance of a fixed [256,64] binary code by MILP."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


DEFAULT_INPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/fixed_nested_pair_depth2_search.json"
)
DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/fixed_nested_pair_depth2_min_distance.json"
)


def solve(source: dict[str, object], time_limit: float) -> dict[str, object]:
    columns = [
        int(value, 16) for value in source["matrices"]["A_columns_hex_256"]
    ]
    dimension = len(columns)
    length = 256
    q_start = 0
    y_start = dimension
    k_start = dimension + length
    variables = dimension + 2 * length

    row_indices: list[int] = []
    column_indices: list[int] = []
    coefficients: list[float] = []
    row_degrees = [0] * length
    for output in range(length):
        for input_index, column in enumerate(columns):
            if (column >> output) & 1:
                row_indices.append(output)
                column_indices.append(q_start + input_index)
                coefficients.append(-1.0)
                row_degrees[output] += 1
        row_indices.extend((output, output))
        column_indices.extend((y_start + output, k_start + output))
        coefficients.extend((1.0, 2.0))

    # Exclude the all-zero information word.
    nonzero_row = length
    for input_index in range(dimension):
        row_indices.append(nonzero_row)
        column_indices.append(q_start + input_index)
        coefficients.append(1.0)

    matrix = coo_matrix(
        (coefficients, (row_indices, column_indices)),
        shape=(length + 1, variables),
    ).tocsr()
    lower_constraints = np.zeros(length + 1)
    upper_constraints = np.zeros(length + 1)
    lower_constraints[-1] = 1.0
    upper_constraints[-1] = float(dimension)

    lower_bounds = np.zeros(variables)
    upper_bounds = np.ones(variables)
    for output, degree in enumerate(row_degrees):
        upper_bounds[k_start + output] = degree // 2

    objective = np.zeros(variables)
    objective[y_start : y_start + length] = 1.0
    integrality = np.ones(variables)
    result = milp(
        c=objective,
        integrality=integrality,
        bounds=Bounds(lower_bounds, upper_bounds),
        constraints=LinearConstraint(
            matrix, lower_constraints, upper_constraints
        ),
        options={
            "time_limit": time_limit,
            "mip_rel_gap": 0.0,
            "presolve": True,
        },
    )

    payload: dict[str, object] = {
        "schema": "binary-linear-code-min-distance-milp-v1",
        "source_schema": source["schema"],
        "source_chosen_seed": source["parameters"]["chosen_seed"],
        "length": length,
        "dimension": dimension,
        "solver_status": int(result.status),
        "solver_message": str(result.message),
        "success": bool(result.success),
        "mip_node_count": getattr(result, "mip_node_count", None),
        "mip_dual_bound": getattr(result, "mip_dual_bound", None),
        "mip_gap": getattr(result, "mip_gap", None),
        "time_limit_seconds": time_limit,
    }
    if result.x is not None:
        information = 0
        for index, value in enumerate(result.x[:dimension]):
            if value > 0.5:
                information |= 1 << index
        codeword = 0
        for index in range(dimension):
            if (information >> index) & 1:
                codeword ^= columns[index]
        payload.update(
            {
                "objective": int(round(float(result.fun))),
                "information_hex": f"{information:016x}",
                "information_weight": information.bit_count(),
                "codeword_hex": f"{codeword:064x}",
                "verified_codeword_weight": codeword.bit_count(),
                "optimality_certificate": (
                    "Exact minimum distance" if result.status == 0 else
                    "Feasible upper bound only; solver did not prove optimality"
                ),
            }
        )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--time-limit", type=float, default=120.0)
    args = parser.parse_args()
    source = json.loads(args.input.read_text(encoding="utf-8"))
    payload = solve(source, args.time_limit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"status,{payload['solver_status']},{payload['solver_message']}")
    if "objective" in payload:
        print(
            f"distance,{payload['objective']},"
            f"input_weight,{payload['information_weight']},"
            f"verified,{payload['verified_codeword_weight']}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
