#!/usr/bin/env python3
"""Verify a QSopt_ex certificate for an exactly scaled BCH sandwich LP.

The checker reconstructs the rational LP from ``coupled_lp_exact.json`` rather
than trusting decimal solver output.  It verifies primal feasibility, dual
feasibility, dual signs, and exact equality of the primal and dual objectives.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from export_lp import independent_oa_constraints
from export_scaled_rational_lp import variable_scales


ROOT = Path(__file__).resolve().parents[1]


def parse_assignments(text: str, start: str, end: str) -> dict[str, Fraction]:
    section = text.split(start, 1)[1].split(end, 1)[0]
    values: dict[str, Fraction] = {}
    for line in section.splitlines():
        if " = " not in line:
            continue
        name, value = line.split(" = ", 1)
        values[name.strip()] = Fraction(value.strip())
    return values


def normalized_rows(model: dict, scales: dict[str, int]) -> list[dict]:
    rows = []
    for row in independent_oa_constraints(model):
        if row["name"].startswith(("q_nonneg_", "h_nonneg_")):
            continue
        coefficients = {
            name: Fraction(int(value) * scales[name])
            for name, value in row["coeffs"].items()
        }
        rhs = Fraction(int(row["rhs"]))
        row_scale = max(
            1,
            abs(rhs),
            *(abs(value) for value in coefficients.values()),
        )
        rows.append(
            {
                "name": row["name"],
                "coeffs": {
                    name: value / row_scale
                    for name, value in coefficients.items()
                },
                "sense": row["sense"],
                "rhs": rhs / row_scale,
            }
        )
    return rows


def verify(model_path: Path, lp_path: Path, solution_path: Path) -> dict:
    model = json.loads(model_path.read_text(encoding="utf-8"))
    receipt = json.loads(lp_path.with_suffix(".json").read_text(encoding="utf-8"))
    objective = receipt.get(
        "objective_label", receipt.get("objective_scaled_variable")
    )
    objective_coefficients = receipt.get(
        "objective_scaled_coefficients", {objective: 1}
    )
    scales = variable_scales(model)
    rows = normalized_rows(model, scales)

    solution_text = solution_path.read_text(encoding="ascii")
    if "status = OPTIMAL" not in solution_text:
        raise AssertionError("QSopt_ex solution is not marked OPTIMAL")
    primal = parse_assignments(solution_text, "VARS:", "REDUCED COST:")
    dual_sparse = parse_assignments(solution_text, "PI:", "SLACK:")
    variables = model["metadata"]["variables"]
    x = {name: primal.get(name, Fraction()) for name in variables}
    dual = [
        dual_sparse.get(f"c{index}", Fraction())
        for index in range(1, len(rows) + 1)
    ]

    if any(value < 0 for value in x.values()):
        raise AssertionError("primal nonnegativity failed")
    for row in rows:
        lhs = sum(
            (coefficient * x[name] for name, coefficient in row["coeffs"].items()),
            Fraction(),
        )
        relation = {
            "eq": lhs == row["rhs"],
            "le": lhs <= row["rhs"],
            "ge": lhs >= row["rhs"],
        }[row["sense"]]
        if not relation:
            raise AssertionError(f"primal row failed: {row['name']}")

    for multiplier, row in zip(dual, rows):
        if row["sense"] == "le" and multiplier < 0:
            raise AssertionError(f"negative <= multiplier: {row['name']}")
        if row["sense"] == "ge" and multiplier > 0:
            raise AssertionError(f"positive >= multiplier: {row['name']}")

    for name in variables:
        lhs = sum(
            (
                multiplier * row["coeffs"].get(name, Fraction())
                for multiplier, row in zip(dual, rows)
            ),
            Fraction(),
        )
        required = Fraction(objective_coefficients.get(name, 0))
        if lhs < required:
            raise AssertionError(f"dual column failed: {name}")

    primal_objective = sum(
        (
            Fraction(coefficient) * x[name]
            for name, coefficient in objective_coefficients.items()
        ),
        Fraction(),
    )
    dual_objective = sum(
        (multiplier * row["rhs"] for multiplier, row in zip(dual, rows)),
        Fraction(),
    )
    if primal_objective != dual_objective:
        raise AssertionError("primal and dual objectives differ")

    objective_multiplier = int(receipt["objective_physical_multiplier"])
    physical_optimum = primal_objective * objective_multiplier
    integer_cap = physical_optimum.numerator // physical_optimum.denominator
    lattice_modulus = 128 if objective == "h_38" else 1
    lattice_cap = integer_cap - integer_cap % lattice_modulus
    derived_c38_cap = 31 * lattice_cap if objective == "h_38" else None

    def digest(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    result = {
        "classification": "independently verified exact rational LP certificate",
        "model": str(model_path.resolve()),
        "lp": str(lp_path.resolve()),
        "solution": str(solution_path.resolve()),
        "objective_label": objective,
        "objective_scaled_coefficients": objective_coefficients,
        "objective_physical_multiplier": objective_multiplier,
        "scaled_optimum": {
            "numerator": str(primal_objective.numerator),
            "denominator": str(primal_objective.denominator),
        },
        "physical_optimum": {
            "numerator": str(physical_optimum.numerator),
            "denominator": str(physical_optimum.denominator),
        },
        "physical_integer_cap": integer_cap,
        "lattice_modulus": lattice_modulus,
        "physical_lattice_cap": lattice_cap,
        "derived_A38_C_cap": derived_c38_cap,
        "derived_A38_C_log2": (
            math.log2(derived_c38_cap) if derived_c38_cap is not None else None
        ),
        "rows_verified": len(rows),
        "variables_verified": len(variables),
        "nonzero_dual_multipliers": sum(value != 0 for value in dual),
        "checks": [
            "primal nonnegativity",
            "all primal rows",
            "dual row-sign restrictions",
            "all dual column inequalities",
            "exact primal-dual objective equality",
        ],
        "sha256": {"lp": digest(lp_path), "solution": digest(solution_path)},
    }
    certificate_path = solution_path.with_suffix(".certificate.json")
    certificate_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model",
        type=Path,
        default=ROOT / "generated" / "coupled_lp_exact.json",
    )
    parser.add_argument(
        "--lp",
        type=Path,
        default=ROOT / "generated" / "max_h_38_scaled_rational.lp",
    )
    parser.add_argument(
        "--solution",
        type=Path,
        default=ROOT / "generated" / "max_h_38_scaled_rational.sol",
    )
    args = parser.parse_args()
    print(json.dumps(verify(args.model, args.lp, args.solution), indent=2))


if __name__ == "__main__":
    main()
