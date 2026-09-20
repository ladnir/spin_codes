#!/usr/bin/env python3
"""Export an exactly equivalent, binomial-scaled BCH sandwich LP.

Each half-spectrum count x_w is replaced by S_w y_w, where S_w is an exact
integer near its random-code scale.  Each constraint is then divided by its
largest exact coefficient or right-hand side.  Both transformations preserve
the feasible set over the rationals and improve the numerical bootstrap used
by exact simplex solvers.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from export_lp import independent_oa_constraints


ROOT = Path(__file__).resolve().parents[1]


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def expression(terms: list[tuple[Fraction, str]]) -> str:
    pieces: list[str] = []
    for coefficient, variable in terms:
        if not coefficient:
            continue
        sign = "+" if coefficient > 0 else "-"
        magnitude = abs(coefficient)
        scalar = "" if magnitude == 1 else fraction_text(magnitude) + " "
        pieces.append(f"{sign} {scalar}{variable}")
    return " ".join(pieces) if pieces else "0"


def variable_scales(model: dict) -> dict[str, int]:
    metadata = model["metadata"]
    n = int(metadata["n"])
    q_dimension = int(metadata["Q_dimension"])
    denominator = 1 << (n - 1 - q_dimension)
    scales: dict[str, int] = {}
    for name in metadata["variables"]:
        weight = int(name.split("_", 1)[1])
        scales[name] = max(1, math.comb(n, weight) // denominator)
    return scales


def objective_spec(objective: str, variables: list[str], scales: dict[str, int]):
    if objective in variables:
        return {objective: 1}, scales[objective]
    if objective.startswith("c_"):
        weight = objective.split("_", 1)[1]
        q_name, h_name = f"q_{weight}", f"h_{weight}"
        if q_name in variables and h_name in variables:
            if scales[q_name] != scales[h_name]:
                raise AssertionError("q and h shell scales differ")
            return {q_name: 1, h_name: 31}, scales[q_name]
    raise ValueError(f"unknown objective {objective}")


def export_scaled(model_path: Path, objective: str, output: Path) -> dict:
    model = json.loads(model_path.read_text(encoding="utf-8"))
    variables = model["metadata"]["variables"]
    scales = variable_scales(model)
    objective_coefficients, objective_multiplier = objective_spec(
        objective, variables, scales
    )
    constraints = independent_oa_constraints(model)

    objective_terms = [
        (Fraction(coefficient), name)
        for name, coefficient in objective_coefficients.items()
    ]
    lines = ["Maximize", f" obj: {expression(objective_terms)}", "Subject To"]
    emitted = 0
    for row in constraints:
        if row["name"].startswith(("q_nonneg_", "h_nonneg_")):
            continue
        scaled_terms = [
            (int(coefficient) * scales[name], name)
            for name, coefficient in row["coeffs"].items()
        ]
        rhs = int(row["rhs"])
        row_scale = max(
            1,
            abs(rhs),
            max((abs(coefficient) for coefficient, _ in scaled_terms), default=0),
        )
        normalized = [
            (Fraction(coefficient, row_scale), name)
            for coefficient, name in scaled_terms
        ]
        operator = {"eq": "=", "le": "<=", "ge": ">="}[row["sense"]]
        lines.append(
            f" {expression(normalized)} {operator} {fraction_text(Fraction(rhs, row_scale))}"
        )
        emitted += 1

    lines.append("Bounds")
    lines.extend(f" 0 <= {name}" for name in variables)
    lines.append("End")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="ascii")

    receipt = {
        "classification": "exact rational column- and row-scaled LP",
        "source_model": str(model_path.resolve()),
        "objective_label": objective,
        "objective_scaled_coefficients": objective_coefficients,
        "objective_physical_multiplier": objective_multiplier,
        "variables": len(variables),
        "linear_rows": emitted,
        "variable_scales": scales,
        "lp_path": str(output.resolve()),
        "equivalence": (
            "physical variable x_name = variable_scales[name] * LP variable name; "
            "every row was divided by a positive exact integer"
        ),
    }
    output.with_suffix(".json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model",
        type=Path,
        default=ROOT / "generated" / "coupled_lp_exact.json",
    )
    parser.add_argument("--objective", default="h_38")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "generated" / "max_h_38_scaled_rational.lp",
    )
    args = parser.parse_args()
    print(json.dumps(export_scaled(args.model, args.objective, args.output), indent=2))


if __name__ == "__main__":
    main()
