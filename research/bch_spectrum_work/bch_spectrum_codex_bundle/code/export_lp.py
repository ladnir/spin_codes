"""Export the exact continuous BCH-sandwich relaxation in CPLEX LP syntax.

The generated file contains integer coefficients only and is suitable for an exact
rational LP solver such as QSopt_ex.  Congruence and integer-variable metadata from
the JSON interchange model are not part of the continuous LP file.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from spectrum_exact import symmetric_kraw_coeff


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = ROOT / "generated" / "coupled_lp_exact.json"


def linear_expression(coefficients: dict[str, str], order: dict[str, int]) -> str:
    terms = []
    for name, coefficient_text in sorted(
        coefficients.items(), key=lambda item: order[item[0]]
    ):
        coefficient = int(coefficient_text)
        if coefficient > 0:
            terms.append(f"+ {coefficient} {name}")
        elif coefficient < 0:
            terms.append(f"- {-coefficient} {name}")
    return " ".join(terms) if terms else "0"


def independent_oa_constraints(model):
    """Replace moment OA rows by their independent even Krawtchouk basis."""
    metadata = model["metadata"]
    n = metadata["n"]
    q_dimension = metadata["Q_dimension"]
    weights = [
        int(name.split("_", 1)[1])
        for name in metadata["variables"]
        if name.startswith("q_")
    ]
    constraints = [
        constraint
        for constraint in model["constraints"]
        if not constraint["name"].startswith(("Q_OA_t", "H_OA_t"))
    ]
    for label, prefix in (("Q", "q_"), ("H", "h_")):
        for degree in range(0, 16, 2):
            coefficients = {}
            for weight in weights:
                coefficient = symmetric_kraw_coeff(n, degree, weight)
                if coefficient:
                    coefficients[f"{prefix}{weight}"] = str(coefficient)
            constraints.append({
                "name": f"{label}_OA_Kraw_{degree}",
                "coeffs": coefficients,
                "sense": "eq",
                "rhs": str(1 << q_dimension if degree == 0 else 0),
            })
    return constraints


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--objective", default="h_38")
    parser.add_argument("--sense", choices=("min", "max"), default="max")
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--oa-basis",
        choices=("independent-krawtchouk", "moments"),
        default="independent-krawtchouk",
    )
    parser.add_argument(
        "--omit-prefix",
        action="append",
        default=[],
        help="omit constraints whose names start with this prefix (diagnostics)",
    )
    args = parser.parse_args()

    model = json.loads(args.model.read_text())
    variables = model["metadata"]["variables"]
    if args.objective not in variables:
        raise SystemExit(f"unknown objective variable: {args.objective}")
    order = {name: position for position, name in enumerate(variables)}

    output = args.output
    if output is None:
        output = ROOT / "generated" / f"{args.sense}_{args.objective}.lp"

    lines = [
        "Maximize" if args.sense == "max" else "Minimize",
        f" obj: {args.objective}",
        "Subject To",
    ]
    constraints = (
        independent_oa_constraints(model)
        if args.oa_basis == "independent-krawtchouk"
        else model["constraints"]
    )
    if args.omit_prefix:
        constraints = [
            constraint
            for constraint in constraints
            if not constraint["name"].startswith(tuple(args.omit_prefix))
        ]
    skipped_nonnegativity = 0
    for constraint in constraints:
        # LP variables are given explicit nonnegative bounds below.
        if constraint["name"].startswith(("q_nonneg_", "h_nonneg_")):
            skipped_nonnegativity += 1
            continue
        operator = {"eq": "=", "le": "<=", "ge": ">="}[constraint["sense"]]
        expression = linear_expression(constraint["coeffs"], order)
        # QSopt_ex's LP reader treats CPLEX-style row labels as columns, so rows
        # are deliberately emitted without labels.
        lines.append(f" {expression} {operator} {constraint['rhs']}")

    lines.append("Bounds")
    lines.extend(f" 0 <= {name}" for name in variables)
    lines.append("End")
    output.write_text("\n".join(lines) + "\n")

    print(f"Wrote {output}")
    print(f"{len(variables)} variables")
    print(f"{len(constraints) - skipped_nonnegativity} linear rows")
    print(
        f"Skipped {len(model['divisibility_constraints'])} non-LP divisibility constraints"
    )


if __name__ == "__main__":
    main()
