"""Solve a scaled floating relaxation of the exact BCH-sandwich model.

This script is diagnostic only.  A successful HiGHS result is useful for locating
candidate optima and assessing the strength of the model, but it is not a rigorous
bound.  Final claims require an independently verified exact rational certificate.
"""

from __future__ import annotations

import argparse
import json
from math import comb, ldexp, log2
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from spectrum_exact import symmetric_kraw_coeff


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = ROOT / "generated" / "coupled_lp_exact.json"


def variable_scale(name: str, n: int, q_dimension: int) -> float:
    """Return a binomial-scale estimate for one half-spectrum coefficient."""
    weight = int(name.split("_", 1)[1])
    # Q and a Q coset both contain 2^q_dimension words in the even-weight space.
    estimate = ldexp(float(comb(n, weight)), q_dimension - (n - 1))
    return max(1.0, estimate)


def scaled_matrix(rows, variables, scales):
    """Build a row-normalized sparse matrix and right-hand side."""
    index = {name: i for i, name in enumerate(variables)}
    data = []
    row_indices = []
    column_indices = []
    rhs_values = []

    for row_number, constraint in enumerate(rows):
        entries = []
        for name, coefficient_text in constraint["coeffs"].items():
            column = index[name]
            coefficient = float(int(coefficient_text)) * scales[column]
            entries.append((column, coefficient))

        rhs = float(int(constraint["rhs"]))
        normalization = max(
            1.0,
            abs(rhs),
            max((abs(value) for _, value in entries), default=0.0),
        )
        for column, value in entries:
            row_indices.append(row_number)
            column_indices.append(column)
            data.append(value / normalization)
        rhs_values.append(rhs / normalization)

    matrix = coo_matrix(
        (data, (row_indices, column_indices)),
        shape=(len(rows), len(variables)),
    ).tocsr()
    return matrix, np.asarray(rhs_values)


def build_relaxation(model, scale_overrides=None):
    metadata = model["metadata"]
    variables = metadata["variables"]
    scales = np.asarray([
        variable_scale(name, metadata["n"], metadata["Q_dimension"])
        for name in variables
    ])
    if scale_overrides:
        index = {name: i for i, name in enumerate(variables)}
        for name, scale in scale_overrides.items():
            if scale <= 0 or not np.isfinite(scale):
                raise ValueError(f"invalid scale for {name}: {scale}")
            scales[index[name]] = scale

    equalities = []
    inequalities = []
    for constraint in model["constraints"]:
        sense = constraint["sense"]
        if sense == "eq":
            equalities.append(constraint)
        elif sense == "le":
            inequalities.append(constraint)
        elif sense == "ge":
            inequalities.append({
                **constraint,
                "coeffs": {
                    name: str(-int(value))
                    for name, value in constraint["coeffs"].items()
                },
                "rhs": str(-int(constraint["rhs"])),
            })
        else:
            raise ValueError(f"unknown constraint sense: {sense}")

    a_eq, b_eq = scaled_matrix(equalities, variables, scales)
    a_ub, b_ub = scaled_matrix(inequalities, variables, scales)
    return variables, scales, a_eq, b_eq, a_ub, b_ub


def use_krawtchouk_oa_basis(model):
    """Replace the exact OA moment rows by an equivalent orthogonal basis."""
    metadata = model["metadata"]
    n = metadata["n"]
    q_dimension = metadata["Q_dimension"]
    weights = [int(name.split("_", 1)[1]) for name in metadata["variables"] if name.startswith("q_")]

    constraints = [
        constraint
        for constraint in model["constraints"]
        if not constraint["name"].startswith(("Q_OA_t", "H_OA_t"))
    ]
    for label, prefix in (("Q", "q_"), ("H", "h_")):
        for degree in range(16):
            coefficients = {
                f"{prefix}{weight}": str(symmetric_kraw_coeff(n, degree, weight))
                for weight in weights
                if symmetric_kraw_coeff(n, degree, weight)
            }
            constraints.append({
                "name": f"{label}_OA_Kraw_{degree}",
                "coeffs": coefficients,
                "sense": "eq",
                "rhs": str(1 << q_dimension if degree == 0 else 0),
            })
    return {**model, "constraints": constraints}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--objective", default="h_38")
    parser.add_argument("--sense", choices=("min", "max"), default="max")
    parser.add_argument(
        "--oa-basis",
        choices=("krawtchouk", "moment"),
        default="krawtchouk",
        help="floating representation of the exact OA row space",
    )
    parser.add_argument(
        "--objective-scale-log2",
        type=float,
        help="override the objective variable's column scale with 2^VALUE",
    )
    args = parser.parse_args()

    model = json.loads(args.model.read_text())
    if args.oa_basis == "krawtchouk":
        model = use_krawtchouk_oa_basis(model)
    scale_overrides = None
    if args.objective_scale_log2 is not None:
        scale_overrides = {args.objective: 2.0 ** args.objective_scale_log2}
    variables, scales, a_eq, b_eq, a_ub, b_ub = build_relaxation(
        model, scale_overrides
    )
    try:
        objective_index = variables.index(args.objective)
    except ValueError as error:
        raise SystemExit(f"unknown objective variable: {args.objective}") from error

    objective = np.zeros(len(variables))
    objective[objective_index] = 1.0 if args.sense == "min" else -1.0
    result = linprog(
        objective,
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=a_eq,
        b_eq=b_eq,
        bounds=(0.0, None),
        method="highs-ipm",
        # HiGHS presolve loses the low-shell information in this exceptionally
        # ill-conditioned model before the scaled IPM sees it.
        options={"presolve": False},
    )

    print("WARNING: floating-point diagnostic; this is not a rigorous bound")
    print(f"column scale for {args.objective} = 2^{log2(scales[objective_index]):.6f}")
    print(f"HiGHS status: {result.status} ({result.message})")
    if not result.success:
        raise SystemExit(1)

    value = result.x[objective_index] * scales[objective_index]
    print(f"{args.sense} {args.objective} = {value:.17g}")
    if value > 0:
        print(f"log2({args.objective}) = {log2(value):.12f}")
    print(f"max normalized equality residual = {np.max(np.abs(result.eqlin.residual)):.3e}")
    print(f"min normalized inequality slack = {np.min(result.ineqlin.residual):.3e}")


if __name__ == "__main__":
    main()
