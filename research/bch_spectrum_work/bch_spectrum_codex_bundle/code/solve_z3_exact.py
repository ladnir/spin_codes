"""Check or optimize the BCH-sandwich LP with exact Z3 rational arithmetic."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import z3

from export_lp import independent_oa_constraints


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = ROOT / "generated" / "coupled_lp_exact.json"


def fixed_variables(constraints):
    """Extract variables fixed by one-term equality rows."""
    fixed = {}
    for constraint in constraints:
        if constraint["sense"] != "eq" or len(constraint["coeffs"]) != 1:
            continue
        name, coefficient_text = next(iter(constraint["coeffs"].items()))
        coefficient = int(coefficient_text)
        rhs = int(constraint["rhs"])
        value, remainder = divmod(rhs, coefficient)
        if remainder == 0:
            if name in fixed and fixed[name] != value:
                raise ArithmeticError(f"conflicting fixed values for {name}")
            fixed[name] = value
    return fixed


def z3_expression(coefficients, variables, fixed):
    terms = []
    for name, coefficient_text in coefficients.items():
        coefficient = int(coefficient_text)
        if name in fixed:
            terms.append(z3.IntVal(coefficient * fixed[name]))
        else:
            terms.append(z3.IntVal(coefficient) * variables[name])
    return z3.Sum(terms) if terms else z3.IntVal(0)


def add_continuous_constraints(solver, constraints, variables, fixed) -> int:
    for constraint in constraints:
        expression = z3_expression(constraint["coeffs"], variables, fixed)
        rhs = z3.IntVal(int(constraint["rhs"]))
        sense = constraint["sense"]
        if sense == "eq":
            solver.add(expression == rhs)
        elif sense == "le":
            solver.add(expression <= rhs)
        elif sense == "ge":
            solver.add(expression >= rhs)
        else:
            raise ValueError(f"unknown constraint sense: {sense}")
    return len(constraints)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--mode", choices=("feasible", "optimize"), default="feasible")
    parser.add_argument("--objective", default="h_38")
    parser.add_argument("--sense", choices=("min", "max"), default="max")
    parser.add_argument("--timeout-seconds", type=int, default=600)
    parser.add_argument(
        "--omit-prefix",
        action="append",
        default=[],
        help="omit constraints whose names start with this prefix (diagnostics)",
    )
    args = parser.parse_args()

    model = json.loads(args.model.read_text())
    names = model["metadata"]["variables"]
    if args.objective not in names:
        raise SystemExit(f"unknown objective variable: {args.objective}")
    constraints = independent_oa_constraints(model)
    if args.omit_prefix:
        constraints = [
            constraint
            for constraint in constraints
            if not constraint["name"].startswith(tuple(args.omit_prefix))
        ]
    fixed = fixed_variables(constraints)
    variables = {name: z3.Real(name) for name in names if name not in fixed}
    z3.set_option(timeout=args.timeout_seconds * 1000)

    if args.mode == "feasible":
        solver = z3.Then("simplify", "solve-eqs", "qflra").solver()
        solver.set(timeout=args.timeout_seconds * 1000)
        row_count = add_continuous_constraints(solver, constraints, variables, fixed)
        result = solver.check()
        print(f"Exact Z3 feasibility status: {result}")
        print(
            f"{len(variables)} rational variables, {len(fixed)} fixed variables, "
            f"{row_count} linear rows"
        )
        if result == z3.unknown:
            print(f"Reason: {solver.reason_unknown()}")
            raise SystemExit(2)
        if result == z3.unsat:
            raise SystemExit(1)
        value = (
            z3.IntVal(fixed[args.objective])
            if args.objective in fixed
            else solver.model().eval(variables[args.objective], model_completion=True)
        )
        print(f"One feasible model has {args.objective} = {value}")
        return

    optimizer = z3.Optimize()
    optimizer.set(timeout=args.timeout_seconds * 1000)
    row_count = add_continuous_constraints(
        optimizer, constraints, variables, fixed
    )
    objective = (
        z3.IntVal(fixed[args.objective])
        if args.objective in fixed
        else variables[args.objective]
    )
    handle = optimizer.maximize(objective) if args.sense == "max" else optimizer.minimize(objective)
    result = optimizer.check()
    print(f"Exact Z3 optimization status: {result}")
    print(
        f"{len(variables)} rational variables, {len(fixed)} fixed variables, "
        f"{row_count} linear rows"
    )
    if result == z3.unknown:
        print(f"Reason: {optimizer.reason_unknown()}")
        raise SystemExit(2)
    if result == z3.unsat:
        raise SystemExit(1)
    bound = handle.upper() if args.sense == "max" else handle.lower()
    value = optimizer.model().eval(objective, model_completion=True)
    print(f"Exact {args.sense} bound for {args.objective} = {bound}")
    print(f"Objective value in returned model = {value}")


if __name__ == "__main__":
    main()
