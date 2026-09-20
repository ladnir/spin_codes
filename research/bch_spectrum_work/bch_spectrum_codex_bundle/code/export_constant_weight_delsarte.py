#!/usr/bin/env python3
"""Export the Johnson-scheme Delsarte LP for a constant-weight shell.

For a binary constant-weight code of length n, weight w, and minimum distance
2*delta, let a_j denote its average number of neighbors at distance 2*j.
The LP maximizes 1 + sum(a_j) subject to the exact MacWilliams positivity
constraints of the Johnson association scheme.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def eberlein(n: int, w: int, relation: int, eigenspace: int) -> int:
    """Return P[eigenspace, relation] for the Johnson scheme J(n,w)."""

    result = 0
    for index in range(relation + 1):
        if index > w - eigenspace:
            continue
        if n - w + index - eigenspace < index:
            continue
        result += (
            (-1) ** (relation - index)
            * math.comb(w - index, relation - index)
            * math.comb(w - eigenspace, index)
            * math.comb(n - w + index - eigenspace, index)
        )
    return result


def relation_valency(n: int, w: int, relation: int) -> int:
    return math.comb(w, relation) * math.comb(n - w, relation)


def coefficient(n: int, w: int, relation: int, eigenspace: int) -> Fraction:
    """Return Q[relation,eigenspace] divided by its positive multiplicity."""

    return Fraction(
        eberlein(n, w, relation, eigenspace),
        relation_valency(n, w, relation),
    )


def term_text(value: Fraction, variable: str, first: bool) -> str:
    sign = "-" if value < 0 else "+"
    magnitude = abs(value)
    if magnitude.denominator == 1:
        scalar = str(magnitude.numerator)
    else:
        scalar = f"{magnitude.numerator}/{magnitude.denominator}"
    body = variable if magnitude == 1 else f"{scalar} {variable}"
    if first:
        return body if value >= 0 else f"- {body}"
    return f" {sign} {body}"


def expression(terms: list[tuple[Fraction, str]]) -> str:
    nonzero = [(value, variable) for value, variable in terms if value]
    if not nonzero:
        return "0"
    return "".join(
        term_text(value, variable, position == 0)
        for position, (value, variable) in enumerate(nonzero)
    )


def export_lp(n: int, w: int, delta: int, output: Path) -> dict:
    rank = min(w, n - w)
    relations = list(range(delta, rank + 1))
    variables = [f"a_{relation}" for relation in relations]
    lines = [
        "Maximize",
        " obj: " + " + ".join(variables),
        "Subject To",
    ]
    for eigenspace in range(1, rank + 1):
        terms = [
            (coefficient(n, w, relation, eigenspace), f"a_{relation}")
            for relation in relations
        ]
        # The omitted relation-zero term is a_0=1 and has coefficient one.
        lines.append(f" pos_{eigenspace}: {expression(terms)} >= -1")

    lines.append("Bounds")
    for relation in relations:
        upper = relation_valency(n, w, relation)
        lines.append(f" 0 <= a_{relation} <= {upper}")
    lines.append("End")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="ascii")

    packing_upper = math.comb(n, w - delta + 1) // math.comb(
        w, w - delta + 1
    )
    metadata = {
        "classification": "exact rational Johnson-scheme Delsarte LP",
        "n": n,
        "weight": w,
        "minimum_distance": 2 * delta,
        "relations": relations,
        "variables": len(variables),
        "positivity_constraints": rank,
        "objective_value_excludes_a0": True,
        "shell_bound_is_one_plus_objective": True,
        "simple_packing_upper": packing_upper,
        "simple_packing_upper_log2": math.log2(packing_upper),
        "lp_path": str(output),
    }
    metadata_path = output.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def attach_qsopt_certificate(
    metadata: dict, n: int, w: int, delta: int, solution: Path
) -> dict:
    """Parse and independently check the exact primal solution from QSopt_ex."""

    lines = solution.read_text(encoding="ascii").splitlines()
    if not any(line.strip() in {"status = OPTIMAL", "status OPTIMAL"} for line in lines):
        raise ValueError("QSopt_ex solution is not marked OPTIMAL")
    objective_line = next(
        line for line in lines if line.strip().startswith("Value =")
    )
    objective = Fraction(objective_line.split("=", 1)[1].strip())

    variables: dict[int, Fraction] = {}
    in_variables = False
    for line in lines:
        stripped = line.strip()
        if stripped == "VARS:":
            in_variables = True
            continue
        if in_variables and stripped.endswith(":"):
            break
        if in_variables and stripped:
            name, value = stripped.split("=", 1)
            variables[int(name.strip().removeprefix("a_"))] = Fraction(value.strip())

    rank = min(w, n - w)
    relations = range(delta, rank + 1)
    values = {relation: variables.get(relation, Fraction(0)) for relation in relations}
    objective_matches_variables = sum(values.values(), Fraction(0)) == objective
    bounds_hold = all(
        0 <= values[relation] <= relation_valency(n, w, relation)
        for relation in relations
    )
    positivity_slacks = [
        Fraction(1)
        + sum(
            coefficient(n, w, relation, eigenspace) * values[relation]
            for relation in relations
        )
        for eigenspace in range(1, rank + 1)
    ]
    positivity_holds = all(slack >= 0 for slack in positivity_slacks)
    if not (objective_matches_variables and bounds_hold and positivity_holds):
        raise ValueError("the reported exact primal solution failed independent checks")

    # For M a >= -1, write -M a <= 1.  A nonnegative vector y with
    # (-M)^T y >= 1 gives the upper bound sum(y) on sum(a_j).  QSopt_ex reports
    # the row prices with the opposite sign for these >= constraints.
    row_prices: dict[int, Fraction] = {}
    in_prices = False
    for line in lines:
        stripped = line.strip()
        if stripped == "PI:":
            in_prices = True
            continue
        if in_prices and stripped.endswith(":"):
            break
        if in_prices and stripped:
            name, value = stripped.split("=", 1)
            row_prices[int(name.strip().removeprefix("pos_"))] = Fraction(
                value.strip()
            )
    dual = {
        eigenspace: -row_prices.get(eigenspace, Fraction(0))
        for eigenspace in range(1, rank + 1)
    }
    dual_nonnegative = all(value >= 0 for value in dual.values())
    dual_covers_variables = all(
        sum(
            -coefficient(n, w, relation, eigenspace) * dual[eigenspace]
            for eigenspace in range(1, rank + 1)
        )
        >= 1
        for relation in relations
    )
    dual_objective_matches = sum(dual.values(), Fraction(0)) == objective
    if not (dual_nonnegative and dual_covers_variables and dual_objective_matches):
        raise ValueError("the reported exact dual solution failed independent checks")

    shell_upper = objective + 1
    integer_upper = shell_upper.numerator // shell_upper.denominator
    metadata.update(
        {
            "qsopt_solution_path": str(solution.resolve()),
            "qsopt_status": "OPTIMAL",
            "exact_objective": str(objective),
            "exact_shell_upper": str(shell_upper),
            "integer_shell_upper": integer_upper,
            "integer_shell_upper_log2": math.log2(integer_upper),
            "independent_exact_primal_checks": {
                "objective_matches_variables": objective_matches_variables,
                "variable_bounds_hold": bounds_hold,
                "johnson_positivity_holds": positivity_holds,
            },
            "independent_exact_dual_checks": {
                "multipliers_nonnegative": dual_nonnegative,
                "all_distance_variables_covered": dual_covers_variables,
                "dual_objective_matches_primal": dual_objective_matches,
            },
            "improvement_over_simple_packing_bits": (
                metadata["simple_packing_upper_log2"] - math.log2(integer_upper)
            ),
        }
    )
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=256)
    parser.add_argument("--weight", type=int, default=52)
    parser.add_argument("--minimum-distance", type=int, default=38)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "generated" / "johnson_n256_w52_d38.lp",
    )
    parser.add_argument(
        "--solution",
        type=Path,
        help="optional exact QSopt_ex .sol file to check and attach to the receipt",
    )
    args = parser.parse_args()
    if args.minimum_distance % 2:
        raise SystemExit("minimum distance must be even")
    delta = args.minimum_distance // 2
    if not (0 < delta <= min(args.weight, args.n - args.weight)):
        raise SystemExit("invalid Johnson-scheme parameters")
    metadata = export_lp(args.n, args.weight, delta, args.output)
    solution = args.solution
    if solution is None and args.output.with_suffix(".sol").is_file():
        solution = args.output.with_suffix(".sol")
    if solution is not None:
        metadata = attach_qsopt_certificate(
            metadata, args.n, args.weight, delta, solution
        )
        args.output.with_suffix(".json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
