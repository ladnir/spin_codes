#!/usr/bin/env python3
"""Optimize the 24th agreement moment from the exact RS coset moments."""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

import z3


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "generated" / "rs_coset_moment_lp.json"
Q = 256
EVALUATION_POINTS = 255
DIMENSION = 18
MAXIMUM_AGREEMENTS = 37
TARGET_MOMENT = 24


def as_fraction(value: z3.ArithRef) -> Fraction:
    if z3.is_int_value(value):
        return Fraction(value.as_long(), 1)
    if z3.is_rational_value(value):
        return Fraction(value.numerator_as_long(), value.denominator_as_long())
    raise TypeError(f"expected exact rational, got {value!r}")


def log2_fraction(value: Fraction) -> float:
    return math.log2(value.numerator) - math.log2(value.denominator)


def main() -> None:
    multiplicities = [z3.Real(f"B_{agreements}") for agreements in range(38)]
    optimizer = z3.Optimize()
    for multiplicity in multiplicities:
        optimizer.add(multiplicity >= 0)

    moment_rhs: list[int] = []
    for order in range(DIMENSION + 1):
        rhs = math.comb(EVALUATION_POINTS, order) * Q ** (DIMENSION - order)
        moment_rhs.append(rhs)
        optimizer.add(
            z3.Sum(
                [
                    math.comb(agreements, order) * multiplicities[agreements]
                    for agreements in range(order, MAXIMUM_AGREEMENTS + 1)
                ]
            )
            == rhs
        )

    objective_expression = z3.Sum(
        [
            math.comb(agreements, TARGET_MOMENT) * multiplicities[agreements]
            for agreements in range(TARGET_MOMENT, MAXIMUM_AGREEMENTS + 1)
        ]
    )
    handle = optimizer.maximize(objective_expression)
    status = optimizer.check()
    if status != z3.sat:
        raise RuntimeError(f"exact moment LP status: {status}")
    model = optimizer.model()
    objective = as_fraction(model.eval(objective_expression))
    assert as_fraction(handle.value()) == objective

    solution: dict[int, Fraction] = {}
    for agreements, variable in enumerate(multiplicities):
        value = as_fraction(model.eval(variable, model_completion=True))
        if value:
            solution[agreements] = value

    for order, rhs in enumerate(moment_rhs):
        lhs = sum(
            value * math.comb(agreements, order)
            for agreements, value in solution.items()
        )
        assert lhs == rhs

    dual_variables = [z3.Real(f"y_{order}") for order in range(DIMENSION + 1)]
    dual_optimizer = z3.Optimize()
    for agreements in range(MAXIMUM_AGREEMENTS + 1):
        dual_optimizer.add(
            z3.Sum(
                [
                    math.comb(agreements, order) * dual_variables[order]
                    for order in range(min(agreements, DIMENSION) + 1)
                ]
            )
            >= math.comb(agreements, TARGET_MOMENT)
        )
    dual_expression = z3.Sum(
        [moment_rhs[order] * dual_variables[order] for order in range(DIMENSION + 1)]
    )
    dual_handle = dual_optimizer.minimize(dual_expression)
    dual_status = dual_optimizer.check()
    if dual_status != z3.sat:
        raise RuntimeError(f"exact dual moment LP status: {dual_status}")
    dual_model = dual_optimizer.model()
    dual_objective = as_fraction(dual_model.eval(dual_expression))
    assert as_fraction(dual_handle.value()) == dual_objective == objective
    dual_solution = [
        as_fraction(dual_model.eval(variable, model_completion=True))
        for variable in dual_variables
    ]
    for agreements in range(MAXIMUM_AGREEMENTS + 1):
        majorant = sum(
            dual_solution[order] * math.comb(agreements, order)
            for order in range(min(agreements, DIMENSION) + 1)
        )
        assert majorant >= math.comb(agreements, TARGET_MOMENT)

    random_reference = Fraction(math.comb(255, 24), Q**6)
    ratio = objective / random_reference
    locator_upper = objective / math.comb(37, 24)
    result = {
        "classification": "exact rational linear-program optimum",
        "model": {
            "variables": "B_n for 0 <= n <= 37",
            "meaning": "number of degree-at-most-18 fixed-constant polynomials with n agreements",
            "constraints": (
                "B_n >= 0 and sum_n B_n*C(n,s)=C(255,s)*256^(18-s) "
                "for 0 <= s <= 18"
            ),
            "objective": "maximize sum_n B_n*C(n,24)",
        },
        "objective_I24_upper": {
            "numerator": str(objective.numerator),
            "denominator": str(objective.denominator),
            "log2": log2_fraction(objective),
        },
        "incidence_ratio_R6_upper": {
            "numerator": str(ratio.numerator),
            "denominator": str(ratio.denominator),
            "decimal": float(ratio),
            "log2": log2_fraction(ratio),
        },
        "locator_count_upper": {
            "numerator": str(locator_upper.numerator),
            "denominator": str(locator_upper.denominator),
            "log2": log2_fraction(locator_upper),
        },
        "nonzero_extremizer": {
            str(agreements): {
                "numerator": str(value.numerator),
                "denominator": str(value.denominator),
                "log2": log2_fraction(value),
            }
            for agreements, value in solution.items()
        },
        "dual_binomial_basis_coefficients": {
            str(order): {
                "numerator": str(value.numerator),
                "denominator": str(value.denominator),
            }
            for order, value in enumerate(dual_solution)
            if value
        },
        "exact_certificate_checks": {
            "all_primal_moment_equalities": True,
            "all_primal_variables_nonnegative": True,
            "all_38_dual_majorant_inequalities": True,
            "primal_dual_objectives_equal": True,
        },
        "interpretation": (
            "These eighteen universal MDS/coset moments do not use the special "
            "monomial beyond the proven 37-agreement endpoint. The extremizer "
            "identifies which additional algebraic constraints are needed."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "R6_upper": float(ratio),
        "R6_upper_log2": log2_fraction(ratio),
        "locator_upper_log2": log2_fraction(locator_upper),
        "extremizer_support": sorted(solution),
    }, indent=2))


if __name__ == "__main__":
    main()
