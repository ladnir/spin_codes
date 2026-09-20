#!/usr/bin/env python3
"""Build an exact/sampled ladder for normalized BCH locator families.

For a punctured BCH code with designed distance d=2t+1, normalized locators
have the form z^d+g(z^2), where deg(g)<=t and g(0)=1.  The substitution
u=z^2 turns the roots into agreements g(u)=u^e.  This script enumerates every
t-point interpolant for q=8 and q=32 and imports Monte Carlo receipts for
q=128 and q=256.
"""

from __future__ import annotations

import itertools
import json
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path

from prepare_wambach_bases import gf_mul_generic, gf_pow_generic


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "generated"


def evaluate(coefficients: list[int], point: int, m: int, modulus: int) -> int:
    value = 0
    for coefficient in reversed(coefficients):
        value = gf_mul_generic(value, point, m, modulus) ^ coefficient
    return value


def interpolate(
    nodes: tuple[int, ...], values: list[int], m: int, modulus: int
) -> list[int]:
    divided = values[:]
    degree = len(nodes) - 1
    inverses = {
        value: gf_pow_generic(value, (1 << m) - 2, m, modulus)
        for value in range(1, 1 << m)
    }
    for order in range(1, degree + 1):
        for index in range(degree, order - 1, -1):
            divided[index] = gf_mul_generic(
                divided[index] ^ divided[index - 1],
                inverses[nodes[index] ^ nodes[index - order]],
                m,
                modulus,
            )

    polynomial = [divided[degree]]
    for index in range(degree - 1, -1, -1):
        updated = [0] * (len(polynomial) + 1)
        for j, coefficient in enumerate(polynomial):
            updated[j] ^= gf_mul_generic(coefficient, nodes[index], m, modulus)
            updated[j + 1] ^= coefficient
        updated[0] ^= divided[index]
        polynomial = updated
    return polynomial


def exact_family(m: int, modulus: int, t: int, d: int, expected_list: int) -> dict:
    q = 1 << m
    exponent = (d * (q // 2)) % (q - 1)
    target = [gf_pow_generic(u, exponent, m, modulus) for u in range(q)]
    histogram: Counter[int] = Counter()
    moment_sums = [0] * (d - t + 1)

    for base in itertools.combinations(range(1, q), t):
        nodes = (0,) + base
        values = [1] + [target[u] for u in base]
        polynomial = interpolate(nodes, values, m, modulus)
        agreements = sum(
            evaluate(polynomial, u, m, modulus) == target[u]
            for u in range(1, q)
        )
        assert t <= agreements <= d
        extra = agreements - t
        histogram[extra] += 1
        for order in range(d - t + 1):
            moment_sums[order] += comb(extra, order)

    bases = comb(q - 1, t)
    list_size, remainder = divmod(histogram[d - t], comb(d, t))
    assert remainder == 0 and list_size == expected_list
    moments = []
    for order, moment_sum in enumerate(moment_sums):
        reference = Fraction(comb(q - 1 - t, order), q**order)
        empirical = Fraction(moment_sum, bases)
        ratio = empirical / reference
        moments.append(
            {
                "order": order,
                "moment_exact": f"{empirical.numerator}/{empirical.denominator}",
                "binomial_reference_exact": (
                    f"{reference.numerator}/{reference.denominator}"
                ),
                "ratio_exact": f"{ratio.numerator}/{ratio.denominator}",
                "ratio": float(ratio),
            }
        )
    endpoint_reference = Fraction(comb(q - 1, d), q ** (d - t))
    endpoint_ratio = Fraction(list_size, 1) / endpoint_reference
    return {
        "classification": "exact enumeration of every base interpolant",
        "field_size": q,
        "field_modulus": hex(modulus),
        "base_points": t,
        "locator_degree": d,
        "target_exponent": exponent,
        "base_subsets": bases,
        "extra_agreement_histogram": dict(sorted(histogram.items())),
        "normalized_locator_list_size": list_size,
        "factorial_moment_comparison": moments,
        "endpoint_full_splitting_ratio": float(endpoint_ratio),
        "endpoint_full_splitting_ratio_exact": (
            f"{endpoint_ratio.numerator}/{endpoint_ratio.denominator}"
        ),
    }


def sampled_family(path: Path, exact_list_size: int | None = None) -> dict:
    receipt = json.loads(path.read_text())
    result = {
        "classification": receipt["classification"],
        "field_size": receipt["field_size"],
        "field_modulus": hex(receipt["field_modulus"]),
        "base_points": receipt["base_points_per_interpolant"],
        "locator_degree": receipt["agreement_count_upper"],
        "target_exponent": receipt["target_exponent"],
        "samples": receipt["samples"],
        "extra_agreement_histogram": receipt["extra_agreement_histogram"],
        "factorial_moment_comparison": receipt["binomial_moment_comparison"],
    }
    if exact_list_size is not None:
        q = result["field_size"]
        t = result["base_points"]
        d = result["locator_degree"]
        endpoint_reference = Fraction(comb(q - 1, d), q ** (d - t))
        endpoint_ratio = Fraction(exact_list_size, 1) / endpoint_reference
        result["normalized_locator_list_size"] = exact_list_size
        result["endpoint_full_splitting_ratio"] = float(endpoint_ratio)
        result["endpoint_full_splitting_ratio_exact"] = (
            f"{endpoint_ratio.numerator}/{endpoint_ratio.denominator}"
        )
        endpoint_incidence = []
        for order in range(min(8, d - t) + 1):
            contribution = Fraction(
                exact_list_size * comb(d, t + order) * q**order,
                comb(q - 1, t + order),
            )
            endpoint_incidence.append(
                {
                    "order": order,
                    "reference_ratio_contribution": float(contribution),
                    "reference_ratio_contribution_exact": (
                        f"{contribution.numerator}/{contribution.denominator}"
                    ),
                }
            )
        result["endpoint_incidence_contribution"] = endpoint_incidence
    if "middle_subfield" in receipt:
        result["middle_subfield"] = receipt["middle_subfield"]
    return result


def main() -> None:
    result = {
        "classification": {
            "q8_q32": "exact",
            "q128_q256_incidence_moments": "diagnostic Monte Carlo",
            "q128_endpoint_list": "exact from the published weight enumerator",
        },
        "families": [
            exact_family(3, 0xB, 1, 3, 1),
            exact_family(5, 0x25, 3, 7, 5),
            sampled_family(GENERATED / "locator_extension_monte_carlo_q128.json", 330),
            sampled_family(GENERATED / "locator_extension_monte_carlo.json"),
        ],
    }
    output = GENERATED / "locator_size_ladder.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
