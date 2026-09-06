#!/usr/bin/env python3
"""Verify the symmetric-function form of the 24-point incidence equations."""

from __future__ import annotations

import json
import random
from pathlib import Path

from affine_wambach import MODULUS, gf_mul
from analyze_locator_incidence import multiplication_table, polynomial_product_roots


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "generated" / "symmetric_incidence_identity.json"


def complete_homogeneous(points: list[int], maximum: int, mul: list[bytes]) -> list[int]:
    coefficients = [1] + [0] * maximum
    for point in points:
        row = mul[point]
        for degree in range(1, maximum + 1):
            coefficients[degree] ^= row[coefficients[degree - 1]]
    return coefficients


def monomial_remainder(exponent: int, modulus: list[int], mul: list[bytes]) -> list[int]:
    degree = len(modulus) - 1
    remainder = [1] + [0] * (degree - 1)
    for _ in range(exponent):
        shifted = [0, *remainder]
        leading = shifted[-1]
        if leading:
            row = mul[leading]
            for index, coefficient in enumerate(modulus):
                shifted[index] ^= row[coefficient]
        assert shifted[-1] == 0
        remainder = shifted[:-1]
    return remainder


def check_set(points: list[int], mul: list[bytes]) -> bool:
    assert len(points) == 24 and len(set(points)) == 24 and all(points)
    modulus = polynomial_product_roots(points, mul)
    elementary = [modulus[24 - index] for index in range(25)]
    homogeneous = complete_homogeneous(points, 127, mul)
    remainder = monomial_remainder(146, modulus, mul)

    assert remainder[0] == mul[elementary[24]][homogeneous[122]]
    for remainder_degree in range(19, 24):
        predicted = 0
        for elementary_degree in range(24 - remainder_degree):
            predicted ^= mul[elementary[elementary_degree]][
                homogeneous[146 - remainder_degree - elementary_degree]
            ]
        assert remainder[remainder_degree] == predicted

    # In characteristic two, H(t)=E(t)H(t)^2.  This halves the indices in
    # the six high-degree symmetric coordinates.
    for degree in range(122, 128):
        halved = 0
        for elementary_degree in range(degree & 1, 25, 2):
            if elementary_degree > degree:
                break
            lower = homogeneous[(degree - elementary_degree) // 2]
            halved ^= mul[elementary[elementary_degree]][mul[lower][lower]]
        assert homogeneous[degree] == halved

    remainder_condition = remainder[0] == 1 and all(
        remainder[degree] == 0 for degree in range(19, 24)
    )
    symmetric_condition = (
        mul[elementary[24]][homogeneous[122]] == 1
        and all(homogeneous[degree] == 0 for degree in range(123, 128))
    )
    assert remainder_condition == symmetric_condition
    return symmetric_condition


def main() -> None:
    mul = multiplication_table()
    locator = json.loads(
        (ROOT / "generated" / "locator_reduction_wambach.json").read_text()
    )
    agreements = [int(value, 16) for value in locator["agreement_points_hex"]]

    known_subsets = []
    for offset in range(len(agreements)):
        rotated = agreements[offset:] + agreements[:offset]
        known_subsets.append(rotated[:24])
    assert all(check_set(points, mul) for points in known_subsets)

    rng = random.Random(0x14624)
    random_trials = 1000
    random_hits = 0
    domain = list(range(1, 256))
    for _ in range(random_trials):
        random_hits += check_set(rng.sample(domain, 24), mul)

    result = {
        "classification": "exact finite-field identity verification",
        "field_modulus": hex(MODULUS),
        "set_size": 24,
        "target_monomial_exponent": 146,
        "remainder_condition": (
            "constant coefficient of x^146 mod P_T is 1 and coefficients "
            "19 through 23 vanish"
        ),
        "equivalent_symmetric_condition": (
            "e_24(T)*h_122(T)=1 and h_123(T)=...=h_127(T)=0"
        ),
        "coefficient_identity": (
            "[x^j](x^146 mod P_T)="
            "sum_{r=0}^{23-j} e_r(T) h_{146-j-r}(T) for 19<=j<=23"
        ),
        "frobenius_halving_identity": (
            "h_n=sum_{0<=i<=24, i congruent n mod 2} e_i*h_((n-i)/2)^2"
        ),
        "known_locator_subsets_checked": len(known_subsets),
        "random_sets_checked": random_trials,
        "random_admissible_sets_seen": random_hits,
        "scope": (
            "The identities are algebraic; the finite checks guard the field "
            "and coefficient conventions. They do not count their common zeros."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
