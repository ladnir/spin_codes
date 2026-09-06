#!/usr/bin/env python3
"""Verify the normalized locator reduction for BCH(255,131,37).

For a weight-37 support, the BCH syndromes through degree 36 vanish.  Newton
identities then put its normalized locator in the form

    Lambda(z) = g(z^2) + z^37,  deg(g) <= 18,  g(0) = 1.

Writing u=z^2 and using z=u^128 in GF(256), the roots correspond exactly to
the agreements g(u)=u^146 on GF(256)^*.  The script checks every identity on
Wambach's published minimum word and emits a compact reproducibility receipt.
"""

from __future__ import annotations

import json
from pathlib import Path

from affine_wambach import ALPHA, P_WORD_OCTAL, gf_mul, gf_pow, octal_bits


ROOT = Path(__file__).resolve().parents[1]
BETA_EXPONENT = 13


def locator_coefficients(support: list[int]) -> list[int]:
    """Ascending coefficients of product_{x in support} (1 + x*z)."""
    coefficients = [1]
    for x in support:
        updated = [0] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            updated[degree] ^= coefficient
            updated[degree + 1] ^= gf_mul(coefficient, x)
        coefficients = updated
    return coefficients


def evaluate(coefficients: list[int], x: int) -> int:
    value = 0
    for coefficient in reversed(coefficients):
        value = gf_mul(value, x) ^ coefficient
    return value


def wambach_support() -> list[int]:
    beta = gf_pow(ALPHA, BETA_EXPONENT)
    coordinate = 1
    support = []
    for bit in octal_bits(P_WORD_OCTAL):
        if bit:
            support.append(coordinate)
        coordinate = gf_mul(coordinate, beta)
    assert len(support) == 37
    return support


def verify_locator_reduction() -> dict:
    beta = gf_pow(ALPHA, BETA_EXPONENT)
    support = wambach_support()

    syndromes = []
    for degree in range(1, 37):
        syndrome = 0
        for x in support:
            syndrome ^= gf_pow(x, degree)
        syndromes.append(syndrome)
    assert syndromes == [0] * 36

    original_locator = locator_coefficients(support)
    leading = original_locator[37]
    shifts = [
        shift
        for shift in range(255)
        if gf_mul(leading, gf_pow(beta, 37 * shift)) == 1
    ]
    assert len(shifts) == 1
    shift = shifts[0]
    normalized_support = [gf_mul(x, gf_pow(beta, shift)) for x in support]
    locator = locator_coefficients(normalized_support)
    assert locator[0] == locator[37] == 1
    assert all(locator[degree] == 0 for degree in range(1, 37, 2))

    g = [locator[2 * degree] for degree in range(19)]
    assert g[0] == 1
    agreements = [
        u
        for u in range(1, 256)
        if evaluate(g, u) == gf_pow(u, 146)
    ]
    expected_agreements = sorted(gf_pow(x, 253) for x in normalized_support)
    assert sorted(agreements) == expected_agreements
    assert len(agreements) == 37

    result = {
        "classification": "exact finite-field verification of locator reduction",
        "field_modulus": "0x14d",
        "coordinate_primitive_element": "alpha^13",
        "punctured_code": "BCH(255,131,37)",
        "syndromes_1_through_36_zero": True,
        "unique_cyclic_normalizing_shift": shift,
        "normalized_locator": "Lambda(z)=g(z^2)+z^37",
        "g_degree_upper": 18,
        "g_constant": g[0],
        "g_coefficients_hex_ascending": [f"{value:02x}" for value in g],
        "agreement_equation": "g(u)=u^146 on GF(256)^*",
        "agreement_count": len(agreements),
        "agreement_points_hex": [f"{value:02x}" for value in agreements],
        "counting_identity": (
            "Let L_37 be the normalized-polynomial list size. Then "
            "L_37=A_37(P_punctured)/255=19*h_38/128. Indeed, the affine "
            "2-design gives A_37(P_punctured)=38*A_38(P)/256, while "
            "A_38(P)=255*h_38. Each weight-37 punctured support has a full "
            "cyclic orbit because no nontrivial divisor of 255 divides 37."
        ),
        "application_relation": "A_38(C)=3968*L_37/19",
    }
    output = ROOT / "generated" / "locator_reduction_wambach.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    result = verify_locator_reduction()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
