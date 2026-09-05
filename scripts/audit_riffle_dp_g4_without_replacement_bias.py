#!/usr/bin/env python3
"""Record the exact counterexample to the retired conditional-bias product."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = ROOT / "explorations" / "riffle_dp_g4_g2_component_mixing_probe.json"
OUTPUT = ROOT / "explorations" / "riffle_dp_g4_without_replacement_bias_audit.json"
M = 524_352


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    component = json.loads(COMPONENT.read_text())
    degree_one = next(
        row for row in component["component_rows"]
        if row["component_degrees"] == [1]
    )
    value_four = next(
        row for row in degree_one["value_rows"] if row["packet_value"] == 4
    )
    if value_four["maximum_absolute_bias_numerator"] != 0:
        raise RuntimeError("bias audit: degree-one value-four character is not balanced")

    exact = Fraction(1, M - 1)
    retired = Fraction(1, (M - 1) ** 2)
    payload = {
        "schema": "riffle-dp-g4-without-replacement-bias-audit-v1",
        "evidence_label": "EXACT_REFUTATION_OF_BOUNDING_STEP",
        "component_receipt_sha256": digest(COMPONENT),
        "packet_positions": M,
        "component_degree": 1,
        "packet_value": 4,
        "character_sum": 0,
        "exact_two_draw_character_expectation": {
            "sign": -1,
            "magnitude_numerator": str(exact.numerator),
            "magnitude_denominator": str(exact.denominator),
        },
        "retired_product_magnitude": {
            "numerator": str(retired.numerator),
            "denominator": str(retired.denominator),
        },
        "underestimate_factor": M - 1,
        "derivation": (
            "For signs s_p in {+1,-1} with sum_p s_p=0, the ordered sum over "
            "distinct p,q of s_p s_q equals -M. Division by M(M-1) gives "
            "-1/(M-1). The retired argument multiplied two marginal caps instead."
        ),
        "replacement": (
            "Sample labeled cells independently and uniformly, then condition on "
            "all cells being distinct. The conditioned law is the uniform ordered "
            "injection, so P_perm(E) <= P_iid(E)/P_iid(distinct)."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"exact_magnitude=1/{M - 1}")
    print(f"retired_magnitude=1/{(M - 1) ** 2}")
    print(f"underestimate_factor={M - 1}")
    print(f"output={OUTPUT}")
    print("status=EXACT_REFUTATION_OF_BOUNDING_STEP")


if __name__ == "__main__":
    main()
