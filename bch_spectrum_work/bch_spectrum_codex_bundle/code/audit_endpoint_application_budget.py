"""Audit what the weight-38 test does and does not close downstream.

Use the directed transfer logarithms with Arb enclosures, then sum exactly.
The historical A38 target is an allocation from a joint low-shell envelope;
it is not by itself a sufficient bound on the full Q1 functional.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

from certify_random_inner_threshold_outward import pow2_upper, fraction_log2

ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "generated"


def exact(row):
    return Fraction(int(row["numerator"]), int(row["denominator"]))


def record(x):
    return {"numerator": str(x.numerator), "denominator": str(x.denominator),
            "log2_diagnostic": fraction_log2(x)}


def main():
    path = GENERATED / "random_inner_threshold_outward.json"
    lp_path = GENERATED / "exact_bch_lp_caps.json"
    threshold = json.loads(path.read_text())
    lp = json.loads(lp_path.read_text())
    target = Fraction(threshold["parameters"]["target"])
    high = exact(threshold["high_functional_upper"])
    logs = threshold["coefficient_log2_upper"]
    coefficient = {int(w): Fraction.from_float(pow2_upper(x)) for w, x in logs.items()}
    rows = []
    total = high
    without38 = high
    for row in lp["shells"]:
        w = row["weight"]
        cap = (int(threshold["a38_sufficient_cap"]["certified_integer_cap"])
               if w == 38 else row["exact_Aw_C_cap"])
        pair = coefficient[w] + coefficient[256-w]
        term = cap*pair
        total += term
        if w != 38:
            without38 += term
        rows.append({"weight": w, "multiplicity_cap": str(cap),
                     "cap_source": "conditional endpoint test target" if w == 38 else "exact BCH-sandwich LP",
                     "pair_coefficient_upper": record(pair), "functional_term_upper": record(term)})
    # These assertions certify failure of this sufficient upper-bound check,
    # not failure of the code's true distance or true Q1 functional.
    assert total > target and without38 > target
    proposed_caps = {40: 50_000_000_000_000, 42: 5_000_000_000_000_000}
    proposed_total = high
    for row in rows:
        cap = proposed_caps.get(row["weight"], int(row["multiplicity_cap"]))
        proposed_total += cap*exact(row["pair_coefficient_upper"])
    assert proposed_total < target
    result = {
        "classification": "exact rational aggregation of directed transfer coefficients; audit of insufficient combined upper bounds",
        "a38_bound_alone_closes_Q1": False,
        "reason": "The historical threshold allocates the residual over a joint inflated low-shell model. Its A38 allocation is not a standalone Q1 certificate with the current other-shell LP caps.",
        "target": str(target), "rows": rows, "high_shell_upper": record(high),
        "combined_upper_with_endpoint_cap": record(total),
        "other_shell_upper_even_if_A38_were_zero": record(without38),
        "combined_margin_bits_diagnostic": -fraction_log2(total),
        "other_shell_margin_bits_diagnostic": -fraction_log2(without38),
        "proposed_two_shell_targets": {
            "classification": "unproved multiplicity targets; conditional arithmetic implication only",
            "caps": {str(k): str(v) for k, v in proposed_caps.items()},
            "combined_upper_if_targets_hold": record(proposed_total),
            "margin_bits_diagnostic": -fraction_log2(proposed_total),
            "sufficient_for_audited_interior_shell_sum": True},
        "needed_next": "Improve other low-shell caps, led by weight 40; reallocate and outward-check the full shell budget before claiming Q1 closure.",
        "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (path, lp_path, Path(__file__), ROOT / "code" / "certify_random_inner_threshold_outward.py")}
    }
    output = GENERATED / "endpoint_application_budget_audit.json"
    output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"combined_margin_bits": result["combined_margin_bits_diagnostic"],
                      "other_shell_margin_bits": result["other_shell_margin_bits_diagnostic"],
                      "proposed_40_42_target_margin_bits": -fraction_log2(proposed_total),
                      "a38_alone_closes_Q1": False, "receipt": str(output)}, indent=2))


if __name__ == "__main__":
    main()
