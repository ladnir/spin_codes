#!/usr/bin/env python3
"""Combine exact BCH-sandwich shell caps with the RandomStepConv transfer."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "generated"


def logadd(values: list[float]) -> float:
    largest = max(values)
    return largest + math.log(sum(math.exp(value - largest) for value in values))


def main() -> None:
    diagnostic = json.loads(
        (GENERATED / "random_inner_oa15_majorant_diagnostic.json").read_text()
    )
    log_coefficients = diagnostic["log_shell_coefficients_natural"]
    application_caps = {
        row["weight"]: row["sufficient_log2_multiplicity_cap"]
        for row in diagnostic["proofish_low_prefix_reduction"][
            "sufficient_common_multiplier_shell_caps"
        ]
    }

    h38 = json.loads(
        (GENERATED / "max_h_38_scaled_rational.certificate.json").read_text()
    )
    shell_caps = {38: int(h38["derived_A38_C_cap"])}
    for weight in range(40, 52, 2):
        certificate = json.loads(
            (GENERATED / f"max_c_{weight}_scaled_rational.certificate.json").read_text()
        )
        shell_caps[weight] = int(certificate["physical_lattice_cap"])

    rows = []
    low_terms = []
    for weight, cap in shell_caps.items():
        pair_log = logadd(
            [
                float(log_coefficients[weight]),
                float(log_coefficients[256 - weight]),
            ]
        )
        term_log = math.log(cap) + pair_log
        low_terms.append(term_log)
        log2_cap = math.log2(cap)
        rows.append(
            {
                "weight": weight,
                "exact_Aw_C_cap": cap,
                "exact_Aw_C_cap_log2": log2_cap,
                "sufficient_application_cap_log2": application_caps[weight],
                "excess_over_application_cap_bits": (
                    log2_cap - application_caps[weight]
                ),
                "functional_term_log2": term_log / math.log(2.0),
            }
        )

    high_terms = [
        float(row["term_log_natural"])
        for row in diagnostic["proofish_low_prefix_reduction"][
            "generic_packing_shells"
        ]
    ]
    low_log = logadd(low_terms)
    total_log = logadd(low_terms + high_terms)
    result = {
        "classification": (
            "exact rational BCH-sandwich multiplicity caps evaluated with "
            "binary64 RandomStepConv coefficients"
        ),
        "shells": rows,
        "low_shell_margin_bits": -low_log / math.log(2.0),
        "combined_margin_bits": -total_log / math.log(2.0),
        "weight_38_term_margin_bits": -rows[0]["functional_term_log2"],
        "rigor_boundary": (
            "The multiplicity caps are exact. The transfer coefficients and "
            "reported margins use binary64 arithmetic."
        ),
    }
    output = GENERATED / "exact_bch_lp_caps.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
