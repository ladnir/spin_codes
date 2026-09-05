#!/usr/bin/env python3
"""Outward finite certificate for occupations 32 through 8576.

The certificate is conditional on the pointwise one-row comparison

    mu(v) <= 2^(1/1000) rho

for every nonzero 250-bit vector v.  The separate fanout-spectrum verifier
must establish this premise.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import sys

import numpy as np
from mpmath import iv


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_rm2sub_dense_occupation import (  # noqa: E402
    association_matrix,
    state_code_weights,
    verify_transpose_columns,
)
from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION  # noqa: E402
from certify_rm2sub_dense_compact_interval import (  # noqa: E402
    as_interval,
    collatz_upper,
    interval_three_state_transfer,
    upper_float,
)


iv.dps = 100
B = 250
K = 124
L = 8576
N = B * L
EPOCHS = N // 128
DISTANCE = 235_840
ETA = Fraction(1, 1000)
DEFAULT_DIAGNOSTIC = (
    WORKSTREAM
    / "bch250_124_parityfanout31x33_l256_finite_interval_cover_q32_8576_d11_diagnostic.json"
)


def exact_float_fraction(value: object) -> Fraction:
    return Fraction.from_float(float(value))


def log_integer(value: int):
    if value <= 0:
        raise ValueError("logarithm input must be positive")
    return iv.log(as_interval(Fraction(value)))


def finite_log_interval(
    occupation: int,
    *,
    p: Fraction,
    z: Fraction,
    radius: Fraction,
    prefactor: Fraction,
    eta: Fraction,
):
    choose = math.comb(L, occupation)
    log_choose = log_integer(choose)
    p_iv = as_interval(p)
    log_binomial_mass = (
        log_choose
        + occupation * iv.log(p_iv)
        + (L - occupation) * iv.log(1 - p_iv)
    )
    result = log_choose
    result += occupation * log_integer((1 << K) - 1)
    result -= occupation * iv.log(1 - as_interval(Fraction(1, 1 << B)))
    result -= B * log_binomial_mass
    result -= DISTANCE * iv.log(as_interval(z))
    result += EPOCHS * iv.log(as_interval(radius))
    result += iv.log(as_interval(prefactor))
    result += occupation * as_interval(eta) * iv.log(iv.mpf(2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostic", type=Path, default=DEFAULT_DIAGNOSTIC)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--eta-numerator", type=int, default=1)
    parser.add_argument("--eta-denominator", type=int, default=1000)
    parser.add_argument("--central-minimum-weight", type=int, default=1)
    parser.add_argument("--central-maximum-weight", type=int, default=B)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            WORKSTREAM
            / "bch250_124_parityfanout31x33_l256_dense_q32_8576_outward_d11.json"
        ),
    )
    args = parser.parse_args()
    if args.eta_numerator < 0 or args.eta_denominator <= 0:
        raise ValueError("invalid eta")
    if not 1 <= args.central_minimum_weight <= args.central_maximum_weight <= B:
        raise ValueError("invalid central weight interval")
    eta = Fraction(args.eta_numerator, args.eta_denominator)

    diagnostic = json.loads(args.diagnostic.read_text(encoding="utf-8"))
    if diagnostic["rejected_singletons"]:
        raise ValueError("diagnostic cover has rejected singletons")
    selected = json.loads(args.selection.read_text(encoding="utf-8"))["selected"]
    generator_words = [
        int(value, 16) for value in selected["A_generator_words_hex"]
    ]
    columns = [int(value, 16) for value in selected["B_columns_hex"]]
    verify_transpose_columns(
        generator_words=generator_words,
        columns=columns,
        output_bits=len(columns),
    )
    weights = state_code_weights(generator_words, len(columns))
    classes_array, counts_array, association_array = association_matrix(weights)
    classes = [int(value) for value in classes_array]
    counts = [int(value) for value in counts_array]
    association = [[int(value) for value in row] for row in association_array]
    state_space = 1 << len(generator_words)

    rows = []
    union_upper = iv.mpf(0)
    all_negative = True
    for source in diagnostic["intervals"]:
        lower = int(source["lower_occupation"])
        upper = int(source["upper_occupation"])
        p = exact_float_fraction(source["candidate_probability"])
        z = exact_float_fraction(source["z"])
        vector = [
            exact_float_fraction(value)
            for value in source["collatz_vector_binary64"]
        ]
        transfer = interval_three_state_transfer(
            candidate_probability=p,
            z=z,
            classes=classes,
            counts=counts,
            association=association,
            state_space=state_space,
            step_bits=len(columns),
        )
        radius = collatz_upper(transfer, vector)
        prefactor = vector[0] * max(Fraction(1, 1) / value for value in vector)
        endpoint_intervals = [
            finite_log_interval(
                occupation,
                p=p,
                z=z,
                radius=radius,
                prefactor=prefactor,
                eta=eta,
            )
            for occupation in (lower, upper)
        ]
        endpoint_uppers = [upper_float(value) for value in endpoint_intervals]
        maximum_upper = max(endpoint_uppers)
        accepted = maximum_upper < 0.0
        all_negative = all_negative and accepted
        width = upper - lower + 1
        union_upper += width * iv.exp(as_interval(Fraction.from_float(maximum_upper)))
        rows.append(
            {
                "lower_occupation": lower,
                "upper_occupation": upper,
                "candidate_probability_exact_binary64": str(p),
                "z_exact_binary64": str(z),
                "collatz_vector_exact_binary64": [str(value) for value in vector],
                "radius_upper": str(radius),
                "matrix_prefactor": str(prefactor),
                "endpoint_log_uppers_natural": endpoint_uppers,
                "accepted": accepted,
            }
        )

    union_probability_upper = upper_float(union_upper)
    comparison_to_2_neg_40 = union_upper < as_interval(
        Fraction(1, 1 << 40)
    )
    comparison_to_2_neg_41 = union_upper < as_interval(
        Fraction(1, 1 << 41)
    )
    comparison_to_2_neg_134 = union_upper < as_interval(
        Fraction(1, 1 << 134)
    )
    union_log2_upper = math.nextafter(
        math.log2(union_probability_upper), math.inf
    )
    result = {
        "schema": "bch250-124-rowlocal-fanout-dense-outward-v1",
        "status": (
            "CONDITIONAL_OUTWARD_CERTIFICATE"
            if all_negative and comparison_to_2_neg_40
            else "FAILED_OUTWARD_CERTIFICATE"
        ),
        "conditional_premise": {
            "claim": (
                "For every nonzero v, the expected wrapped-row multiplicity "
                "mu(v) in the stated central weight interval is at most "
                f"2^({args.eta_numerator}/{args.eta_denominator}) times the random-injection "
                "multiplicity rho=(2^124-1)/(2^250-1)."
            ),
            "eta_bits": f"{args.eta_numerator}/{args.eta_denominator}",
            "central_weight_interval": [
                args.central_minimum_weight,
                args.central_maximum_weight,
            ],
        },
        "claim": {
            "occupations": [32, 8576],
            "bad_weight_at_most": DISTANCE,
            "union_log2_upper": union_log2_upper,
            "union_margin_bits": -union_log2_upper,
            "comparison_to_2^-40": comparison_to_2_neg_40,
            "comparison_to_2^-41": comparison_to_2_neg_41,
            "comparison_to_2^-134": comparison_to_2_neg_134,
            "margin_display_is_not_used_for_comparisons": True,
        },
        "parameters": {
            "outer_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "output_bits": N,
            "inner_epochs": EPOCHS,
            "distance": DISTANCE,
            "relative_distance": DISTANCE / N,
            "accepted_intervals": len(rows),
        },
        "convexity_rule": diagnostic["convexity_rule"],
        "arithmetic": {
            "interval_dps": iv.dps,
            "method": (
                "All transfer, Collatz, logarithm, endpoint, exponential, "
                "and union inequalities use 100-digit interval arithmetic. "
                "Binary64 values select exact rational witnesses only."
            ),
        },
        "intervals": rows,
        "limitations": [
            "The receipt is conditional on the stated one-row density premise.",
            "Occupations 1 through 31 require separate outward receipts.",
            "Implementation equivalence is not part of this receipt.",
        ],
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: value for key, value in result.items() if key != "intervals"}, indent=2))
    print(f"output={args.output}")
    if result["status"] != "CONDITIONAL_OUTWARD_CERTIFICATE":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
