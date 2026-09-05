#!/usr/bin/env python3
"""Conditional outward certificate for occupations 1 through 31."""

from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
import sys

from mpmath import iv


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

import certify_golay_ba_rm2sub_finite_q2_64 as base  # noqa: E402


B = 250
K = 124
L = 8576
N = B * L
DISTANCE = 235_840
MAXIMUM_Q = 31
ETA = Fraction(1, 1000)
Q1_DIAGNOSTIC = (
    WORKSTREAM
    / "bch250_124_parityfanout31x33_l256_packing_expected_q1_regular_d11_diagnostic.json"
)
Q2_DIAGNOSTIC = (
    WORKSTREAM
    / "bch250_124_parityfanout31x33_l256_packing_expected_q2_31_d11_diagnostic.json"
)
OUTPUT = (
    WORKSTREAM
    / "bch250_124_parityfanout31x33_l256_q1_31_outward_d11.json"
)


# Reuse the audited positive-arithmetic and RM2Sub recurrence implementation
# with this candidate's finite dimensions.
base.B = B
base.L = L
base.DISTANCE = DISTANCE
base.EPOCHS_PER_REGION = L // base.STEP_BITS


def load_witnesses() -> dict[int, tuple[float, float]]:
    result: dict[int, tuple[float, float]] = {}
    for path in (Q1_DIAGNOSTIC, Q2_DIAGNOSTIC):
        payload = json.loads(path.read_text(encoding="utf-8"))
        for row in payload["occupation_rows"]:
            occupation = int(row["active_regular_outer_blocks"])
            if not 1 <= occupation <= MAXIMUM_Q:
                continue
            p = float(row["best_candidate_probability"])
            z = math.exp(-math.exp(float(row["best_log_surprisal"])))
            if not 0.0 < p < 1.0 or not 0.0 < z < 1.0:
                raise ArithmeticError("invalid diagnostic witness")
            result[occupation] = (p, z)
    if set(result) != set(range(1, MAXIMUM_Q + 1)):
        raise ArithmeticError("diagnostics omit a required occupation")
    return result


def eta_factor_upper() -> float:
    iv.dps = 100
    value = iv.power(iv.mpf(2), iv.mpf(1) / 1000)
    return math.nextafter(float(value.b), math.inf)


def bernoulli_density_envelope_upper(p: float) -> tuple[float, int]:
    q_lower = base.down(1.0 - p)
    rho = base.ratio_upper((1 << K) - 1, (1 << B) - 1)
    numerator = base.mul_up(eta_factor_upper(), rho)
    best = 0.0
    best_weight = -1
    for weight in range(1, B + 1):
        denominator = base.mul_down(
            base.power_down(p, weight),
            base.power_down(q_lower, B - weight),
        )
        candidate = base.div_up(numerator, denominator)
        if candidate > best:
            best = candidate
            best_weight = weight
    return best, best_weight


def scaled_log2_upper(value: base.Scaled) -> float:
    if value[0] <= 0.0:
        return -math.inf
    return math.nextafter(math.log2(value[0]) + value[1], math.inf)


def main() -> None:
    if L % base.STEP_BITS:
        raise ArithmeticError("step bits must divide the row count")
    witnesses = load_witnesses()
    activation, live = base.load_activation_and_live()

    envelopes: dict[str, tuple[float, int]] = {}
    caches: dict[str, list[base.Scaled]] = {}
    for p, z in sorted(set(witnesses.values())):
        p_key = p.hex()
        envelopes.setdefault(p_key, bernoulli_density_envelope_upper(p))
        impulses = base.impulse_matrices_upper(z, activation, live)
        candidates = base.candidate_epoch_matrices(impulses, p)
        regions = base.region_matrices_upper(candidates, MAXIMUM_Q)
        caches[f"{p.hex()}|{z.hex()}"] = base.total_moments_upper(regions)

    rows = []
    aggregate = base.SCALED_ZERO
    for occupation in range(1, MAXIMUM_Q + 1):
        p, z = witnesses[occupation]
        envelope, maximizing_weight = envelopes[p.hex()]
        moment = caches[f"{p.hex()}|{z.hex()}"][occupation]
        chernoff = base.scaled_multiply_up(
            moment,
            base.scaled_power_up(
                base.scaled_from_float(base.div_up(1.0, z)), DISTANCE
            ),
        )
        inner = (
            chernoff
            if base.scaled_at_most_power_of_two(chernoff, 0)
            else base.SCALED_ONE
        )
        contribution = base.scaled_multiply_up(
            base.scaled_from_float(
                base.integer_upper(math.comb(L, occupation))
            ),
            base.scaled_multiply_up(
                base.scaled_power_up(
                    base.scaled_from_float(envelope), occupation
                ),
                inner,
            ),
        )
        aggregate = base.scaled_add_up(aggregate, contribution)
        rows.append(
            {
                "occupation": occupation,
                "candidate_probability_exact_binary64": p.hex(),
                "z_exact_binary64": z.hex(),
                "density_envelope_upper_hex": envelope.hex(),
                "envelope_maximizing_weight": maximizing_weight,
                "inner_probability_upper_scaled": base.scaled_hex(inner),
                "contribution_upper_scaled": base.scaled_hex(contribution),
                "contribution_log2_upper_display": scaled_log2_upper(
                    contribution
                ),
            }
        )

    aggregate_log2 = scaled_log2_upper(aggregate)
    comparison_to_2_neg_40 = base.scaled_at_most_power_of_two(
        aggregate, 40
    )
    comparison_to_2_neg_48 = base.scaled_at_most_power_of_two(
        aggregate, 48
    )
    result = {
        "schema": "bch250-124-rowlocal-fanout-q1-31-outward-v1",
        "status": (
            "CONDITIONAL_OUTWARD_CERTIFICATE"
            if comparison_to_2_neg_40
            else "FAILED_OUTWARD_CERTIFICATE"
        ),
        "conditional_premise": {
            "claim": (
                "For every nonzero v, mu(v) <= 2^(1/1000) rho, where "
                "rho=(2^124-1)/(2^250-1)."
            ),
            "eta_bits": "1/1000",
        },
        "claim": {
            "occupations": [1, MAXIMUM_Q],
            "bad_weight_at_most": DISTANCE,
            "expected_bad_count_upper_scaled": base.scaled_hex(aggregate),
            "union_log2_upper": aggregate_log2,
            "union_margin_bits": -aggregate_log2,
            "comparison_to_2^-40": comparison_to_2_neg_40,
            "comparison_to_2^-48": comparison_to_2_neg_48,
            "margin_display_is_not_used_for_comparisons": True,
        },
        "parameters": {
            "outer_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "epochs_per_region": base.EPOCHS_PER_REGION,
        },
        "arithmetic": {
            "format": "IEEE-754 binary64 with one-ULP outward operations",
            "eta_factor": (
                "The upper endpoint of a 100-digit interval evaluation of "
                "2^(1/1000), advanced one binary64 ULP."
            ),
            "transcendental_scope": (
                "The exponential evaluations select exact binary64 z "
                "witnesses and assert no inequality."
            ),
        },
        "occupation_rows": rows,
        "limitations": [
            "The receipt is conditional on the stated one-row density premise.",
            "Occupations 32 through 8576 use a separate outward receipt.",
            "Implementation equivalence is not part of this receipt.",
        ],
    }
    OUTPUT.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({key: value for key, value in result.items() if key != "occupation_rows"}, indent=2))
    print(f"output={OUTPUT}")
    if result["status"] != "CONDITIONAL_OUTWARD_CERTIFICATE":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
