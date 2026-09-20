#!/usr/bin/env python3
"""Evaluate the full-character certificate needed for terminal-zero bounds."""

from __future__ import annotations

import collections
import hashlib
import json
import math
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import accumulate, build_apply
from probe_riffle_dp_g4_component_mixing import contribution_states


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
OUTER = EXPLORATIONS / "riffle_dp_g4_g1_support33_refutation.json"
BIAS_PROBE = EXPLORATIONS / "riffle_dp_g4_g2_full_character_bias_probe.json"
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_terminal_zero_parseval_gate.json"
PACKET_POSITIONS = 524_352
PACKET_SUPPORT = 33
STATE_BITS = 64


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log2_interval(value: Fraction, places: int = 12) -> list[str]:
    with localcontext() as context:
        context.prec = 100
        estimate = (
            Decimal(value.numerator).ln() - Decimal(value.denominator).ln()
        ) / Decimal(2).ln()
        unit = Decimal(1).scaleb(-places)
        return [
            format(estimate.quantize(unit, rounding=ROUND_FLOOR), "f"),
            format(estimate.quantize(unit, rounding=ROUND_CEILING), "f"),
        ]


def packet_profile(receipt: dict) -> tuple[int, ...]:
    counts: collections.Counter[int] = collections.Counter()
    for local_word in receipt["local_words"]:
        word = int(local_word["codeword_hex"], 16)
        counts.update(
            (word >> shift) & 15
            for shift in range(0, 128, 4)
            if (word >> shift) & 15
        )
    profile = tuple(counts[value] for value in range(1, 16))
    if sum(profile) != PACKET_SUPPORT:
        raise RuntimeError("terminal-zero gate: support-33 profile changed")
    return profile


def main() -> None:
    outer = json.loads(OUTER.read_text())
    bias_probe = json.loads(BIAS_PROBE.read_text())
    bias_numerators = {
        row["packet_value"]: row["maximum_found_absolute_bias_numerator"]
        for row in bias_probe["value_rows"]
    }
    if set(bias_numerators) != set(range(1, 16)):
        raise RuntimeError("terminal-zero gate: incomplete bias probe")

    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    collision_rows = []
    for value in range(1, 16):
        states = contribution_states(step, value)
        _, multiplicities = np.unique(states, return_counts=True)
        collision_sum = int(
            np.dot(multiplicities.astype(object), multiplicities.astype(object))
        )
        maximum_multiplicity = int(multiplicities.max())
        if len(multiplicities) != PACKET_POSITIONS or maximum_multiplicity != 1:
            raise RuntimeError("terminal-zero gate: contribution map is not injective")
        collision_rows.append(
            {
                "packet_value": value,
                "distinct_contribution_states": len(multiplicities),
                "collision_sum": collision_sum,
                "maximum_multiplicity": maximum_multiplicity,
            }
        )

    distinct_probability = Fraction(
        math.prod(range(PACKET_POSITIONS - PACKET_SUPPORT + 1, PACKET_POSITIONS + 1)),
        PACKET_POSITIONS**PACKET_SUPPORT,
    )
    nontrivial_parseval_mass = Fraction(1 << STATE_BITS, PACKET_POSITIONS) - 1

    profile_families: dict[tuple[int, ...], list[str]] = collections.defaultdict(list)
    for receipt in outer["receipts"]:
        profile_families[packet_profile(receipt)].append(receipt["family_id"])

    rows = []
    aggregate = Fraction()
    for profile, family_ids in sorted(profile_families.items()):
        best = None
        for anchor_index, count in enumerate(profile):
            if count < 2:
                continue
            product = Fraction(1)
            for value_index, value_count in enumerate(profile):
                exponent = value_count - (2 if value_index == anchor_index else 0)
                product *= Fraction(
                    bias_numerators[value_index + 1], PACKET_POSITIONS
                ) ** exponent
            iid_bound = Fraction(1, 1 << STATE_BITS) * (
                1 + nontrivial_parseval_mass * product
            )
            conditioned_bound = min(Fraction(1), iid_bound / distinct_probability)
            candidate = (conditioned_bound, anchor_index + 1, product, iid_bound)
            if best is None or candidate[0] < best[0]:
                best = candidate
        if best is None:
            raise RuntimeError("terminal-zero gate: profile has no repeated value")
        conditioned_bound, anchor_value, product, iid_bound = best
        aggregate += len(family_ids) * conditioned_bound
        rows.append(
            {
                "packet_value_counts_1_through_15": list(profile),
                "family_ids": family_ids,
                "parseval_anchor_value": anchor_value,
                "nonanchor_bias_product_numerator": str(product.numerator),
                "nonanchor_bias_product_denominator": str(product.denominator),
                "iid_terminal_zero_bound_numerator": str(iid_bound.numerator),
                "iid_terminal_zero_bound_denominator": str(iid_bound.denominator),
                "conditioned_terminal_zero_bound_numerator": str(
                    conditioned_bound.numerator
                ),
                "conditioned_terminal_zero_bound_denominator": str(
                    conditioned_bound.denominator
                ),
                "conditioned_terminal_zero_log2_interval": log2_interval(
                    conditioned_bound
                ),
            }
        )

    def aggregate_for_caps(caps: dict[int, Fraction]) -> Fraction:
        result = Fraction()
        for profile, family_ids in profile_families.items():
            best_bound = None
            for anchor_index, count in enumerate(profile):
                if count < 2:
                    continue
                product = Fraction(1)
                for value_index, value_count in enumerate(profile):
                    exponent = value_count - (2 if value_index == anchor_index else 0)
                    product *= caps[value_index + 1] ** exponent
                iid_bound = Fraction(1, 1 << STATE_BITS) * (
                    1 + nontrivial_parseval_mass * product
                )
                conditioned_bound = min(
                    Fraction(1), iid_bound / distinct_probability
                )
                if best_bound is None or conditioned_bound < best_bound:
                    best_bound = conditioned_bound
            if best_bound is None:
                raise RuntimeError("terminal-zero gate: profile has no coarse anchor")
            result += len(family_ids) * best_bound
        return result

    coarse_caps = {
        value: Fraction(5, 8) if value == 15 else Fraction(1, 2)
        for value in range(1, 16)
    }
    coarse_aggregate = aggregate_for_caps(coarse_caps)

    payload = {
        "schema": "riffle-dp-2lap-g4-terminal-zero-parseval-gate-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "DIAGNOSTIC_CONDITIONAL_GATE",
        "outer_receipt_sha256": digest(OUTER),
        "full_character_bias_probe_sha256": digest(BIAS_PROBE),
        "packet_positions": PACKET_POSITIONS,
        "packet_support": PACKET_SUPPORT,
        "exact_distinctness_probability_numerator": str(
            distinct_probability.numerator
        ),
        "exact_distinctness_probability_denominator": str(
            distinct_probability.denominator
        ),
        "exact_contribution_injectivity": True,
        "collision_rows": collision_rows,
        "conditional_theorem": (
            "Suppose B_v is a proved upper bound on every nontrivial full-state "
            "character sum for packet value v. Analyze independent uniform packet "
            "cells by Fourier inversion and Parseval. Then divide the resulting "
            "point-probability bound by the exact probability that all 33 cells "
            "are distinct. The displayed rows instantiate B_v with the maxima "
            "found by the diagnostic hill search."
        ),
        "parseval_identity": (
            "Injectivity gives sum_chi |W_v(chi)/M|^2 = 2^64/M. "
            "Two copies of the anchor value pay this exact sum; every other "
            "factor pays its assumed full-character maximum."
        ),
        "profile_rows": rows,
        "conditional_aggregate_over_26_words": {
            "numerator": str(aggregate.numerator),
            "denominator": str(aggregate.denominator),
            "log2_interval": log2_interval(aggregate),
            "below_2^-40": aggregate < Fraction(1, 1 << 40),
        },
        "coarse_certificate_target": {
            "assumed_full_character_caps": {
                str(value): (
                    "5/8" if value == 15 else "1/2"
                )
                for value in range(1, 16)
            },
            "aggregate_numerator": str(coarse_aggregate.numerator),
            "aggregate_denominator": str(coarse_aggregate.denominator),
            "aggregate_log2_interval": log2_interval(coarse_aggregate),
            "below_2^-40": coarse_aggregate < Fraction(1, 1 << 40),
            "interpretation": (
                "It is unnecessary to certify every diagnostic maximum. The "
                "coarse caps 5/8 for value 15 and 1/2 for every other nonzero "
                "nibble already close the support-33 terminal-zero row."
            ),
        },
        "scope_limitation": (
            "The contribution-state injectivity, Parseval calculation, and "
            "conditioning reduction are exact. The full-character maxima are "
            "diagnostic lower bounds on the required caps, not certificates. "
            "Therefore the displayed terminal-zero values are proof targets, "
            "not rigorous upper bounds."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    worst = max(rows, key=lambda row: Fraction(
        int(row["conditioned_terminal_zero_bound_numerator"]),
        int(row["conditioned_terminal_zero_bound_denominator"]),
    ))
    print(f"unique_support33_profiles={len(rows)}")
    print("contribution_maps_injective=True")
    print(
        "worst_conditional_profile_log2="
        f"{worst['conditioned_terminal_zero_log2_interval']}"
    )
    print(
        "conditional_aggregate_log2="
        f"{payload['conditional_aggregate_over_26_words']['log2_interval']}"
    )
    print(
        "conditional_aggregate_below_2^-40="
        f"{payload['conditional_aggregate_over_26_words']['below_2^-40']}"
    )
    print(f"output={OUTPUT}")
    print("status=DIAGNOSTIC_FULL_CHARACTER_CERTIFICATE_GATE")


if __name__ == "__main__":
    main()
