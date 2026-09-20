#!/usr/bin/env python3
"""Derive the exact support-33 PacketMul terminal-zero Fourier gate."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
MANIFEST = CANDIDATE / "manifest.json"
PRIMARY = CANDIDATE / "receipts" / "goal01_zero_symbol_primary.json"
INDEPENDENT = CANDIDATE / "receipts" / "goal01_zero_symbol_independent.json"
SEARCH = CANDIDATE / "receipts" / "goal01_mixed_character_search.json"
OUTER = ROOT / "explorations" / "riffle_dp_g4_g1_support33_refutation.json"
OUTPUT = CANDIDATE / "receipts" / "goal01_terminal_zero_gate.json"
M = 524_352
H = 33
WORD_COUNT = 26
STATE_BITS = 64


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fraction_row(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


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


def root_interval(value: Fraction, degree: int, decimal_places: int = 15) -> tuple[Fraction, Fraction]:
    """Return adjacent decimal rationals that rigorously bracket value^(1/degree)."""
    scale = 10**decimal_places
    low = 0
    high = scale + 1
    while low + 1 < high:
        middle = (low + high) // 2
        candidate = Fraction(middle, scale)
        if candidate**degree <= value:
            low = middle
        else:
            high = middle
    return Fraction(low, scale), Fraction(high, scale)


def aggregate_bound(parseval: Fraction, distinct_probability: Fraction, cap: Fraction) -> Fraction:
    iid = Fraction(1, 1 << STATE_BITS) * (
        1 + (parseval - 1) * cap ** (H - 2)
    )
    return WORD_COUNT * iid / distinct_probability


def packet_support(receipt: dict) -> int:
    support = 0
    for local_word in receipt["local_words"]:
        word = int(local_word["codeword_hex"], 16)
        support += sum(
            ((word >> shift) & 15) != 0 for shift in range(0, 128, 4)
        )
    return support


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    independent = json.loads(INDEPENDENT.read_text(encoding="utf-8"))
    search = json.loads(SEARCH.read_text(encoding="utf-8"))
    outer = json.loads(OUTER.read_text(encoding="utf-8"))
    if manifest["candidate_id"] != "riffle_packetmul_2lap_g4":
        raise RuntimeError("candidate manifest changed")
    if not independent["all_primary_checks_passed"]:
        raise RuntimeError("independent algebraic replay did not pass")
    if len(outer["receipts"]) != WORD_COUNT:
        raise RuntimeError("authenticated outer word count changed")
    support_histogram = Counter(packet_support(receipt) for receipt in outer["receipts"])
    if support_histogram != Counter({H: WORD_COUNT}):
        raise RuntimeError("authenticated outer support row changed")

    collision = primary["cross_value_collisions"]
    parseval_row = collision["normalized_parseval_sum"]
    parseval = Fraction(
        int(parseval_row["numerator"]), int(parseval_row["denominator"])
    )
    distinct_probability = Fraction(
        math.prod(range(M - H + 1, M + 1)),
        M**H,
    )

    # Solve the exact equality at aggregate probability 2^-40.
    required_power = (
        distinct_probability * (1 << (STATE_BITS - 40)) / WORD_COUNT - 1
    ) / (parseval - 1)
    if not 0 < required_power < 1:
        raise RuntimeError("required cap is outside (0,1)")
    cap_lower, cap_upper = root_interval(required_power, H - 2)
    if not cap_lower ** (H - 2) <= required_power < cap_upper ** (H - 2):
        raise RuntimeError("required cap interval does not bracket the exact root")
    q_lower = (15 * cap_lower + 1) / 16
    q_upper = (15 * cap_upper + 1) / 16

    simple_q_cap = Fraction(5, 8)
    simple_beta_cap = Fraction(3, 5)
    if (16 * simple_q_cap - 1) / 15 != simple_beta_cap:
        raise RuntimeError("simple q/beta cap conversion failed")
    simple_aggregate = aggregate_bound(parseval, distinct_probability, simple_beta_cap)
    if not simple_aggregate < Fraction(1, 1 << 40):
        raise RuntimeError("simple mixed-character lemma does not close the row")
    half_q_aggregate = aggregate_bound(parseval, distinct_probability, Fraction(7, 15))

    best = search["best_exact_character_found"]
    observed_cap = Fraction(
        int(best["absolute_beta_numerator"]), int(best["beta_denominator"])
    )
    observed_aggregate = aggregate_bound(parseval, distinct_probability, observed_cap)

    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal01-terminal-zero-gate-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "candidate_id": manifest["candidate_id"],
        "evidence_label": "EXACT_REDUCTION_CONDITIONAL_MIXED_LEMMA",
        "active_manifest_sha256": digest(MANIFEST),
        "primary_receipt_sha256": digest(PRIMARY),
        "independent_receipt_sha256": digest(INDEPENDENT),
        "mixed_search_receipt_sha256": digest(SEARCH),
        "authenticated_outer_receipt_sha256": digest(OUTER),
        "authenticated_support_histogram": {
            str(support): count for support, count in sorted(support_histogram.items())
        },
        "packet_cells": M,
        "packet_support": H,
        "authenticated_word_count": WORD_COUNT,
        "exact_distinctness_probability": fraction_row(distinct_probability),
        "exact_packetmul_parseval_sum": fraction_row(parseval),
        "reduction": (
            "For each supported source packet, the independent nonzero multiplier and "
            "iid cell draw make z_y(p) iid uniform on the 15*M labeled samples. Fourier "
            "inversion gives 2^-64 sum_chi beta_chi^33. Taking absolute values, paying "
            "two powers with exact Parseval, and paying the other 31 with a uniform cap "
            "B gives 2^-64[1+(Parseval-1)B^31]. Division by Pr[D] conditions the iid "
            "cells on being distinct, recovering the ordered packet-permutation law."
        ),
        "aggregate_bound_formula": (
            "26 * 2^-64 * [1 + (Parseval-1)*B^31] / Pr[D]"
        ),
        "required_uniform_absolute_beta_cap": {
            "exact_power_identity": (
                "B_req^31 = (Pr[D]*2^24/26 - 1)/(Parseval-1)"
            ),
            "lower_decimal": format(float(cap_lower), ".15f"),
            "upper_decimal": format(float(cap_upper), ".15f"),
            "rigorous_lower_rational": fraction_row(cap_lower),
            "rigorous_upper_rational": fraction_row(cap_upper),
            "equivalent_q_upper_lower_decimal": format(float(q_lower), ".15f"),
            "equivalent_q_upper_upper_decimal": format(float(q_upper), ".15f"),
        },
        "simple_sufficient_lemma_instantiation": {
            "q_cap": "5/8",
            "implied_absolute_beta_cap": "3/5",
            "aggregate_bound": fraction_row(simple_aggregate),
            "aggregate_log2_interval": log2_interval(simple_aggregate),
            "below_2^-40": simple_aggregate < Fraction(1, 1 << 40),
            "security_margin_bits_interval": [
                format(Decimal(-40) - Decimal(log2_interval(simple_aggregate)[1]), "f"),
                format(Decimal(-40) - Decimal(log2_interval(simple_aggregate)[0]), "f"),
            ],
        },
        "stronger_half_q_reference": {
            "q_cap": "1/2",
            "implied_absolute_beta_cap": "7/15",
            "aggregate_log2_interval": log2_interval(half_q_aggregate),
        },
        "diagnostic_search_comparison": {
            "best_character_hex": best["character_hex"],
            "best_zero_symbol_count": best["zero_symbol_count"],
            "observed_absolute_beta": fraction_row(observed_cap),
            "aggregate_if_observed_value_were_a_proved_cap_log2_interval": log2_interval(
                observed_aggregate
            ),
            "interpretation": (
                "This is refutation evidence only. It is not a global upper bound."
            ),
        },
        "remaining_mixed_character_lemma": (
            "For every nonzero terminal character chi, at most 5M/8 packet cells p "
            "have a_chi(p)=0. Equivalently, q_chi <= 5/8. Since beta_chi is at "
            "least -1/15 and at most 3/5, this implies |beta_chi| <= 3/5 and "
            "closes the authenticated support-33 terminal-zero row with the displayed margin."
        ),
        "scope_limitation": (
            "The Fourier/Parseval/conditioning reduction and numerical threshold are "
            "exact. The remaining mixed-character lemma is explicit but unproved. This "
            "receipt does not prove the full construction or other placement strata."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        "B_req_interval="
        f"[{payload['required_uniform_absolute_beta_cap']['lower_decimal']},"
        f"{payload['required_uniform_absolute_beta_cap']['upper_decimal']}]"
    )
    print(
        "simple_q_cap_5_over_8_aggregate_log2="
        f"{payload['simple_sufficient_lemma_instantiation']['aggregate_log2_interval']}"
    )
    print(f"output={OUTPUT}")
    print("status=EXACT_TERMINAL_ZERO_REDUCTION_OPEN_MIXED_LEMMA")


if __name__ == "__main__":
    main()
