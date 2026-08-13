#!/usr/bin/env python3
"""Directed-decimal certificate for the PA1 scalar group envelope.

The exact rational Collatz certificate at output pole 997/1000 supplies rho
and the reached-vector prefactor.  For a word occupying g distinct groups,
the uniform gap composition and every zero-state restart give the positive
episode sum used by ``probe_striped_group_low_ledger.py``.

This verifier evaluates the zero- and one-restart terms with Decimal rounding
toward +infinity.  All terms with at least two restarts are bounded by the
corresponding binomial tail after replacing survival by one.  The resulting
finite Decimal is converted exactly to Fraction and compared with
``(181/1000)^g`` for every relevant support g=21..106.
"""

from __future__ import annotations

import json
import math
from decimal import Decimal, ROUND_CEILING, localcontext
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CERTIFICATES = ROOT / "systematic_pa1_geometric_certificates.json"
B = 32768
DISTANCE = 9 * (1 << 21) // 100
FIRST_SUPPORT = 21
LAST_SUPPORT = 106
TURN_OFF_BITS = 39
OUTPUT_POLE = Fraction(997, 1000)
TARGET_GROUP_POLE = Fraction(181, 1000)
PRECISION = 90


def decimal_upper(value: Fraction) -> Decimal:
    with localcontext() as context:
        context.prec = PRECISION
        context.rounding = ROUND_CEILING
        return Decimal(value.numerator) / Decimal(value.denominator)


def load_certificate() -> tuple[Fraction, Fraction]:
    payload = json.loads(CERTIFICATES.read_text(encoding="utf-8"))
    if payload.get("status") != "EXACT_RATIONAL_COLLATZ_CERTIFICATE":
        raise SystemExit("PA1 envelope: invalid Collatz certificate status")
    for row in payload["certificates"]:
        if (
            row["pole_numerator"] == OUTPUT_POLE.numerator
            and row["pole_denominator"] == OUTPUT_POLE.denominator
        ):
            rho = Fraction(
                int(row["rho_numerator"]), int(row["rho_denominator"])
            )
            prefactor = Fraction(
                int(row["prefactor_numerator"]),
                int(row["prefactor_denominator"]),
            )
            if not 0 < rho < 1 or not 0 < prefactor < 1:
                raise SystemExit("PA1 envelope: invalid Collatz constants")
            return rho, prefactor
    raise SystemExit("PA1 envelope: missing pole 997/1000")


def survival_table(rho: Fraction, prefactor: Fraction) -> list[Decimal]:
    with localcontext() as context:
        context.prec = PRECISION
        context.rounding = ROUND_CEILING
        rho_decimal = Decimal(rho.numerator) / Decimal(rho.denominator)
        pole_decimal = Decimal(OUTPUT_POLE.numerator) / Decimal(
            OUTPUT_POLE.denominator
        )
        prefactor_decimal = Decimal(prefactor.numerator) / Decimal(
            prefactor.denominator
        )
        raw = prefactor_decimal * (Decimal(1) / pole_decimal) ** DISTANCE
        result = [Decimal(1)] * (B + 1)
        for live in range(1, B + 1):
            if live > 1:
                raw *= rho_decimal
            result[live] = min(Decimal(1), +raw)
        # Retain the explicit short-live branch of the theorem-facing gap
        # envelope even if the rounded Collatz expression is already >=1.
        for live in range(DISTANCE // 64 + 1):
            result[live] = Decimal(1)
        return result


def weighted_gap_sum(
    group_support: int,
    cancellations: int,
    survival: list[Decimal],
) -> Decimal:
    slack = B - group_support
    remaining_gaps = group_support - cancellations - 1
    first = 1  # C(off+e,e) at off=0
    second = math.comb(slack + remaining_gaps, remaining_gaps)
    total = Decimal(0)
    with localcontext() as context:
        context.prec = PRECISION
        context.rounding = ROUND_CEILING
        for off in range(slack + 1):
            total += Decimal(first * second) * survival[B - off]
            if off == slack:
                break
            first = first * (off + cancellations + 1) // (off + 1)
            denominator = slack - off + remaining_gaps
            second = second * (slack - off) // denominator
        return +total


def group_bound(
    group_support: int,
    survival: list[Decimal],
    restart: Fraction,
) -> Decimal:
    denominator = math.comb(B, group_support)
    target = TARGET_GROUP_POLE**group_support
    first_omitted = group_support
    omitted_tail = Fraction(0)
    for cancellations in range(2, group_support):
        first_term = math.comb(group_support - 1, cancellations) * restart**cancellations
        next_ratio = (
            Fraction(group_support - 1 - cancellations, cancellations + 1)
            * restart
        )
        if next_ratio >= 1:
            continue
        tail_bound = first_term / (1 - next_ratio)
        if tail_bound <= target / (1 << 30):
            first_omitted = cancellations
            omitted_tail = tail_bound
            break
    with localcontext() as context:
        context.prec = PRECISION
        context.rounding = ROUND_CEILING
        restart_decimal = Decimal(restart.numerator) / Decimal(restart.denominator)
        total = Decimal(0)
        for cancellations in range(first_omitted):
            gaps = weighted_gap_sum(group_support, cancellations, survival)
            total += (
                Decimal(math.comb(group_support - 1, cancellations))
                * restart_decimal**cancellations
                * gaps
                / Decimal(denominator)
            )
        # Survival is at most one.  Successive omitted binomial restart terms
        # have decreasing ratios, so the first omitted term divided by
        # one minus its next ratio bounds the complete remaining tail.
        if omitted_tail:
            total += Decimal(omitted_tail.numerator) / Decimal(
                omitted_tail.denominator
            )
        return +total


def main() -> None:
    rho, prefactor = load_certificate()
    restart = Fraction(1, 1 << TURN_OFF_BITS) * prefactor / rho
    survival = survival_table(rho, prefactor)
    worst_ratio = Fraction(0)
    worst_support = 0
    print("directed-decimal PA1 group-envelope certificate")
    print(f"precision={PRECISION} output_pole={OUTPUT_POLE}")
    for support in range(FIRST_SUPPORT, LAST_SUPPORT + 1):
        bound_decimal = group_bound(support, survival, restart)
        bound = Fraction(bound_decimal)
        target = TARGET_GROUP_POLE**support
        if bound > target:
            raise SystemExit(
                f"PA1 envelope: support {support} exceeds (181/1000)^g"
            )
        ratio = bound / target
        if ratio > worst_ratio:
            worst_ratio = ratio
            worst_support = support
    print(f"worst_support={worst_support}")
    print(
        "worst_ratio_log2_upper="
        f"{math.log2(worst_ratio.numerator) - math.log2(worst_ratio.denominator):.12f}"
    )
    print("supports_21_106_le_(181/1000)^g=PASS")
    print("status=EXACT_COLLATZ_DIRECTED_DECIMAL_GAP_SUM")


if __name__ == "__main__":
    main()
