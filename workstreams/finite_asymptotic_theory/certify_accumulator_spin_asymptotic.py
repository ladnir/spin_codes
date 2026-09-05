#!/usr/bin/env python3
"""Certify a Random-block Accumulator SPIN parameter point.

The proof formulas are in ACCUMULATOR_SPIN_ASYMPTOTIC.md.  This script uses
mpmath interval arithmetic for the strict exponent and schedule inequalities.
"""

from __future__ import annotations

import argparse
import json
import math

import mpmath as mp


def interval_argument(value: str) -> mp.ctx_iv.ivmpf:
    if "/" not in value:
        return mp.iv.mpf(value)
    numerator, denominator = value.split("/", maxsplit=1)
    return mp.iv.mpf(numerator) / mp.iv.mpf(denominator)


def endpoints(value: mp.ctx_iv.ivmpf) -> list[float]:
    return [float(value.a), float(value.b)]


def log2(value: mp.ctx_iv.ivmpf) -> mp.ctx_iv.ivmpf:
    return mp.iv.log(value) / mp.iv.log(2)


def rate_half_ceiling(block_constant: float) -> float:
    if block_constant <= 2.0:
        return 0.0
    rho = 2.0 ** (0.5 - 1.0 / block_constant) - 1.0
    return (1.0 - math.sqrt(1.0 - rho * rho)) / 2.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rate", default="1/2")
    parser.add_argument("--delta", default="37/10000")
    parser.add_argument("--block-constant", default="3")
    args = parser.parse_args()

    rate = interval_argument(args.rate)
    delta = interval_argument(args.delta)
    block_constant = interval_argument(args.block_constant)
    one = mp.iv.mpf(1)
    two = mp.iv.mpf(2)

    rho = two * mp.iv.sqrt(delta * (one - delta))
    q = two ** (rate - one) * (one + rho)
    exponent = -log2(q)
    schedule_margin = block_constant * exponent - one

    if float(delta.a) <= 0.0 or float(delta.b) >= 0.5:
        raise AssertionError("delta must lie strictly between zero and one half")
    if float(exponent.a) <= 0.0:
        raise AssertionError("the sparse outer exponent is not positive")
    if float(schedule_margin.a) <= 0.0:
        raise AssertionError("the selected point lacks a strict schedule margin")

    table = []
    for constant in (3, 4, 5, 6, 8, 10, 16, 17, 32):
        table.append(
            {
                "block_constant": constant,
                "rate_half_distance_ceiling_binary64": rate_half_ceiling(constant),
            }
        )

    payload = {
        "schema": "accumulator-spin-asymptotic-certificate-v1",
        "status": "CERTIFIED_STRICT_PARAMETER_POINT",
        "arithmetic": {
            "library": "mpmath.iv",
            "mpmath_version": mp.__version__,
        },
        "parameters": {
            "rate": args.rate,
            "delta": args.delta,
            "block_constant": args.block_constant,
        },
        "certified_intervals": {
            "accumulator_contraction_rho": endpoints(rho),
            "outer_base_q": endpoints(q),
            "outer_decay_log2": endpoints(exponent),
            "schedule_margin": endpoints(schedule_margin),
        },
        "rate_half_ceiling_table": table,
        "table_arithmetic": "binary64 rendering of the proved closed formula",
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
