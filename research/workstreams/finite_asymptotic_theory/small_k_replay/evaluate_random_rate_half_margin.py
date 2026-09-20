#!/usr/bin/env python3
"""Compute the random-linear-code first-moment benchmark at rate one half."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "random_rate_half_d109_margin.json"
LOG2 = math.log(2.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-exponents", type=int, nargs="+", default=list(range(8, 21)))
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    rows = []
    for exponent in args.message_exponents:
        message_bits = 1 << exponent
        output_bits = 2 * message_bits
        cutoff = (
            args.distance_numerator * output_bits + args.distance_denominator - 1
        ) // args.distance_denominator
        weights = np.arange(1, cutoff, dtype=np.float64)
        log_ball = float(
            logsumexp(
                gammaln(output_bits + 1)
                - gammaln(weights + 1)
                - gammaln(output_bits - weights + 1)
            )
        )
        log_expected = (
            message_bits * LOG2
            + math.log1p(-math.ldexp(1.0, -message_bits))
            + log_ball
            - output_bits * LOG2
        )
        rows.append({
            "message_exponent": exponent,
            "message_bits": message_bits,
            "output_bits": output_bits,
            "distance_cutoff": cutoff,
            "relative_distance_lower": cutoff / output_bits,
            "expected_bad_log2_diagnostic": log_expected / LOG2,
            "margin_bits_diagnostic": -log_expected / LOG2,
            "closes_40_bits_diagnostic": log_expected < -40 * LOG2,
        })
    payload = {
        "schema": "random-linear-rate-half-first-moment-margin-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "experiment": (
            "sample a uniform binary linear map from k input bits to N=2k "
            "output bits and union-bound all nonzero messages"
        ),
        "claim": {
            "first_message_exponent_closing_40_bits": next(
                (row["message_exponent"] for row in rows if row["closes_40_bits_diagnostic"]),
                None,
            ),
            "rows": rows,
        },
        "limitations": [
            "Nearest binary64 log-gamma arithmetic is not outward rounded.",
            "This is a first-moment benchmark for a uniform random linear code, not a theorem about BCH, RM, routing, or RandomStepConv.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))


if __name__ == "__main__":
    main()
