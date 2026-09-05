#!/usr/bin/env python3
"""All-occupation envelope diagnostic for one random [256,128] constituent.

The setup samples one binary generator and retains it when it has full row
rank and its complete weight spectrum satisfies a pointwise factor envelope.
The retained constituent is reused in every outer row.  Conditional on this
event, the fixed-code message sum factorizes across active rows, so one
occupation recurrence covers Q=1,...,L.

Nearest binary64 arithmetic makes the distance result diagnostic.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import evaluate_repeated_random512_randomstepconv_g1 as base


WORKSTREAM = Path(__file__).resolve().parent
B = 256
K = 128
L = 8192
N = 1 << 21
D = (109 * N + 999) // 1000
MEMORY = 22
DEFAULT_OUTPUT = (
    WORKSTREAM / "single_random_constituent_B256_allq_factor8_s22_d109.json"
)


def configure(factor: int) -> None:
    if factor <= 1:
        raise ValueError("spectrum factor must exceed one")
    base.B = B
    base.K = K
    base.L = L
    base.ACTIVE_ROWS = L
    base.N = N
    base.D = D
    base.MEMORY = MEMORY
    base.SPECTRUM_FACTOR = factor


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum-factor", type=int, default=8)
    parser.add_argument("--memory-bits", type=int, default=MEMORY)
    parser.add_argument("--grid-min", type=float, default=-9.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--refine-min", type=float, default=0.65)
    parser.add_argument("--refine-max", type=float, default=0.80)
    parser.add_argument("--refine-step", type=float, default=0.025)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    configure(args.spectrum_factor)
    payload = base.evaluate(args)
    payload["schema"] = "single-random-constituent-randomstepconv-allq-v1"
    payload["status"] = "BINARY64_DIAGNOSTIC_CONDITIONAL_SPECTRUM_EVENT"
    payload["parameters"]["outer_code"] = "one random [256,128] linear constituent"
    payload["parameters"]["distance_numerator"] = 109
    payload["parameters"]["distance_denominator"] = 1000
    payload["constituent_event"]["event"] = (
        "full rank and A_w <= "
        f"{args.spectrum_factor} E[A_w] for every 1<=w<=256"
    )
    payload["probability_space"] = (
        "sample candidate generators until retaining one that satisfies the "
        "stated spectrum event; reuse the retained constituent in every row; "
        "then sample independent row-coordinate permutations, region "
        "permutations, and RandomStepConv maps"
    )
    payload["limitations"] = [
        "The distance sum uses nearest binary64 log arithmetic.",
        "The constituent-event probability uses elementary union bounds.",
        "The result is not yet an outward-rounded certificate.",
        "Efficiently testing the complete constituent spectrum is not supplied.",
        "The bounded number of candidate samples is not chosen until the distance margin is known.",
    ]
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(json.dumps(payload["constituent_event"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
