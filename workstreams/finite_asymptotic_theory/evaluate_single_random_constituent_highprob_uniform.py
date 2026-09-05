#!/usr/bin/env python3
"""All-Q Bernoulli-half envelope implied by one-shot spectrum caps.

This deliberately coarse diagnostic is intended for the dense occupations.
Sparse occupations are handled by exact-shell and three-band calculations.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import evaluate_repeated_random512_randomstepconv_g1 as base
from evaluate_single_random_constituent_highprob_renyi import spectrum_caps


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_OUTPUT = WORKSTREAM / "single_random_constituent_B512_highprob_uniform_s22.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-bits", type=int, default=512)
    parser.add_argument("--dimension", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=1 << 21)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--per-shell-failure-exponent", type=int, default=51)
    parser.add_argument("--grid-min", type=float, default=-9.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--refine-min", type=float, default=0.65)
    parser.add_argument("--refine-max", type=float, default=0.80)
    parser.add_argument("--refine-step", type=float, default=0.025)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output_bits % args.block_bits:
        parser.error("block length must divide output length")

    caps, event = spectrum_caps(
        math.ldexp(1.0, -args.per_shell_failure_exponent),
        block_bits=args.block_bits,
        dimension=args.dimension,
    )
    message_count = (1 << args.dimension) - 1
    factor = max(
        int(caps[weight])
        / (message_count * math.comb(args.block_bits, weight) / (1 << args.block_bits))
        for weight in range(1, args.block_bits + 1)
        if int(caps[weight])
    )

    base.B = args.block_bits
    base.K = args.dimension
    base.L = args.output_bits // args.block_bits
    base.ACTIVE_ROWS = base.L
    base.N = args.output_bits
    base.D = (
        args.distance_numerator * args.output_bits
        + args.distance_denominator
        - 1
    ) // args.distance_denominator
    base.MEMORY = args.memory_bits
    base.SPECTRUM_FACTOR = factor
    payload = base.evaluate(args)
    payload["schema"] = "single-random-constituent-highprob-uniform-v1"
    payload["status"] = "BINARY64_DENSE_DIAGNOSTIC"
    payload["parameters"]["outer_code"] = (
        f"one uniform binary [{args.block_bits},{args.dimension}] generator"
    )
    payload["parameters"]["spectrum_factor"] = factor
    payload["parameters"]["spectrum_factor_log2"] = math.log2(factor)
    payload["constituent_event"] = event
    payload["limitations"] = [
        "Nearest binary64 arithmetic is not an outward certificate.",
        "The single Bernoulli-half envelope is intentionally loose in the tails.",
        "A final proof may use this receipt only for a suffix of the occupations.",
    ]
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"spectrum_factor_log2={math.log2(factor):.12f}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
