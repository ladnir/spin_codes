#!/usr/bin/env python3
"""Outward Hilbert--Schmidt certificate for the one-word accumulator.

For nonzero input and output weights h,w, the stationary-weighted transition
matrix has the exact entry

    S[h,w] = N_B(h,w) / sqrt(binomial(B,h) binomial(B,w)),

where N_B is the terminated accumulator input-output enumerator.  The script
uses Arb to sum the exact squared entries of S.  After subtracting the
constant singular direction, this centered Hilbert--Schmidt norm bounds every
nonconstant singular value.  The bound is less sharp than a full eigensolve,
but it is fast and verifier-friendly.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from flint import arb, ctx


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "accumulator_oneword_singular_B512_outward.json"


def accumulator_count(length: int, input_weight: int, output_weight: int) -> int:
    lower = input_weight // 2
    upper = (input_weight + 1) // 2
    if lower > length - output_weight or upper < 1 or upper - 1 > output_weight - 1:
        return 0
    return math.comb(length - output_weight, lower) * math.comb(
        output_weight - 1, upper - 1
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, default=512)
    parser.add_argument("--singular-numerator", type=int, default=11)
    parser.add_argument("--singular-denominator", type=int, default=100)
    parser.add_argument("--precision-bits", type=int, default=256)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.length < 2:
        parser.error("length must be at least two")
    if args.singular_numerator <= 0 or args.singular_denominator <= 0:
        parser.error("singular bound must be positive")

    ctx.prec = args.precision_bits
    length = args.length
    shell_sizes = [math.comb(length, weight) for weight in range(1, length + 1)]
    squared_hilbert_schmidt = arb(0)
    for row, input_weight in enumerate(range(1, length + 1)):
        for column, output_weight in enumerate(range(1, length + 1)):
            count = accumulator_count(length, input_weight, output_weight)
            if count:
                squared_hilbert_schmidt += arb(count * count) / (
                    arb(shell_sizes[row]) * shell_sizes[column]
                )

    # The stationary constant direction has singular value one.  It is
    # orthogonal to every nonconstant singular direction, so removing its
    # squared contribution gives the centered Hilbert--Schmidt norm.
    centered_squared_hilbert_schmidt = squared_hilbert_schmidt - 1
    threshold = (arb(args.singular_numerator) / args.singular_denominator) ** 2
    if not bool(centered_squared_hilbert_schmidt.lower() > 0):
        raise ArithmeticError("centered Hilbert--Schmidt square is not positive")
    if not bool(centered_squared_hilbert_schmidt.upper() < threshold.lower()):
        raise AssertionError("Hilbert--Schmidt bound does not meet the requested square")

    payload = {
        "schema": "accumulator-oneword-hilbert-schmidt-outward-v1",
        "status": "OUTWARD_SINGULAR_VALUE_CERTIFICATE",
        "claim": {
            "length": length,
            "second_singular_value_strictly_below": (
                f"{args.singular_numerator}/{args.singular_denominator}"
            ),
            "centered_squared_hilbert_schmidt_upper": str(
                centered_squared_hilbert_schmidt.upper()
            ),
            "requested_square_lower": str(threshold.lower()),
        },
        "operator": {
            "state_space": "nonzero Hamming weights 1 through B",
            "stationary_measure": "binomial(B,h)/(2^B-1)",
            "stage": "uniform coordinate permutation followed by the zero-state prefix accumulator",
            "centering": "subtract the stationary constant singular direction from the squared Hilbert--Schmidt norm",
        },
        "parameters": {
            "precision_bits": args.precision_bits,
            "arithmetic": "exact integer counts and 256-bit Arb interval summation",
        },
        "limitations": [
            "This certificate concerns the one-word weight operator only.",
            "The ordered-pair interaction-sector bound remains separate.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
