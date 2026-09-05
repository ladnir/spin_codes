#!/usr/bin/env python3
"""Verify the three binary-projection bounds for the pair accumulator kernel."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from analyze_accumulator_pair_chain_small import compositions4, rank_two
from probe_triangle_holder_pair_bound_small import triple_weights
from verify_accumulator_pair_type_kernel import kernel_counts, multinomial


WORKSTREAM = Path(__file__).resolve().parent
OUTPUT = WORKSTREAM / "pair_projection_kernel_B8_exact.json"


def binary_transition(length: int, input_weight: int, output_weight: int) -> Fraction:
    if input_weight == 0:
        return Fraction(int(output_weight == 0))
    half_down = input_weight // 2
    half_up = (input_weight + 1) // 2
    if half_down > length - output_weight or not (1 <= half_up <= output_weight):
        return Fraction(0)
    return Fraction(
        math.comb(length - output_weight, half_down)
        * math.comb(output_weight - 1, half_up - 1),
        math.comb(length, input_weight),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, choices=(8, 16), default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    types = [kind for kind in compositions4(args.length) if rank_two(kind)]
    maximum_ratio = [Fraction(0), Fraction(0), Fraction(0)]
    witnesses: list[dict[str, object] | None] = [None, None, None]
    checked = 0
    for input_type in types:
        input_mass = multinomial(input_type)
        input_weights = triple_weights(input_type)
        for output_type, count in kernel_counts(input_type).items():
            output_weights = triple_weights(output_type)
            pair_probability = Fraction(count, input_mass)
            for omitted in range(3):
                selected = [index for index in range(3) if index != omitted]
                first, second = selected
                overlap_probability = Fraction(
                    input_mass,
                    math.comb(args.length, input_weights[first])
                    * math.comb(args.length, input_weights[second]),
                )
                product = binary_transition(
                    args.length,
                    input_weights[first],
                    output_weights[first],
                ) * binary_transition(
                    args.length,
                    input_weights[second],
                    output_weights[second],
                )
                left = overlap_probability * pair_probability
                if left > product:
                    raise AssertionError(
                        f"projection inequality failed for {input_type} -> {output_type}"
                    )
                if product:
                    ratio = left / product
                    if ratio > maximum_ratio[omitted]:
                        maximum_ratio[omitted] = ratio
                        witnesses[omitted] = {
                            "input_type": list(input_type),
                            "output_type": list(output_type),
                        }
                checked += 1

    payload = {
        "schema": "pair-projection-kernel-exact-v1",
        "status": "EXACT_INTEGER_RATIONAL_VERIFICATION",
        "parameters": {
            "length": args.length,
            "rank_two_input_types": len(types),
            "projection_inequalities_checked": checked,
        },
        "claim": "r_ij(n) P(n,m) <= Q(h_i,w_i) Q(h_j,w_j) for all three binary-character pairs",
        "maximum_ratios": [
            {
                "omitted_character": omitted,
                "numerator": ratio.numerator,
                "denominator": ratio.denominator,
                "witness": witnesses[omitted],
            }
            for omitted, ratio in enumerate(maximum_ratio)
        ],
        "consequence": "Taking the geometric mean gives the three-projection pair-kernel majorant used by probe_pair_projection_holder_bound_small.py.",
        "limitations": [
            "The finite verification is at a small length; the conditioning argument proving the inequality is length independent.",
            "The inequality is valid but its iterated factored envelope is too loose for the BA variance target.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["maximum_ratios"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
