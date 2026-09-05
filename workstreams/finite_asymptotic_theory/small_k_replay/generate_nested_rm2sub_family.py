#!/usr/bin/env python3
"""Generate one exactly audited nested RM2Sub family.

For t=2^m, the family fixes one ordered list of RM(2,m) generators.  The map
A_s uses the constant word, every linear word, and the first remaining
quadratic words.  Thus A_s is a literal subcode of A_{s+1}.  The corresponding
B_s is the transpose of the same generator matrix, so B_s A_s=0.

The script samples one ordered quadratic basis from a declared seed.  It does
not search or select among candidates.  This makes changes in s causal within
the generated chain, but it does not make one chain representative of the
sampling distribution.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

from generate_rm2sub_calibration_constituent import (
    coordinate_columns,
    enumerate_spectrum,
    generator_words,
    kernel_weight_four_count,
    macwilliams_kernel,
    rank,
    sample_independent_masks,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step-bits", type=int, default=128)
    parser.add_argument("--minimum-state-bits", type=int, default=12)
    parser.add_argument("--maximum-state-bits", type=int, default=16)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.step_bits < 32 or args.step_bits & (args.step_bits - 1):
        parser.error("step size must be a power of two at least 32")
    variables = args.step_bits.bit_length() - 1
    base_dimension = 1 + variables
    if args.minimum_state_bits < base_dimension:
        parser.error("minimum state omits the constant or a linear generator")
    monomials = [
        (left, right)
        for left in range(variables)
        for right in range(left + 1, variables)
    ]
    maximum_quadratics = args.maximum_state_bits - base_dimension
    if maximum_quadratics > len(monomials):
        parser.error("maximum state exceeds RM(2,m)")

    quadratic_masks = sample_independent_masks(
        random.Random(args.seed), len(monomials), maximum_quadratics
    )
    full_columns = coordinate_columns(variables, quadratic_masks, monomials)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    members = []
    previous_generators: list[int] | None = None

    for state_bits in range(args.minimum_state_bits, args.maximum_state_bits + 1):
        mask = (1 << state_bits) - 1
        columns = [column & mask for column in full_columns]
        generators = generator_words(columns, state_bits)
        if rank(generators) != state_bits:
            raise ArithmeticError(f"rank failure at s={state_bits}")
        if previous_generators is not None and generators[:-1] != previous_generators:
            raise ArithmeticError(f"prefix nesting failure at s={state_bits}")
        if any((left & right).bit_count() & 1 for left in generators for right in generators):
            raise ArithmeticError(f"BA=0 failure at s={state_bits}")
        spectrum = enumerate_spectrum(generators, state_bits, args.step_bits)
        kernel = macwilliams_kernel(spectrum, state_bits)
        a_distance = next(weight for weight, count in enumerate(spectrum[1:], 1) if count)
        kernel_distance = next(weight for weight, count in enumerate(kernel[1:], 1) if count)
        weight_four = kernel_weight_four_count(columns)
        if kernel[4] != weight_four:
            raise ArithmeticError(f"weight-four audit mismatch at s={state_bits}")

        stem = f"t{args.step_bits}_s{state_bits}"
        selection_path = args.output_dir / f"{stem}_selection.json"
        a_path = args.output_dir / f"{stem}_a_spectrum.json"
        b_path = args.output_dir / f"{stem}_b_kernel_spectrum.json"
        selection = {
            "schema": "rm2sub-nested-constituent-v1",
            "status": "EXACTLY_AUDITED_UNSELECTED_SAMPLE",
            "parameters": {
                "step_bits": args.step_bits,
                "state_bits": state_bits,
                "variables": variables,
                "chain_seed": args.seed,
                "chain_minimum_state_bits": args.minimum_state_bits,
                "chain_maximum_state_bits": args.maximum_state_bits,
            },
            "selected": {
                "seed": args.seed,
                "quadratic_masks_hex": [
                    f"{value:x}" for value in quadratic_masks[: state_bits - base_dimension]
                ],
                "B_columns_hex": [
                    f"{column:0{(state_bits + 3) // 4}x}" for column in columns
                ],
                "A_generator_words_hex": [
                    f"{word:0{args.step_bits // 4}x}" for word in generators
                ],
                "minimum_A_distance": a_distance,
                "minimum_kernel_distance": kernel_distance,
                "weight_four_kernel_words": weight_four,
                "BA_zero": True,
                "prefix_of_next_member": state_bits < args.maximum_state_bits,
            },
            "limitations": [
                "The chain is one deterministic sample from the stated seed.",
                "No candidate search or optimization is performed.",
            ],
        }
        a_payload = {
            "schema": "rm2sub-calibration-a-spectrum-v1",
            "candidate": stem,
            "checks": {
                "spectrum_mass": sum(spectrum),
                "minimum_distance": a_distance,
                "BA_zero": True,
            },
            "spectrum": [
                {"weight": weight, "count": count}
                for weight, count in enumerate(spectrum)
                if count
            ],
        }
        b_payload = {
            "schema": "riffle-splitstate-fixed-b-kernel-spectrum-v1",
            "candidate": stem,
            "parameters": {
                "length": args.step_bits,
                "redundancy": state_bits,
                "kernel_dimension": args.step_bits - state_bits,
            },
            "checks": {
                "dual_mass": sum(spectrum),
                "kernel_mass": sum(kernel),
                "minimum_kernel_distance": kernel_distance,
            },
            "by_total_weight": [
                {
                    "total_weight": weight,
                    "kernel_words": count,
                    "shell_size": math.comb(args.step_bits, weight),
                    "fixed_nonactivation_probability": (
                        count / math.comb(args.step_bits, weight)
                    ),
                    "uniform_support_average_distinct_upper_bound": (
                        count / math.comb(args.step_bits, weight)
                    ),
                }
                for weight, count in enumerate(kernel)
            ],
            "scope": "Exact integer MacWilliams transform of the nested A spectrum.",
        }
        for path, payload in (
            (selection_path, selection),
            (a_path, a_payload),
            (b_path, b_payload),
        ):
            path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        members.append(
            {
                "state_bits": state_bits,
                "a_minimum_distance": a_distance,
                "kernel_minimum_distance": kernel_distance,
                "kernel_weight_four": weight_four,
                "selection": str(selection_path),
                "a_spectrum": str(a_path),
                "kernel_spectrum": str(b_path),
            }
        )
        previous_generators = generators

    manifest = {
        "schema": "rm2sub-nested-family-v1",
        "status": "EXACTLY_AUDITED_UNSELECTED_SAMPLE",
        "construction": (
            "A_s is the prefix of one ordered RM(2,m) generator basis; "
            "B_s is its transpose"
        ),
        "parameters": {
            "step_bits": args.step_bits,
            "minimum_state_bits": args.minimum_state_bits,
            "maximum_state_bits": args.maximum_state_bits,
            "seed": args.seed,
        },
        "members": members,
        "limitations": [
            "Within-chain comparisons isolate state-prefix extension.",
            "Comparisons across chain seeds still include sample variation.",
        ],
    }
    manifest_path = args.output_dir / "NESTED_FAMILY_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"manifest={manifest_path}")
    for member in members:
        print(
            "member,"
            f"s{member['state_bits']},"
            f"dA,{member['a_minimum_distance']},"
            f"dker,{member['kernel_minimum_distance']},"
            f"ker4,{member['kernel_weight_four']}"
        )


if __name__ == "__main__":
    main()
