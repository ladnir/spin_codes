#!/usr/bin/env python3
"""Probe a weight-level covariance-block norm bound at small lengths.

The bound uses exact level energies and binary64 covariance-block entries.
It is a diagnostic for choosing the length-512 certificate interface.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from verify_pair_kernel_small import analytic_moments
from verify_subset_covariance_blocks_small import (
    covariance_blocks,
    covariance_kernel,
    shell_coefficient,
)


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "subset_block_norm_bound_K4_B8_r3_w6_probe.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=4)
    parser.add_argument("--output-bits", type=int, default=8)
    parser.add_argument("--right-degree", type=int, default=3)
    parser.add_argument("--shell-weight", type=int, default=6)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output_bits > 16:
        parser.error("this subset-enumerating probe is limited to output_bits <= 16")

    _returns, kernel = covariance_kernel(
        args.message_bits, args.output_bits, args.right_degree
    )
    blocks = covariance_blocks(args.output_bits, kernel)
    level_energy = np.zeros(args.output_bits + 1)
    for subset in range(1 << args.output_bits):
        coefficient = shell_coefficient(
            args.output_bits, args.shell_weight, subset
        )
        level_energy[subset.bit_count()] += coefficient * coefficient

    entry_norm = np.zeros((args.output_bits + 1, args.output_bits + 1))
    global_norm = 0.0
    for block in blocks:
        matrix = block["matrix"]
        levels = block["levels"]
        global_norm = max(global_norm, float(np.linalg.eigvalsh(matrix)[-1]))
        for row, first_weight in enumerate(levels):
            for column, second_weight in enumerate(levels):
                entry_norm[first_weight, second_weight] = max(
                    entry_norm[first_weight, second_weight],
                    abs(float(matrix[row, column])),
                )

    level_norms = np.sqrt(level_energy)
    level_bound = float(level_norms @ entry_norm @ level_norms)
    global_bound = float(global_norm * level_energy.sum())

    mean, factorial = analytic_moments(
        args.message_bits, args.output_bits, args.right_degree
    )
    exact_mean = mean[args.shell_weight]
    exact_variance = (
        exact_mean
        + factorial[args.shell_weight]
        - exact_mean * exact_mean
    )
    scale = 2.0 ** (2 * (args.message_bits - args.output_bits))
    exact_ratio = float(exact_variance / exact_mean)
    level_ratio_bound = scale * level_bound / float(exact_mean)
    global_ratio_bound = scale * global_bound / float(exact_mean)

    payload = {
        "schema": "pure-ea-subset-block-norm-bound-small-v1",
        "status": "EXACT_KERNEL_WITH_BINARY64_OPERATOR_BOUND_DIAGNOSTIC",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": args.output_bits,
            "right_degree": args.right_degree,
            "shell_weight": args.shell_weight,
        },
        "exact_variance_to_mean": exact_ratio,
        "weight_level_block_variance_to_mean_upper": level_ratio_bound,
        "global_operator_variance_to_mean_upper": global_ratio_bound,
        "weight_level_loss_over_exact": level_ratio_bound / exact_ratio,
        "global_loss_over_exact": global_ratio_bound / exact_ratio,
        "scope": [
            "The covariance kernel and level energies arise from exact finite formulas.",
            "The block eigensystems and displayed bounds use nondirected binary64 arithmetic.",
            "The result is not an outward certificate and does not extrapolate in length.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
