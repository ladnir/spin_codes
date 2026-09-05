#!/usr/bin/env python3
"""First-moment diagnostic for an exact-size block EA or EC constituent.

The default [512,256] ensemble uses fourteen disjoint output regions: eight
of length 37 and six of length 36.  Every input coordinate independently
chooses one neighbor in every region.  An accumulator or memory-15 wrapped
convolution follows the sparse map.  The calculation is binary64.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


DEFAULT_EXPANDER_ROOT = Path(
    r"C:\Users\peter\.codex\worktrees\30ff\permute_conv\expander_codes"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expander-root", type=Path, default=DEFAULT_EXPANDER_ROOT)
    parser.add_argument("--code", choices=("ea", "ec"), default="ec")
    parser.add_argument("--k", type=int, default=256)
    parser.add_argument("--memory", type=int, default=15)
    parser.add_argument(
        "--region-lengths",
        type=int,
        nargs="+",
        default=(37,) * 8 + (36,) * 6,
    )
    parser.add_argument(
        "--markers",
        type=float,
        nargs="+",
        default=(0.03, 0.05, 0.1, 0.2, 0.5),
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def scaled_product_log_moment(
    matrices: list[np.ndarray], start_state: int
) -> float:
    product = np.eye(matrices[0].shape[0], dtype=np.float64)
    log_scale = 0.0
    for matrix in matrices:
        product = product @ matrix
        scale = float(np.max(product))
        if scale == 0.0:
            return -math.inf
        product /= scale
        log_scale += math.log(scale)
    endpoint = float(np.sum(product[start_state, :]))
    return log_scale + math.log(endpoint)


def main() -> None:
    args = parse_args()
    scripts = args.expander_root / "scripts"
    sys.path.insert(0, str(scripts))
    from regular_ec_diagnostic import (  # pylint: disable=import-error,import-outside-toplevel
        regional_shell_distribution,
        uniform_slice_transfer_matrices,
    )
    from regular_ea_diagnostic import (  # pylint: disable=import-error,import-outside-toplevel
        accumulator_uniform_slice_transfers,
    )

    k = args.k
    n = sum(args.region_lengths)
    unique_lengths = sorted(set(args.region_lengths))
    shell_cache = {
        (length, support): regional_shell_distribution(length, support)
        for length in unique_lengths
        for support in range(1, k + 1)
    }

    kernel_support_logs: list[float] = []
    for support in range(1, k + 1):
        log_probability = 0.0
        possible = True
        for length in args.region_lengths:
            probability = shell_cache[(length, support)].get(0, 0.0)
            if probability == 0.0:
                possible = False
                break
            log_probability += math.log(probability)
        if possible:
            kernel_support_logs.append(
                math.lgamma(k + 1)
                - math.lgamma(support + 1)
                - math.lgamma(k - support + 1)
                + log_probability
            )
    kernel_log2 = float(logsumexp(kernel_support_logs) / math.log(2.0))

    rows: list[dict[str, float]] = []
    for marker in args.markers:
        if args.code == "ec":
            slices_by_length = {
                length: uniform_slice_transfer_matrices(
                    length=length,
                    max_weight=length,
                    output_marker=marker,
                    memory=args.memory,
                )
                for length in unique_lengths
            }
            start_state = args.memory
        else:
            slices_by_length = {
                length: accumulator_uniform_slice_transfers(
                    length=length,
                    max_weight=length,
                    output_marker=marker,
                )
                for length in unique_lengths
            }
            start_state = 0
        support_logs: list[float] = []
        for support in range(1, k + 1):
            matrices: list[np.ndarray] = []
            for length in args.region_lengths:
                slices = slices_by_length[length]
                region = np.zeros_like(slices[0])
                for weight, probability in shell_cache[(length, support)].items():
                    region += probability * slices[weight]
                matrices.append(region)
            log_moment = scaled_product_log_moment(matrices, start_state)
            support_logs.append(
                math.lgamma(k + 1)
                - math.lgamma(support + 1)
                - math.lgamma(k - support + 1)
                + log_moment
            )
        ec_log2 = float(logsumexp(support_logs) / math.log(2.0))
        random_log2 = (
            math.log2(2.0**k - 1.0)
            + n * (math.log2(1.0 + marker) - 1.0)
        )
        worst_support = int(np.argmax(np.asarray(support_logs))) + 1
        positive_log2 = math.nan
        if ec_log2 > kernel_log2:
            positive_log2 = kernel_log2 + math.log2(
                2.0 ** (ec_log2 - kernel_log2) - 1.0
            )
        rows.append(
            {
                "marker": marker,
                "ensemble_nonzero_message_mgf_log2": ec_log2,
                "ensemble_positive_output_mgf_log2": positive_log2,
                "random_linear_nonzero_mgf_log2": random_log2,
                "ensemble_minus_random_bits": ec_log2 - random_log2,
                "largest_input_support": worst_support,
                "largest_input_support_log2": support_logs[worst_support - 1]
                / math.log(2.0),
            }
        )

    payload = {
        "schema": "block-expand-recursive-balanced-left-regular-mgf-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "message_bits": k,
            "output_bits": n,
            "recursive_map": args.code,
            "left_degree": len(args.region_lengths),
            "region_lengths_in_order": args.region_lengths,
            "memory": args.memory if args.code == "ec" else None,
        },
        "expected_nonzero_kernel_words_log2": kernel_log2,
        "kernel_failure_probability_markov_bits": -kernel_log2,
        "rows": rows,
        "limitations": [
            "The calculation is binary64 and is not an outward certificate.",
            "A first moment does not prove concentration of the positive spectrum of one reused constituent.",
            "The kernel-word expectation gives a valid Markov target only after outward recomputation.",
        ],
    }
    encoded = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
