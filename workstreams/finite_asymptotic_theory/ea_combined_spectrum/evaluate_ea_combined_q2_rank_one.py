#!/usr/bin/env python3
"""Evaluate the equal-message (rank-one) part of sparse-EA occupation Q=2.

Two distinct active outer rows may carry the same nonzero local message.  The
same sampled constituent then produces the same outer word in both rows, but
the two row-coordinate permutations are independent.  This program averages
that exact rank-one contribution over one sparse-EA constituent.

The result is only a necessary Q=2 gate: distinct local messages form the
rank-two contribution and are not evaluated here.  Arithmetic is binary64 or
long double without directed rounding.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

import evaluate_ebch128_randomstepconv_g1 as transfer  # noqa: E402
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients  # noqa: E402
from ea_combined_spectrum.evaluate_ea_combined_q1 import (  # noqa: E402
    expected_accumulated_spectrum_log2_truncated,
)


DEFAULT_OUTPUT = HERE / "ea_K256_B512_combined_q2_rank_one_degree_sweep.json"


def pair_support_coefficients_truncated(
    region_matrices: np.ndarray,
    block_bits: int,
    maximum_support: int,
) -> np.ndarray:
    """Average two independent supports, retaining degrees <= maximum_support."""
    current = np.full(
        (maximum_support + 1, maximum_support + 1, 2, 2), -math.inf
    )
    current[0, 0] = transfer.log_identity()
    for completed in range(block_bits):
        size = completed + 1
        old_maximum = min(completed, maximum_support)
        new_maximum = min(size, maximum_support)
        old = current[: old_maximum + 1, : old_maximum + 1]
        updated = np.full(
            (new_maximum + 1, new_maximum + 1, 2, 2), -math.inf
        )
        old_degrees = np.arange(old_maximum + 1, dtype=np.float64)
        unselected = np.log((size - old_degrees) / float(size))

        term = transfer.log_matmul(old, region_matrices[0])
        term += unselected[:, None, None, None]
        term += unselected[None, :, None, None]
        updated[: old_maximum + 1, : old_maximum + 1] = term

        selected = np.log(
            np.arange(1, new_maximum + 1, dtype=np.float64) / float(size)
        )
        first_old = current[:new_maximum, : old_maximum + 1]
        term = transfer.log_matmul(first_old, region_matrices[1])
        term += selected[:, None, None, None]
        term += unselected[None, :, None, None]
        updated[1 : new_maximum + 1, : old_maximum + 1] = np.logaddexp(
            updated[1 : new_maximum + 1, : old_maximum + 1], term
        )

        second_old = current[: old_maximum + 1, :new_maximum]
        term = transfer.log_matmul(second_old, region_matrices[1])
        term += unselected[:, None, None, None]
        term += selected[None, :, None, None]
        updated[: old_maximum + 1, 1 : new_maximum + 1] = np.logaddexp(
            updated[: old_maximum + 1, 1 : new_maximum + 1], term
        )

        both_old = current[:new_maximum, :new_maximum]
        term = transfer.log_matmul(both_old, region_matrices[2])
        term += selected[:, None, None, None]
        term += selected[None, :, None, None]
        updated[1 : new_maximum + 1, 1 : new_maximum + 1] = np.logaddexp(
            updated[1 : new_maximum + 1, 1 : new_maximum + 1], term
        )
        current[: new_maximum + 1, : new_maximum + 1] = updated
    return current


def evaluate_transfer(
    *,
    block_bits: int,
    outer_rows: int,
    output_bits: int,
    distance_cutoff: int,
    memory_bits: int,
    log_surprisal: float,
    maximum_support: int,
) -> np.ndarray:
    surprisal = math.exp(log_surprisal)
    z = math.exp(-surprisal)
    zero, active = transfer.step_matrices(z, memory_bits)
    regions = uniform_coefficients(
        transfer.log_entries(zero),
        transfer.log_entries(active),
        outer_rows,
        2,
    )
    pairs = pair_support_coefficients_truncated(
        regions, block_bits, maximum_support
    )
    moments = np.logaddexp(pairs[..., 0, 0], pairs[..., 0, 1])
    return np.minimum(0.0, moments + distance_cutoff * surprisal)


def bits(value: np.longdouble) -> float | None:
    if not np.isfinite(value) or value <= 0:
        return None
    return float(-np.log2(value))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message-bits", type=int, default=256)
    parser.add_argument("--block-bits", type=int, default=512)
    parser.add_argument("--outer-rows", type=int, default=4096)
    parser.add_argument("--output-bits", type=int, default=1 << 21)
    parser.add_argument("--distance-cutoff", type=int, default=228590)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--log-surprisal", type=float, default=-7.5)
    parser.add_argument("--max-support", type=int)
    parser.add_argument("--rank-attempts", type=int, default=16)
    parser.add_argument("--degrees", type=int, nargs="+", default=[19, 21, 23, 25, 29, 33])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.outer_rows * args.block_bits != args.output_bits:
        parser.error("outer_rows * block_bits must equal output_bits")

    maximum_support = args.block_bits if args.max_support is None else args.max_support
    if not 1 <= maximum_support <= args.block_bits:
        parser.error("max-support must lie in [1,block-bits]")
    print("transfer", flush=True)
    q2_transfer = evaluate_transfer(
        block_bits=args.block_bits,
        outer_rows=args.outer_rows,
        output_bits=args.output_bits,
        distance_cutoff=args.distance_cutoff,
        memory_bits=args.memory_bits,
        log_surprisal=args.log_surprisal,
        maximum_support=maximum_support,
    )
    diagonal_log = np.asarray(
        [q2_transfer[weight, weight] for weight in range(maximum_support + 1)],
        dtype=np.float64,
    )
    row_pairs = math.comb(args.outer_rows, 2)
    rows = []
    for degree in args.degrees:
        print(f"degree,{degree}", flush=True)
        spectrum_log2 = expected_accumulated_spectrum_log2_truncated(
            args.message_bits, args.block_bits, degree, maximum_support
        )
        kernel_log = spectrum_log2[0] * math.log(2.0)
        kernel_mean = np.longdouble(math.exp(kernel_log))
        rank_success_lower = max(np.longdouble(0), np.longdouble(1) - kernel_mean)
        abort = min(np.longdouble(1), kernel_mean) ** args.rank_attempts
        log_terms = (
            spectrum_log2[1 : maximum_support + 1] * math.log(2.0)
            + diagonal_log[1:]
        )
        maximum_log_term = float(np.max(log_terms))
        log_sum = maximum_log_term + math.log(
            float(np.sum(np.exp(log_terms - maximum_log_term)))
        )
        log_unconditional = math.log(row_pairs) + log_sum
        unconditional = np.longdouble(math.exp(log_unconditional))
        conditional = (
            np.longdouble(np.inf)
            if rank_success_lower == 0
            else unconditional / rank_success_lower
        )
        total = abort + conditional
        dominant = int(np.argmax(log_terms)) + 1
        rows.append({
            "right_degree": degree,
            "total_xors_per_constituent": args.block_bits * degree - 1,
            "unconditional_kernel_expected_count": float(kernel_mean),
            "rank_test_abort_upper": float(abort),
            "rank_one_q2_unconditional_upper": float(unconditional),
            "rank_one_q2_conditional_upper": float(conditional),
            "rank_one_q2_conditional_margin_bits": bits(conditional),
            "setup_abort_or_rank_one_q2_margin_bits": bits(total),
            "dominant_outer_weight": dominant,
            "dominant_weighted_term_log2": (
                float(log_terms[dominant - 1] / math.log(2.0))
            ),
        })

    payload = {
        "schema": "ea-transfer-weighted-combined-q2-rank-one-v1",
        "status": "BINARY_FLOAT_DIAGNOSTIC",
        "parameters": {
            "message_bits": args.message_bits,
            "block_bits": args.block_bits,
            "outer_rows": args.outer_rows,
            "output_bits": args.output_bits,
            "distance_cutoff": args.distance_cutoff,
            "memory_bits": args.memory_bits,
            "log_surprisal": args.log_surprisal,
            "maximum_support_included": maximum_support,
            "rank_attempts": args.rank_attempts,
        },
        "probability_space": (
            "one sparse-EA constituent is rank-tested and reused; two active "
            "outer rows carry the same nonzero local message; their row "
            "permutations are independent; region permutations and "
            "RandomStepConv maps are averaged by the Q=2 transfer"
        ),
        "rows": rows,
        "limitations": [
            "This covers only the equal-message rank-one sector of occupation Q=2.",
            "The distinct-message rank-two sector remains open.",
            (
                "All rank-one outer weights are included."
                if maximum_support == args.block_bits
                else "Only outer weights through maximum_support are included; the displayed sum is a partial value of this Chernoff proof bound, so it can rule out closure by this witness but cannot prove closure."
            ),
            "Nearest binary64 and long-double arithmetic is not outward rounded.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(rows, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
