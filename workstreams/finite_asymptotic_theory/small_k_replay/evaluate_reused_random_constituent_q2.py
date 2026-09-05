#!/usr/bin/env python3
"""Evaluate Q=2 for one random full-rank constituent reused across rows.

Two nonzero local messages are either equal or distinct.  Equal messages have
one shared constituent image.  Distinct binary messages are linearly
independent, so a uniform full-rank constituent maps them to a uniform ordered
pair of distinct nonzero output vectors.  These two exact laws account for
constituent reuse at occupation two.

The calculation uses nearest binary64 arithmetic and is diagnostic.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import gammaln, logsumexp


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

import evaluate_ebch128_randomstepconv_g1 as inner  # noqa: E402
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients  # noqa: E402
from evaluate_single_random_constituent_q2 import pair_support_coefficients  # noqa: E402


LOG2 = math.log(2.0)
DEFAULT_OUTPUT = HERE / "reused_random_constituent_q2_diagnostic.json"


def log_nonzero_count(bits: int) -> float:
    return bits * LOG2 + math.log1p(-math.ldexp(1.0, -bits))


def log_binomial_rows(bits: int) -> np.ndarray:
    weights = np.arange(bits + 1, dtype=np.float64)
    return (
        gammaln(bits + 1)
        - gammaln(weights + 1)
        - gammaln(bits - weights + 1)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-bits", type=int, default=512)
    parser.add_argument("--message-exponent", type=int, default=13)
    parser.add_argument("--memory-bits", type=int, default=12)
    parser.add_argument(
        "--log-surprisal", type=float, nargs="+", default=[-3.5, -3.0, -2.5]
    )
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    block_bits = args.block_bits
    if block_bits <= 0 or block_bits % 2:
        parser.error("block length must be a positive even integer")
    dimension = block_bits // 2
    message_bits = 1 << args.message_exponent
    if message_bits % dimension:
        parser.error("constituent dimension must divide total message bits")
    outer_rows = message_bits // dimension
    if outer_rows < 2:
        parser.error("occupation Q=2 requires at least two outer rows")
    output_bits = 2 * message_bits
    distance_cutoff = (
        args.distance_numerator * output_bits + args.distance_denominator - 1
    ) // args.distance_denominator

    best = np.full((block_bits + 1, block_bits + 1), math.inf)
    witnesses = np.full_like(best, math.nan)
    for index, log_surprisal in enumerate(args.log_surprisal):
        print(
            f"tilt,{index + 1},{len(args.log_surprisal)},u,{log_surprisal}",
            flush=True,
        )
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        zero, active = inner.step_matrices(z, args.memory_bits)
        regions = uniform_coefficients(
            inner.log_entries(zero),
            inner.log_entries(active),
            outer_rows,
            2,
        )
        pairs = pair_support_coefficients(regions, block_bits)
        moments = np.logaddexp(pairs[..., 0, 0], pairs[..., 0, 1])
        values = np.minimum(0.0, moments + distance_cutoff * surprisal)
        improved = values < best
        best[improved] = values[improved]
        witnesses[improved] = log_surprisal

    log_a = log_nonzero_count(dimension)
    log_a_minus_one = math.log((1 << dimension) - 2)
    log_v = log_nonzero_count(block_bits)
    log_v_minus_one = math.log((1 << block_bits) - 2)
    log_bins = log_binomial_rows(block_bits)
    log_row_pairs = math.log(math.comb(outer_rows, 2))

    equal_terms = []
    distinct_terms = []
    dominant = []
    for first in range(1, block_bits + 1):
        equal_value = (
            log_row_pairs
            + log_a
            + float(log_bins[first])
            - log_v
            + float(best[first, first])
        )
        equal_terms.append(equal_value)
        dominant.append((equal_value, "equal", first, first, witnesses[first, first]))
        for second in range(1, block_bits + 1):
            log_pair_count = float(log_bins[first] + log_bins[second])
            if first == second:
                # Remove the forbidden pair in which the two output vectors
                # are identical.  At weights 0 and B the shell has size one.
                if log_bins[first] == 0.0:
                    continue
                log_pair_count += math.log1p(-math.exp(-float(log_bins[first])))
            value = (
                log_row_pairs
                + log_a
                + log_a_minus_one
                + log_pair_count
                - log_v
                - log_v_minus_one
                + float(best[first, second])
            )
            distinct_terms.append(value)
            dominant.append((value, "distinct", first, second, witnesses[first, second]))

    equal_log = float(logsumexp(equal_terms))
    distinct_log = float(logsumexp(distinct_terms))
    total_log = float(np.logaddexp(equal_log, distinct_log))
    dominant.sort(reverse=True)
    payload = {
        "schema": "reused-random-full-rank-constituent-q2-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "constituent": f"random full-rank [{block_bits},{dimension}]",
            "message_exponent": args.message_exponent,
            "message_bits": message_bits,
            "output_bits": output_bits,
            "outer_rows": outer_rows,
            "memory_bits": args.memory_bits,
            "distance_cutoff": distance_cutoff,
            "relative_distance_lower": distance_cutoff / output_bits,
            "log_surprisal_witnesses": args.log_surprisal,
        },
        "probability_space": {
            "outer": "one uniform full-rank constituent sampled once and reused in every row",
            "routing": "independent uniform row-coordinate and transposed-region permutations",
            "inner": "independent RandomStepConv maps sampled once and shared by all messages",
        },
        "claim": {
            "equal_local_message_log2_upper_diagnostic": equal_log / LOG2,
            "distinct_local_message_log2_upper_diagnostic": distinct_log / LOG2,
            "q2_log2_upper_diagnostic": total_log / LOG2,
            "q2_margin_bits_diagnostic": -total_log / LOG2,
            "q2_closes_40_bits_diagnostic": total_log < -40 * LOG2,
            "dominant_terms": [
                {
                    "sector": sector,
                    "first_weight": first,
                    "second_weight": second,
                    "pointwise_log2": value / LOG2,
                    "log_surprisal": float(witness),
                }
                for value, sector, first, second, witness in dominant[:20]
            ],
        },
        "proof_reduction": [
            "Split ordered pairs of nonzero local messages into equal and distinct sectors.",
            "An equal pair has one uniform nonzero output vector under the sampled full-rank map.",
            "A distinct binary pair is linearly independent and has a uniform ordered pair of distinct nonzero output vectors.",
            "Average the exact shell-pair law against independent row permutations and the exact Q=2 regional transfer.",
        ],
        "limitations": [
            "Nearest binary64 log arithmetic is not outward rounded.",
            "Only occupation Q=2 is covered.",
            "The finite witness list is not asserted optimal.",
            "For Q>=3, the rank and linear-relation type of the active local-message tuple must also be tracked.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
