#!/usr/bin/env python3
"""Evaluate occupation Q=2 for a fixed repeated constituent spectrum.

Two active outer rows select codewords independently from the same fixed
linear constituent.  Independent row-coordinate permutations make their
supports independent conditional on their weights.  Therefore the ordinary
weight spectrum determines the complete Q=2 message sum; no pair enumerator
of the constituent is required.

The calculation uses exact uniform-support transfer coefficients in the proof
model and nearest binary64 log arithmetic.  It is diagnostic.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

import evaluate_ebch128_randomstepconv_g1 as inner  # noqa: E402
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients  # noqa: E402
from evaluate_single_random_constituent_q2 import pair_support_coefficients  # noqa: E402
from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)


LOG2 = math.log(2.0)
DEFAULT_OUTPUT = HERE / "fixed_spectrum_q2_diagnostic.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constituent", choices=sorted(CONSTITUENTS), default="rm49")
    parser.add_argument("--message-exponent", type=int, default=11)
    parser.add_argument("--memory-bits", type=int, default=18)
    parser.add_argument("--log-surprisal", type=float, nargs="+", default=[-3.0, -2.5, -2.0])
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    constituent = CONSTITUENTS[args.constituent]
    spectrum = load_spectrum(constituent)
    message_bits = 1 << args.message_exponent
    if message_bits % constituent.dimension:
        parser.error("constituent dimension must divide total message bits")
    outer_rows = message_bits // constituent.dimension
    if outer_rows < 2:
        parser.error("occupation Q=2 requires at least two outer rows")
    output_bits = 2 * message_bits
    distance_cutoff = (
        args.distance_numerator * output_bits + args.distance_denominator - 1
    ) // args.distance_denominator
    block_bits = constituent.block_bits
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

    weights = sorted(weight for weight, count in spectrum.items() if weight and count)
    terms = []
    rows = []
    log_row_pairs = math.log(math.comb(outer_rows, 2))
    for first in weights:
        for second in weights:
            value = (
                log_row_pairs
                + math.log(spectrum[first])
                + math.log(spectrum[second])
                + float(best[first, second])
            )
            terms.append(value)
            rows.append((value, first, second, float(witnesses[first, second])))
    aggregate_log2 = float(logsumexp(terms) / LOG2)
    rows.sort(reverse=True)
    payload = {
        "schema": "fixed-repeated-exact-spectrum-randomstepconv-q2-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "constituent": constituent.name,
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
            "outer": "one fixed deterministic constituent, repeated in every outer row",
            "routing": "independent uniform coordinate permutation per row and independent uniform position permutation per transposed region",
            "inner": "independent RandomStepConv map at every position, sampled once and shared by all messages",
        },
        "claim": {
            "q2_log2_upper_diagnostic": aggregate_log2,
            "q2_margin_bits_diagnostic": -aggregate_log2,
            "q2_closes_40_bits_diagnostic": aggregate_log2 < -40,
            "dominant_pairs": [
                {
                    "first_weight": first,
                    "second_weight": second,
                    "pointwise_log2": value / LOG2,
                    "log_surprisal": witness,
                }
                for value, first, second, witness in rows[:20]
            ],
        },
        "limitations": [
            "Nearest binary64 log arithmetic is not outward rounded.",
            "Only occupation Q=2 is covered.",
            "The finite witness list is not asserted optimal.",
            "The inner is RandomStepConv, not RM2Sub.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
