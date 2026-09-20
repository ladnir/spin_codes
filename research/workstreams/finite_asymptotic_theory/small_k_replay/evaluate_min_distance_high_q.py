#!/usr/bin/env python3
"""Bound Q>=3 by thinning every nonzero outer word to minimum distance.

After an independent uniform row permutation, any fixed subset of d active
coordinates becomes a uniform d-subset.  RandomStepConv's tilted moment is
nonincreasing when input ones are added, so deleting all other ones gives an
upper bound.  A Bernoulli-p row conditioned to weight d is a uniform d-subset.

The program pays the exact conditioning probability for every active row and
uses the exact fixed-Q region coefficient.  Arithmetic is nearest binary64.
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
from small_k_replay.evaluate_exact_spectra_q1_phase import CONSTITUENTS  # noqa: E402


LOG2 = math.log(2.0)
DEFAULT_OUTPUT = HERE / "minimum_distance_high_q_diagnostic.json"


def log_choose(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constituent", choices=sorted(CONSTITUENTS), default="rm49")
    parser.add_argument("--message-exponent", type=int, default=11)
    parser.add_argument("--memory-bits", type=int, default=18)
    parser.add_argument("--minimum-occupation", type=int, default=3)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument(
        "--probabilities",
        type=float,
        nargs="+",
        default=[0.04, 0.05, 0.0625, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3],
    )
    parser.add_argument(
        "--log-surprisals",
        type=float,
        nargs="+",
        default=[-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 0.7, 0.8, 1.0],
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    constituent = CONSTITUENTS[args.constituent]
    message_bits = 1 << args.message_exponent
    if message_bits % constituent.dimension:
        parser.error("constituent dimension must divide total message bits")
    outer_rows = message_bits // constituent.dimension
    if not 1 <= args.minimum_occupation <= outer_rows:
        parser.error("minimum occupation must lie in [1,outer_rows]")
    if any(not 0 < p < 1 for p in args.probabilities):
        parser.error("every Bernoulli probability must lie in (0,1)")
    output_bits = 2 * message_bits
    distance_cutoff = (
        args.distance_numerator * output_bits + args.distance_denominator - 1
    ) // args.distance_denominator
    block_bits = constituent.block_bits
    local_dimension = constituent.dimension
    minimum_distance = constituent.minimum_distance
    log_nonzero_messages = local_dimension * LOG2 + math.log1p(
        -math.ldexp(1.0, -local_dimension)
    )

    best = np.full(outer_rows + 1, math.inf)
    witness_p = np.full(outer_rows + 1, math.nan)
    witness_u = np.full(outer_rows + 1, math.nan)
    for probability in args.probabilities:
        log_condition = (
            log_choose(block_bits, minimum_distance)
            + minimum_distance * math.log(probability)
            + (block_bits - minimum_distance) * math.log1p(-probability)
        )
        for log_surprisal in args.log_surprisals:
            surprisal = math.exp(log_surprisal)
            z = math.exp(-surprisal)
            zero, active = inner.step_matrices(z, args.memory_bits)
            candidate = (1.0 - probability) * zero + probability * active
            regions = uniform_coefficients(
                inner.log_entries(zero),
                inner.log_entries(candidate),
                outer_rows,
                outer_rows,
            )
            selected = regions[args.minimum_occupation :]
            powered = inner.log_power(selected, block_bits)
            moments = np.logaddexp(powered[..., 0, 0], powered[..., 0, 1])
            for offset, moment in enumerate(moments):
                occupation = args.minimum_occupation + offset
                outer = (
                    log_choose(outer_rows, occupation)
                    + occupation * log_nonzero_messages
                    - occupation * log_condition
                )
                value = outer + min(0.0, float(moment) + distance_cutoff * surprisal)
                if value < best[occupation]:
                    best[occupation] = value
                    witness_p[occupation] = probability
                    witness_u[occupation] = log_surprisal

    rows = [
        {
            "occupation": occupation,
            "log2_upper_diagnostic": float(best[occupation] / LOG2),
            "margin_bits_diagnostic": float(-best[occupation] / LOG2),
            "bernoulli_probability": float(witness_p[occupation]),
            "log_surprisal": float(witness_u[occupation]),
        }
        for occupation in range(args.minimum_occupation, outer_rows + 1)
    ]
    aggregate_log2 = float(
        logsumexp(best[args.minimum_occupation :]) / LOG2
    )
    payload = {
        "schema": "fixed-outer-minimum-distance-randomstepconv-high-q-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "constituent": constituent.name,
            "message_exponent": args.message_exponent,
            "message_bits": message_bits,
            "output_bits": output_bits,
            "outer_rows": outer_rows,
            "memory_bits": args.memory_bits,
            "minimum_distance": minimum_distance,
            "distance_cutoff": distance_cutoff,
            "minimum_occupation": args.minimum_occupation,
        },
        "claim": {
            "covered_occupations": [args.minimum_occupation, outer_rows],
            "aggregate_log2_upper_diagnostic": aggregate_log2,
            "aggregate_margin_bits_diagnostic": -aggregate_log2,
            "closes_40_bits_diagnostic": aggregate_log2 < -40,
            "rows": rows,
        },
        "proof_reduction": [
            "For each nonzero outer word, fix any minimum_distance active coordinates and delete the remaining active coordinates.",
            "RandomStepConv monotonicity makes deletion an upper bound on the tilted low-output moment.",
            "The uniform row permutation maps the retained coordinates to a uniform minimum_distance-subset.",
            "Represent that subset as Bernoulli-p coordinates conditioned on exact weight minimum_distance and divide by the exact conditioning probability.",
        ],
        "limitations": [
            "Nearest binary64 log arithmetic is not outward rounded.",
            "The finite probability and tilt grids are not asserted optimal.",
            "The inner is RandomStepConv, not RM2Sub.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
