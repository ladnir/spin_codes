#!/usr/bin/env python3
"""Evaluate Q=1 for one random full-rank constituent reused across rows.

For each block length B, sample one binary [B,B/2] full-rank linear map
uniformly and reuse it in every outer row.  Sampling may be implemented by
rejection from iid generator matrices.  For a fixed nonzero local message,
its image is uniform over the 2^B-1 nonzero B-bit vectors.  Consequently the
expected nonzero spectrum is exact and suffices for the joint outer/routing/
RandomStepConv Q=1 union bound.

The calculation uses nearest binary64 arithmetic and is diagnostic.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from scipy.special import gammaln


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
sys.path.insert(0, str(WORKSTREAM))

from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    Constituent,
    evaluate_case,
    tilt_grid,
)


LOG2 = math.log(2.0)
DEFAULT_OUTPUT = HERE / "reused_random_constituent_q1_phase.json"


def log_nonzero_count(bits: int) -> float:
    return bits * LOG2 + math.log1p(-math.ldexp(1.0, -bits))


def expected_spectrum(block_bits: int) -> dict[int, float]:
    dimension = block_bits // 2
    log_messages = log_nonzero_count(dimension)
    log_outputs = log_nonzero_count(block_bits)
    result = {}
    for weight in range(1, block_bits + 1):
        log_binomial = (
            gammaln(block_bits + 1)
            - gammaln(weight + 1)
            - gammaln(block_bits - weight + 1)
        )
        result[weight] = math.exp(log_messages + log_binomial - log_outputs)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-bits", type=int, nargs="+", default=[32, 128, 512])
    parser.add_argument(
        "--message-exponents", type=int, nargs="+", default=list(range(8, 17))
    )
    parser.add_argument(
        "--memory-bits", type=int, nargs="+", default=[8, 12, 16, 18, 20, 22, 24]
    )
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--grid-min", type=float, default=-12.0)
    parser.add_argument("--grid-max", type=float, default=1.0)
    parser.add_argument("--grid-step", type=float, default=0.5)
    parser.add_argument("--fine-min", type=float, default=-10.0)
    parser.add_argument("--fine-max", type=float, default=-5.0)
    parser.add_argument("--fine-step", type=float, default=0.1)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    for block_bits in args.block_bits:
        if block_bits <= 0 or block_bits % 2:
            parser.error("every block length must be a positive even integer")
    if not 0 < args.distance_numerator < args.distance_denominator:
        parser.error("distance fraction must lie in (0,1)")

    tilts = tilt_grid(args)
    cases = []
    for block_bits in args.block_bits:
        dimension = block_bits // 2
        constituent = Constituent(
            name=f"random full-rank [{block_bits},{dimension}]",
            block_bits=block_bits,
            dimension=dimension,
            minimum_distance=1,
            spectrum_path=Path(),
        )
        spectrum = expected_spectrum(block_bits)
        for message_exponent in args.message_exponents:
            if (1 << message_exponent) < dimension:
                continue
            for memory_bits in args.memory_bits:
                print(
                    f"case,random-full-rank-{block_bits},kexp,{message_exponent},memory,{memory_bits}",
                    flush=True,
                )
                cases.append(
                    evaluate_case(
                        constituent=constituent,
                        spectrum=spectrum,
                        message_exponent=message_exponent,
                        memory_bits=memory_bits,
                        distance_numerator=args.distance_numerator,
                        distance_denominator=args.distance_denominator,
                        tilts=tilts,
                    )
                )

    payload = {
        "schema": "reused-random-full-rank-constituent-q1-phase-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "probability_space": {
            "outer": (
                "for each case, one uniform full-rank binary [B,B/2] linear map "
                "is sampled once and reused in every outer row"
            ),
            "routing": (
                "independent uniform coordinate permutation per row and independent "
                "uniform position permutation per transposed region"
            ),
            "inner": (
                "independent RandomStepConv map at every position, sampled once and "
                "shared by all messages"
            ),
        },
        "parameters": {
            "block_bits": args.block_bits,
            "message_exponents": args.message_exponents,
            "memory_bits": args.memory_bits,
            "distance_numerator": args.distance_numerator,
            "distance_denominator": args.distance_denominator,
            "tilt_count": len(tilts),
        },
        "cases": cases,
        "proof_reduction": [
            "Sample one full-rank constituent uniformly; rejection sampling removes any rank-failure event.",
            "For each fixed nonzero local message, its constituent image is uniform over all nonzero B-bit vectors.",
            "Average the exact binomial shell counts against the same Q=1 uniform-routing and RandomStepConv transfer used for BCH and RM.",
        ],
        "limitations": [
            "Nearest binary64 log arithmetic is not outward rounded.",
            "Only occupation Q=1 is covered.",
            "The finite tilt grid is not asserted optimal.",
            "This is one random constituent reused across rows, not a fresh constituent per row and not one global random code.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
