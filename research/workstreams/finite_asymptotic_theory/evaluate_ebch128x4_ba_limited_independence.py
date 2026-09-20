#!/usr/bin/env python3
"""Test a primal/dual-distance wrapper for one repeated EBCH128x4--BA code.

The wrapper uses two high-probability events over the BA interleavers.  The
primal event supplies a minimum distance.  The dual event makes a uniform
codeword r-wise independent.  An exact even central moment then bounds the
number of codewords below a chosen weight threshold.  Row-local coordinate
permutations and monotonicity reduce the two resulting categories to uniform
minimum-weight shells.

The transfer calculations are nearest-binary64 diagnostics.  They determine
whether this proof interface is strong enough before an outward verifier is
built.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np

from evaluate_single_random_constituent_shared_two_band_dense import (
    LOG2,
    optimize_counts,
)


WORKSTREAM = Path(__file__).resolve().parent
SPECTRA = WORKSTREAM / "ebch128x4_ba0_16_B512_expected_spectra.json"
OUTPUT = WORKSTREAM / "ebch128x4_ba_limited_independence_probe.json"
B = 512
K = 256
L = 4096
N = 1 << 21
D = (109 * N + 999) // 1000


def logadd2(left: float, right: float) -> float:
    if left == -math.inf:
        return right
    if right == -math.inf:
        return left
    maximum = max(left, right)
    return maximum + math.log2(2.0 ** (left - maximum) + 2.0 ** (right - maximum))


def spectrum_logs(stage: dict[str, object]) -> list[float]:
    result = [-math.inf] * (B + 1)
    for row in stage["spectrum"]:
        value = row["log2_expected_multiplicity"]
        if value is not None:
            result[int(row["weight"])] = float(value)
    return result


def event_distance(log_spectrum: list[float], failure_bits: float) -> tuple[int, float]:
    cumulative = -math.inf
    accepted_cumulative = -math.inf
    cutoff = 0
    for weight in range(1, B + 1):
        cumulative = logadd2(cumulative, log_spectrum[weight])
        if cumulative <= -failure_bits:
            cutoff = weight
            accepted_cumulative = cumulative
    return cutoff + 1, accepted_cumulative


def even_binomial_central_moment(order: int) -> Fraction:
    if order <= 0 or order % 2:
        raise ValueError("central-moment order must be positive and even")
    # X-256 is integral for X~Bin(512,1/2), so this is exact.
    numerator = sum(
        math.comb(B, weight) * (weight - B // 2) ** order
        for weight in range(B + 1)
    )
    return Fraction(numerator, 1 << B)


def tail_log2_upper(order: int, cutoff: int) -> float:
    if not 0 <= cutoff < B // 2:
        raise ValueError("tail cutoff must be below B/2")
    moment = even_binomial_central_moment(order)
    return math.log2(moment.numerator) - math.log2(moment.denominator) - order * math.log2(B // 2 - cutoff)


def shell_majorant_log(
    count_log2: float,
    retained_weight: int,
) -> tuple[float, float]:
    probability = retained_weight / B
    shell_log_probability = (
        math.lgamma(B + 1)
        - math.lgamma(retained_weight + 1)
        - math.lgamma(B - retained_weight + 1)
        + retained_weight * math.log(probability)
        + (B - retained_weight) * math.log1p(-probability)
    )
    return probability, count_log2 * LOG2 - shell_log_probability


def sampled_defect_counts(occupation: int) -> list[int]:
    candidates = {
        0,
        1,
        min(2, occupation),
        occupation // 8,
        occupation // 4,
        occupation // 2,
        occupation,
    }
    return sorted(value for value in candidates if 0 <= value <= occupation)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stages", type=int, nargs="+", default=[10, 11, 12, 13])
    parser.add_argument("--event-failure-bits", type=float, default=44.0)
    parser.add_argument("--cutoffs", type=int, nargs="+", default=[128, 144, 160, 176, 192])
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
        default=[3, 8, 16, 32, 64, 128, 160, 256, 512, 1024, 2048, 4096],
    )
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    payload = json.loads(SPECTRA.read_text(encoding="utf-8"))
    primal = {int(row["accumulators"]): spectrum_logs(row) for row in payload["stages"]}
    dual = {int(row["accumulators"]): spectrum_logs(row) for row in payload["dual_stages"]}
    rows = []
    for stage in args.stages:
        primal_distance, primal_failure = event_distance(primal[stage], args.event_failure_bits)
        dual_distance, dual_failure = event_distance(dual[stage], args.event_failure_bits)
        order = dual_distance - 1
        if order % 2:
            order -= 1
        if order < 2:
            continue
        for cutoff in args.cutoffs:
            tail_log2 = tail_log2_upper(order, cutoff)
            low_count_log2 = min(float(K), K + tail_log2)
            low_probability, low_majorant = shell_majorant_log(
                low_count_log2, primal_distance
            )
            regular_probability, regular_majorant = shell_majorant_log(
                float(K), cutoff + 1
            )
            values = np.asarray([0.0, low_probability, regular_probability])
            majorants = np.asarray([0.0, low_majorant, regular_majorant])
            probes = []
            for occupation in args.occupations:
                for defects in sampled_defect_counts(occupation):
                    counts = (L - occupation, defects, occupation - defects)
                    result = optimize_counts(
                        counts,
                        values=values,
                        log_majorants=majorants,
                        block_bits=B,
                        distance_cutoff=D,
                        memory_bits=args.memory_bits,
                        output_bits=N,
                        temperatures=(0.6, 1.0),
                        scales=(1.4, 3.0),
                        maximum_iterations=600,
                    )
                    probes.append(result)
            worst = min(probes, key=lambda row: float(row["margin_bits"]))
            row = {
                "accumulator_stages": stage,
                "primal_minimum_distance": primal_distance,
                "dual_minimum_distance": dual_distance,
                "independence_order": order,
                "primal_event_failure_log2_upper": primal_failure,
                "dual_event_failure_log2_upper": dual_failure,
                "tail_cutoff": cutoff,
                "tail_probability_log2_upper": tail_log2,
                "low_count_log2_upper": low_count_log2,
                "low_retained_weight": primal_distance,
                "regular_retained_weight": cutoff + 1,
                "low_majorant_log2": low_majorant / LOG2,
                "regular_majorant_log2": regular_majorant / LOG2,
                "worst_sampled": worst,
                "probes": probes,
            }
            rows.append(row)
            print(
                f"stage,{stage},r,{order},cutoff,{cutoff},tail,{tail_log2:.6f},"
                f"worst_margin,{worst['margin_bits']:.6f},"
                f"q,{worst['occupation']},defects,{worst['defect_rows']}",
                flush=True,
            )

    result = {
        "schema": "ebch128x4-ba-limited-independence-v1",
        "status": "BINARY64_SAMPLED_TRANSFER_DIAGNOSTIC",
        "parameters": {
            "outer": "four EBCH [128,64,22] blocks followed by uniform-interleaved terminated accumulators",
            "outer_block_bits": B,
            "outer_dimension": K,
            "outer_rows": L,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": args.memory_bits,
            "event_failure_bits_per_primal_or_dual_event": args.event_failure_bits,
            "sampled_occupations": args.occupations,
        },
        "rows": rows,
        "proof_interface": {
            "tail_bound": "dual distance gives r-wise independent coordinates; exact even Binomial(512,1/2) central moment plus Markov bounds the low-weight fraction",
            "category_reduction": "primal distance retains d bits in low rows; rows above the cutoff retain cutoff+1 bits; row-local permutations make both retained supports uniform",
            "transfer": "two-category Bernoulli change of measure and RandomStepConv-M22",
        },
        "limitations": [
            "The occupation and defect-composition sets are sampled, not complete covers.",
            "The transfer and expected spectra use nearest binary64 arithmetic.",
            "A passing diagnostic would still require outward primal/dual event checks and complete sparse/dense covers.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
