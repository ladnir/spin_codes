#!/usr/bin/env python3
"""Categorical dense diagnostic after merging the two spectrum tails."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_single_random_constituent_highprob_bands import (
    band_majorant,
    log_multinomial,
    softmax,
)
from evaluate_single_random_constituent_highprob_renyi import LOG2, spectrum_caps


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_OUTPUT = WORKSTREAM / "single_random_constituent_B512_shared_two_band_dense_probe_s22.json"


def inner_log_moment(
    bit_probability: float,
    surprisal: float,
    *,
    memory_bits: int,
    output_bits: int,
) -> float:
    z = math.exp(-surprisal)
    zero, active = transfer.step_matrices(z, memory_bits)
    mixed = (1.0 - bit_probability) * zero + bit_probability * active
    powered = transfer.log_power(transfer.log_entries(mixed), output_bits)
    return float(np.logaddexp(powered[0, 0], powered[0, 1]))


def evaluate_with_witness(
    counts: tuple[float, float, float],
    probabilities: np.ndarray,
    surprisal: float,
    *,
    values: np.ndarray,
    log_majorants: np.ndarray,
    block_bits: int,
    distance_cutoff: int,
    memory_bits: int,
    output_bits: int,
) -> tuple[float, float, dict[str, object]]:
    active = np.asarray(counts) > 0.0
    if np.any(probabilities[active] <= 0.0):
        return math.inf, math.inf, {}
    log_type_count = (
        math.lgamma(sum(counts) + 1.0)
        - sum(math.lgamma(count + 1.0) for count in counts)
    )
    outer = log_type_count + float(np.dot(np.asarray(counts), log_majorants))
    bit_probability = float(np.dot(probabilities, values))
    log_conditioning = log_type_count + float(
        np.dot(np.asarray(counts)[active], np.log(probabilities[active]))
    )
    moment = inner_log_moment(
        bit_probability,
        surprisal,
        memory_bits=memory_bits,
        output_bits=output_bits,
    )
    raw = moment + distance_cutoff * surprisal - block_bits * log_conditioning
    return outer + min(0.0, raw), raw, {
        "reference_type_probabilities": probabilities.tolist(),
        "reference_bit_probability": bit_probability,
        "surprisal": surprisal,
        "log_surprisal": math.log(surprisal),
        "log_conditioning_probability": log_conditioning,
        "log_inner_moment": moment,
        "raw_reference_bad_log_upper": raw,
    }


def optimize_counts(
    counts: tuple[int, int, int],
    *,
    values: np.ndarray,
    log_majorants: np.ndarray,
    block_bits: int,
    distance_cutoff: int,
    memory_bits: int,
    output_bits: int,
    temperatures: tuple[float, ...] = (0.5, 0.75, 1.0, 1.25),
    scales: tuple[float, ...] = (0.7, 1.4, 3.0, 6.0),
    maximum_iterations: int = 1600,
) -> dict[str, object]:
    active = np.flatnonzero(np.asarray(counts) > 0)
    frequencies = np.asarray([counts[index] for index in active], dtype=float)
    frequencies /= np.sum(frequencies)

    def unpack(point: np.ndarray) -> tuple[np.ndarray, float]:
        probabilities = np.zeros(3)
        probabilities[active] = softmax(point[:-1])
        surprisal = math.exp(min(4.0, max(-16.0, float(point[-1]))))
        return probabilities, surprisal

    def objective(point: np.ndarray) -> float:
        probabilities, surprisal = unpack(point)
        return evaluate_with_witness(
            counts,
            probabilities,
            surprisal,
            values=values,
            log_majorants=log_majorants,
            block_bits=block_bits,
            distance_cutoff=distance_cutoff,
            memory_bits=memory_bits,
            output_bits=output_bits,
        )[0]

    density = (counts[1] + counts[2]) / sum(counts)
    bit_probability = density * float(
        np.dot(frequencies, values[active])
    )
    starts = []
    for temperature in temperatures:
        start_probabilities = softmax(temperature * np.log(frequencies))
        start_bit_probability = float(
            np.dot(start_probabilities, values[active])
        )
        for scale in scales:
            starts.append(
                np.concatenate(
                    (
                        np.log(start_probabilities),
                        [
                            math.log(
                                max(1e-6, scale * 2.2 * start_bit_probability)
                            )
                        ],
                    )
                )
            )
    best = None
    for start in starts:
        result = minimize(
            objective,
            start,
            method="Nelder-Mead",
            options={
                "maxiter": maximum_iterations,
                "xatol": 1e-9,
                "fatol": 1e-6,
            },
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    probabilities, surprisal = unpack(best.x)
    value, raw, details = evaluate_with_witness(
        counts,
        probabilities,
        surprisal,
        values=values,
        log_majorants=log_majorants,
        block_bits=block_bits,
        distance_cutoff=distance_cutoff,
        memory_bits=memory_bits,
        output_bits=output_bits,
    )
    return {
        "counts": list(counts),
        "occupation": counts[1] + counts[2],
        "defect_rows": counts[1],
        "central_rows": counts[2],
        "log2_upper": value / LOG2,
        "margin_bits": -value / LOG2,
        "raw_negative": raw < 0.0,
        "optimizer_success": bool(best.success),
        **details,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--occupations", type=int, nargs="+", default=[129])
    parser.add_argument("--block-bits", type=int, default=512)
    parser.add_argument("--dimension", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=1 << 21)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--per-shell-failure-exponent", type=int, default=51)
    parser.add_argument("--low-upper", type=int, default=79)
    parser.add_argument("--central-upper", type=int, default=432)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output_bits % args.block_bits:
        parser.error("block length must divide output length")
    outer_rows = args.output_bits // args.block_bits
    if any(not 1 <= q <= outer_rows for q in args.occupations):
        parser.error("occupation outside valid range")
    distance_cutoff = (
        args.distance_numerator * args.output_bits
        + args.distance_denominator
        - 1
    ) // args.distance_denominator

    caps, event = spectrum_caps(
        math.ldexp(1.0, -args.per_shell_failure_exponent),
        block_bits=args.block_bits,
        dimension=args.dimension,
    )
    support = [weight for weight in range(1, args.block_bits + 1) if int(caps[weight])]
    low = band_majorant(caps, min(support), args.low_upper, block_bits=args.block_bits)
    central = band_majorant(
        caps, args.low_upper + 1, args.central_upper, block_bits=args.block_bits
    )
    high = band_majorant(
        caps, args.central_upper + 1, max(support), block_bits=args.block_bits
    )
    defect_probability = min(
        float(low["value_probability"]),
        1.0 - float(high["value_probability"]),
    )
    defect_log_majorant = math.log(2.0) + max(
        float(low["log_majorant"]), float(high["log_majorant"])
    )
    values = np.asarray([0.0, defect_probability, float(central["value_probability"])])
    log_majorants = np.asarray(
        [0.0, defect_log_majorant, float(central["log_majorant"])]
    )

    rows = []
    total_rows = sum(q + 1 for q in args.occupations)
    completed = 0
    for occupation in args.occupations:
        for defects in range(occupation + 1):
            counts = (outer_rows - occupation, defects, occupation - defects)
            row = optimize_counts(
                counts,
                values=values,
                log_majorants=log_majorants,
                block_bits=args.block_bits,
                distance_cutoff=distance_cutoff,
                memory_bits=args.memory_bits,
                output_bits=args.output_bits,
            )
            rows.append(row)
            completed += 1
            print(
                f"composition,{completed},{total_rows},q,{occupation},d,{defects},margin,{row['margin_bits']:.9f}",
                flush=True,
            )

    aggregate = float(logsumexp([float(row["log2_upper"]) * LOG2 for row in rows]))
    payload = {
        "schema": "single-random-constituent-shared-two-band-dense-v1",
        "status": "BINARY64_CATEGORICAL_DIAGNOSTIC",
        "parameters": {
            "outer_code": f"one uniform binary [{args.block_bits},{args.dimension}] generator",
            "outer_rows": outer_rows,
            "output_bits": args.output_bits,
            "distance_cutoff": distance_cutoff,
            "memory_bits": args.memory_bits,
            "occupations": args.occupations,
        },
        "spectrum_event": event,
        "bands": {"low": low, "central": central, "high": high},
        "merged_defect": {
            "value_probability": defect_probability,
            "log2_majorant": defect_log_majorant / LOG2,
        },
        "claim": {
            "aggregate_log2_upper": aggregate / LOG2,
            "aggregate_margin_bits": -aggregate / LOG2,
            "worst": max(rows, key=lambda row: float(row["log2_upper"])),
        },
        "rows": rows,
        "limitations": [
            "Only the listed occupation slices are evaluated.",
            "Nearest binary64 arithmetic is not an outward certificate.",
            "A complete dense convex cell cover remains open.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
