#!/usr/bin/env python3
"""Sample a three-band all-Q bound for the high-probability random spectrum."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_single_random_constituent_highprob_renyi import (
    B,
    D,
    K,
    L,
    LOG2,
    MEMORY,
    N,
    log_choose,
    spectrum_caps,
)


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    WORKSTREAM / "single_random_constituent_B256_highprob_band_probe.json"
)
BAND_LIMITS = ((15, 51), (52, 204), (205, 241))


def log_multinomial(counts: tuple[int, ...]) -> float:
    return math.lgamma(sum(counts) + 1) - sum(
        math.lgamma(count + 1) for count in counts
    )


def softmax(values: np.ndarray) -> np.ndarray:
    shifted = values - np.max(values)
    result = np.exp(shifted)
    return result / np.sum(result)


def band_majorant(
    caps: np.ndarray,
    lower: int,
    upper: int,
    *,
    block_bits: int = B,
) -> dict[str, float | int]:
    def objective(probability: float) -> float:
        return max(
            math.log(int(caps[weight]))
            - log_choose(block_bits, weight)
            - weight * math.log(probability)
            - (block_bits - weight) * math.log1p(-probability)
            for weight in range(lower, upper + 1)
            if int(caps[weight])
        )

    result = minimize_scalar(
        objective,
        bounds=(1e-8, 1.0 - 1e-8),
        method="bounded",
        options={"xatol": 1e-14},
    )
    probability = float(result.x)
    values = {
        weight: (
            math.log(int(caps[weight]))
            - log_choose(block_bits, weight)
            - weight * math.log(probability)
            - (block_bits - weight) * math.log1p(-probability)
        )
        for weight in range(lower, upper + 1)
        if int(caps[weight])
    }
    maximizing_weight = max(values, key=values.get)
    return {
        "lower_weight": lower,
        "upper_weight": upper,
        "value_probability": probability,
        "log_majorant": values[maximizing_weight],
        "log2_majorant": values[maximizing_weight] / LOG2,
        "maximizing_weight": maximizing_weight,
    }


def inner_log_moment(bit_probability: float, surprisal: float) -> float:
    z = math.exp(-surprisal)
    zero, active = transfer.step_matrices(z, MEMORY)
    mixed = (1.0 - bit_probability) * zero + bit_probability * active
    powered = transfer.log_power(transfer.log_entries(mixed), N)
    return float(np.logaddexp(powered[0, 0], powered[0, 1]))


def optimize_composition(
    bands: list[dict[str, float | int]],
    active_counts: tuple[int, int, int],
) -> dict[str, object]:
    occupation = sum(active_counts)
    if not 1 <= occupation <= L:
        raise ValueError("composition occupation must lie in [1,L]")
    counts = (L - occupation, *active_counts)
    values = np.asarray(
        [0.0, *(float(band["value_probability"]) for band in bands)]
    )
    log_majorants = np.asarray(
        [0.0, *(float(band["log_majorant"]) for band in bands)]
    )
    active_types = np.flatnonzero(np.asarray(counts) > 0)
    active_type_counts = np.asarray([counts[index] for index in active_types])
    frequencies = active_type_counts / float(L)
    log_type_count = log_multinomial(counts)
    outer_log = float(np.dot(np.asarray(counts), log_majorants))

    def evaluate(point: np.ndarray) -> tuple[float, dict[str, object]]:
        probabilities = softmax(point[:-1])
        surprisal = math.exp(min(4.0, max(-16.0, float(point[-1]))))
        bit_probability = float(
            np.dot(probabilities, values[active_types])
        )
        if not 0.0 < bit_probability < 1.0:
            return math.inf, {}
        log_conditioning = log_type_count + float(
            np.dot(active_type_counts, np.log(probabilities))
        )
        log_moment = inner_log_moment(bit_probability, surprisal)
        raw = log_moment + D * surprisal - B * log_conditioning
        result = log_type_count + outer_log + min(0.0, raw)
        full_probabilities = np.zeros(4)
        full_probabilities[active_types] = probabilities
        return result, {
            "reference_type_probabilities": full_probabilities.tolist(),
            "reference_bit_probability": bit_probability,
            "surprisal": surprisal,
            "log_surprisal": math.log(surprisal),
            "log_conditioning_probability": log_conditioning,
            "log_inner_moment": log_moment,
            "raw_reference_bad_log_upper": raw,
        }

    starts = []
    bit_probability = float(np.dot(frequencies, values[active_types]))
    for scale in (0.7, 1.0, 1.4):
        surprisal = min(3.0, max(1e-6, scale * 2.2 * bit_probability))
        starts.append(np.concatenate((np.log(frequencies), [math.log(surprisal)])))
    best = None
    for start in starts:
        result = minimize(
            lambda point: evaluate(point)[0],
            start,
            method="Nelder-Mead",
            options={"maxiter": 1800, "xatol": 1e-10, "fatol": 1e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    value, details = evaluate(best.x)
    return {
        "occupation": occupation,
        "active_counts": list(active_counts),
        "inactive_count": L - occupation,
        "log2_upper": value / LOG2,
        "margin_bits": -value / LOG2,
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
        **details,
    }


def qL_compositions() -> list[tuple[int, int, int]]:
    result = {(L, 0, 0), (0, L, 0), (0, 0, L)}
    tail_counts = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)
    for tail in tail_counts:
        if tail <= L:
            result.add((tail, L - tail, 0))
            result.add((0, L - tail, tail))
            result.add((tail // 2, L - tail, tail - tail // 2))
    return sorted(result)


def compositions_for_occupation(occupation: int) -> list[tuple[int, int, int]]:
    if occupation <= 32:
        return [
            (low, central, occupation - low - central)
            for low in range(occupation + 1)
            for central in range(occupation - low + 1)
        ]
    result = {
        (occupation, 0, 0),
        (0, occupation, 0),
        (0, 0, occupation),
    }
    probes = sorted(
        {
            1,
            2,
            4,
            8,
            16,
            32,
            64,
            128,
            256,
            512,
            1024,
            occupation // 4,
            occupation // 2,
            3 * occupation // 4,
        }
    )
    for tail in probes:
        if 0 < tail < occupation:
            result.add((tail, occupation - tail, 0))
            result.add((0, occupation - tail, tail))
            result.add(
                (tail // 2, occupation - tail, tail - tail // 2)
            )
    return sorted(result)


def main() -> None:
    global B, D, K, L, MEMORY, N
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-bits", type=int, default=B)
    parser.add_argument("--dimension", type=int, default=K)
    parser.add_argument("--output-bits", type=int, default=N)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--memory-bits", type=int, default=MEMORY)
    parser.add_argument(
        "--band-limits",
        type=int,
        nargs=6,
        metavar=("LO1", "HI1", "LO2", "HI2", "LO3", "HI3"),
    )
    parser.add_argument("--per-shell-failure-exponent", type=int, default=51)
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output_bits % args.block_bits:
        parser.error("block length must divide output length")
    B = args.block_bits
    K = args.dimension
    N = args.output_bits
    L = N // B
    D = (
        args.distance_numerator * N + args.distance_denominator - 1
    ) // args.distance_denominator
    MEMORY = args.memory_bits
    if args.occupations is None:
        args.occupations = [L]
    if args.band_limits is None:
        if B != 256:
            parser.error("nondefault block length requires --band-limits")
        band_limits = BAND_LIMITS
    else:
        band_limits = tuple(
            (args.band_limits[index], args.band_limits[index + 1])
            for index in range(0, 6, 2)
        )
    caps, event = spectrum_caps(
        math.ldexp(1.0, -args.per_shell_failure_exponent),
        block_bits=B,
        dimension=K,
    )
    bands = [band_majorant(caps, *limits, block_bits=B) for limits in band_limits]
    if any(not 1 <= occupation <= L for occupation in args.occupations):
        parser.error(f"occupations must lie in [1,{L}]")
    rows = []
    compositions = sorted(
        {
            composition
            for occupation in args.occupations
            for composition in (
                qL_compositions()
                if occupation == L
                else compositions_for_occupation(occupation)
            )
        }
    )
    for index, composition in enumerate(compositions, start=1):
        row = optimize_composition(bands, composition)
        rows.append(row)
        print(
            f"composition={index}/{len(compositions)},counts={composition},"
            f"margin={row['margin_bits']:.9f}",
            flush=True,
        )
    payload = {
        "schema": "single-random-constituent-highprob-three-band-v1",
        "status": "BINARY64_SAMPLED_QL_COMPOSITIONS",
        "parameters": {
            "outer_code": f"one uniform binary [{B},{K}] generator",
            "outer_rows": L,
            "output_bits": N,
            "distance_cutoff": D,
            "memory_bits": MEMORY,
            "occupations": args.occupations,
        },
        "spectrum_event": event,
        "bands": bands,
        "worst_sampled": min(rows, key=lambda row: float(row["margin_bits"])),
        "rows": rows,
        "limitations": [
            "Only the listed occupations and compositions are evaluated.",
            "The optimizer and arithmetic use nearest binary64.",
            "A complete composition cover and outward verifier remain open.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"bands": bands, "worst": payload["worst_sampled"]}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
