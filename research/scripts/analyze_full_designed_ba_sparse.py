#!/usr/bin/env python3
"""Finite-N sparse-spectrum evaluator for full-length designed BA-3.

The code is the direct sum of a fixed injective binary constituent C, followed
by two independent uniform interleavers and two length-N accumulators.  This
script computes the exact low-weight spectrum of the direct sum and an upper
bound on its contribution to final weight at most D.

For the second accumulator, the cumulative probability is bounded by its last
term times a geometric series.  The bound is valid below the mode because the
ratio p(t-1)/p(t) increases with t.  Inputs of weights one and two are handled
in closed form.  At small N, --cdf exact evaluates the complete cumulative sum
and provides a validation oracle.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.special import gammaln, logsumexp


LN2 = math.log(2.0)


def log2_positive_integer(value: int) -> float:
    if value <= 0:
        raise ValueError("value must be positive")
    high_bit = value.bit_length() - 1
    shift = max(0, high_bit - 52)
    return math.log2(value >> shift) + shift


def log2_comb_scalar(n: int, k: int) -> float:
    if n < 0 or k < 0 or k > n:
        return -math.inf
    return float((gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)) / LN2)


def log2_comb_vector(n: np.ndarray, k: int) -> np.ndarray:
    result = np.full(n.shape, -np.inf, dtype=np.float64)
    valid = (n >= 0) & (k >= 0) & (k <= n)
    nv = n[valid]
    result[valid] = (
        gammaln(nv + 1) - gammaln(k + 1) - gammaln(nv - k + 1)
    ) / LN2
    return result


def read_spectrum(path: Path, local_length: int, local_dimension: int) -> list[int]:
    spectrum = [0] * (local_length + 1)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            weight = int(row["weight"])
            count = int(row["count"])
            if weight < 0 or weight > local_length or count < 0:
                raise ValueError(f"invalid spectrum row: {row}")
            spectrum[weight] += count
    if spectrum[0] != 1 or sum(spectrum) != 1 << local_dimension:
        raise ValueError("spectrum is not an injective linear code of the declared dimension")
    return spectrum


def convolve_truncated(left: list[int], right: list[int], cap: int) -> list[int]:
    result = [0] * (cap + 1)
    left_support = [(i, value) for i, value in enumerate(left) if value]
    right_support = [(i, value) for i, value in enumerate(right) if value]
    for i, a in left_support:
        for j, b in right_support:
            if i + j > cap:
                break
            result[i + j] += a * b
    return result


def direct_sum_low_spectrum(
    local: list[int], copies: int, cap: int
) -> list[int]:
    nonzero = local.copy()
    nonzero[0] = 0
    minimum_weight = next(weight for weight, count in enumerate(nonzero) if count)
    maximum_active_blocks = min(copies, cap // minimum_weight)
    result = [0] * (cap + 1)
    result[0] = 1
    active_power = [0] * (cap + 1)
    active_power[0] = 1
    for active_blocks in range(1, maximum_active_blocks + 1):
        active_power = convolve_truncated(active_power, nonzero, cap)
        placement_count = math.comb(copies, active_blocks)
        for weight, count in enumerate(active_power):
            if count:
                result[weight] += placement_count * count
    return result


def accumulator_logpmf(input_weight: int, outputs: np.ndarray, length: int) -> np.ndarray:
    if input_weight <= 0:
        raise ValueError("positive input weight required")
    half_down = input_weight // 2
    half_up = (input_weight + 1) // 2
    return (
        log2_comb_vector(length - outputs, half_down)
        + log2_comb_vector(outputs - 1, half_up - 1)
        - log2_comb_scalar(length, input_weight)
    )


def second_accumulator_cdf_upper(length: int, distance: int) -> np.ndarray:
    max_input = min(length, 2 * distance)
    result = np.full(max_input + 1, -np.inf, dtype=np.float64)
    result[0] = 0.0
    if max_input >= 1:
        result[1] = math.log2(distance / length)
    if max_input >= 2:
        numerator = distance * length - distance * (distance + 1) / 2
        result[2] = math.log2(numerator / math.comb(length, 2))

    weights = np.arange(3, max_input + 1, dtype=np.int64)
    half_down = weights // 2
    half_up = (weights + 1) // 2
    feasible = half_up <= distance
    log_pmf = np.full(weights.shape, -np.inf, dtype=np.float64)
    if np.any(feasible):
        wf = weights[feasible]
        ad = half_down[feasible]
        au = half_up[feasible]
        log_pmf[feasible] = (
            np.array([log2_comb_scalar(length - distance, int(a)) for a in ad])
            + np.array([log2_comb_scalar(distance - 1, int(b - 1)) for b in au])
            - np.array([log2_comb_scalar(length, int(w)) for w in wf])
        )
        ratio = (
            (length - distance + 1) * (distance - au)
            / ((length - distance + 1 - ad) * (distance - 1))
        )
        valid_geometric = ratio < 1.0
        correction = np.zeros(ratio.shape, dtype=np.float64)
        correction[valid_geometric] = -np.log2(1.0 - ratio[valid_geometric])
        correction[~valid_geometric] = -log_pmf[feasible][~valid_geometric]
        log_pmf[feasible] = np.minimum(0.0, log_pmf[feasible] + correction)
    result[3:] = log_pmf
    return result


def second_accumulator_cdf_exact(length: int, distance: int) -> np.ndarray:
    max_input = min(length, 2 * distance)
    result = np.full(max_input + 1, -np.inf, dtype=np.float64)
    result[0] = 0.0
    for input_weight in range(1, max_input + 1):
        low = (input_weight + 1) // 2
        high = min(distance, length - input_weight // 2)
        if low <= high:
            outputs = np.arange(low, high + 1, dtype=np.int64)
            logs = accumulator_logpmf(input_weight, outputs, length)
            result[input_weight] = float(logsumexp(logs * LN2) / LN2)
    return result


def sparse_failure_contributions(
    initial: list[int], length: int, distance: int, second_cdf: np.ndarray
) -> tuple[list[dict[str, float | int]], float]:
    max_middle = len(second_cdf) - 1
    middle = np.arange(1, max_middle + 1, dtype=np.int64)
    rows: list[dict[str, float | int]] = []
    total = -math.inf
    for outer_weight, multiplicity in enumerate(initial):
        if outer_weight == 0 or multiplicity == 0:
            continue
        low = (outer_weight + 1) // 2
        high = min(max_middle, length - outer_weight // 2)
        if low > high:
            continue
        selected = middle[low - 1 : high]
        first_logpmf = accumulator_logpmf(outer_weight, selected, length)
        conditional = float(
            logsumexp((first_logpmf + second_cdf[selected]) * LN2) / LN2
        )
        log_multiplicity = log2_positive_integer(multiplicity)
        contribution = log_multiplicity + conditional
        total = float(np.logaddexp2(total, contribution))
        rows.append(
            {
                "outer_weight": outer_weight,
                "log2_multiplicity": log_multiplicity,
                "log2_conditional_bad_probability": conditional,
                "log2_expected_bad_contribution": contribution,
            }
        )
    return rows, total


def analyze(
    spectrum_path: Path,
    code_name: str,
    local_length: int,
    local_dimension: int,
    length: int,
    distance: int,
    outer_weight_cap: int,
    cdf_mode: str,
) -> dict[str, object]:
    if length <= 0 or length % local_length:
        raise ValueError("full length must be a positive multiple of local length")
    if 2 * local_dimension != local_length:
        raise ValueError("this evaluator expects a rate-1/2 constituent")
    if distance <= 0 or distance >= length // 2:
        raise ValueError("distance must lie in (0,N/2)")

    local = read_spectrum(spectrum_path, local_length, local_dimension)
    copies = length // local_length
    begin = perf_counter()
    initial = direct_sum_low_spectrum(local, copies, outer_weight_cap)
    initial_seconds = perf_counter() - begin

    begin = perf_counter()
    if cdf_mode == "exact":
        second_cdf = second_accumulator_cdf_exact(length, distance)
        cdf_guarantee = "exact"
    elif cdf_mode == "upper":
        second_cdf = second_accumulator_cdf_upper(length, distance)
        cdf_guarantee = "upper bound by monotone geometric tail"
    else:
        raise ValueError(f"unknown cdf mode: {cdf_mode}")
    cdf_seconds = perf_counter() - begin

    begin = perf_counter()
    rows, total = sparse_failure_contributions(
        initial, length, distance, second_cdf
    )
    contribution_seconds = perf_counter() - begin
    dominant = max(rows, key=lambda row: row["log2_expected_bad_contribution"])

    return {
        "schema": "full-designed-ba3-sparse-finite-v1",
        "construction": {
            "name": f"DBA-3({code_name})",
            "length": length,
            "dimension": length // 2,
            "rate": 0.5,
            "constituent": code_name,
            "local_length": local_length,
            "local_dimension": local_dimension,
            "copies": copies,
            "global_interleavers": 2,
            "global_accumulators": 2,
        },
        "target": {
            "distance": distance,
            "relative_distance": distance / length,
            "outer_weight_cap": outer_weight_cap,
        },
        "method": {
            "direct_sum_low_spectrum": "exact integer active-block expansion",
            "second_accumulator_cdf": cdf_guarantee,
            "scope": "all messages whose constituent-output weight is at most outer_weight_cap",
            "limitation": "a separate bridge is required for larger constituent-output weights",
        },
        "result": {
            "sparse_log2_expected_bad_mass": total,
            "sparse_margin_bits": -total,
            "dominant_outer_weight": dominant["outer_weight"],
            "dominant_log2_contribution": dominant["log2_expected_bad_contribution"],
        },
        "timing_seconds": {
            "initial_spectrum": initial_seconds,
            "second_cdf": cdf_seconds,
            "contributions": contribution_seconds,
        },
        "contributions": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spectrum", type=Path, required=True)
    parser.add_argument("--code-name", required=True)
    parser.add_argument("--local-length", type=int, required=True)
    parser.add_argument("--local-dimension", type=int, required=True)
    parser.add_argument("--length", type=int, default=1 << 21)
    parser.add_argument("--distance", type=int)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument("--outer-weight-cap", type=int, default=256)
    parser.add_argument("--cdf", choices=("exact", "upper"), default="upper")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    distance = args.distance
    if distance is None:
        distance = math.floor(args.relative_distance * args.length)
    result = analyze(
        args.spectrum,
        args.code_name,
        args.local_length,
        args.local_dimension,
        args.length,
        distance,
        args.outer_weight_cap,
        args.cdf,
    )
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
