#!/usr/bin/env python3
"""Expected spectrum of a systematic expander followed by accumulators.

The rate-1/2 expander maps a k-bit message x to (x, p), where each of the k
parity bits is the XOR of a uniformly sampled size-d subset of x.  Checks are
sampled independently.  Independent uniform interleavers precede accumulator
stages.  All sums are nonnegative, so binary64 log-domain composition is
stable.  Use the exact ParityMix analyzer at small sizes to validate changes.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.special import gammaln, logsumexp


LN2 = math.log(2.0)


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


def odd_intersection_count(message_size: int, degree: int, weight: int) -> int:
    return sum(
        math.comb(weight, intersection)
        * math.comb(message_size - weight, degree - intersection)
        for intersection in range(1, degree + 1, 2)
        if intersection <= weight and degree - intersection <= message_size - weight
    )


def systematic_expander_spectrum(message_size: int, degree: int) -> np.ndarray:
    block_size = 2 * message_size
    denominator = math.comb(message_size, degree)
    spectrum = np.full(block_size + 1, -np.inf, dtype=np.float64)
    parity_weights = np.arange(message_size + 1)
    log_parity_binomials = np.array(
        [log2_comb_scalar(message_size, h) for h in parity_weights]
    )

    for input_weight in range(message_size + 1):
        odd = odd_intersection_count(message_size, degree, input_weight)
        input_log_count = log2_comb_scalar(message_size, input_weight)
        if odd == 0:
            local = np.array([input_log_count])
        elif odd == denominator:
            local = np.full(message_size + 1, -np.inf)
            local[-1] = input_log_count
        else:
            probability = odd / denominator
            local = (
                input_log_count
                + log_parity_binomials
                + parity_weights * math.log2(probability)
                + (message_size - parity_weights) * math.log2(1.0 - probability)
            )
        indices = slice(input_weight, input_weight + len(local))
        spectrum[indices] = np.logaddexp2(spectrum[indices], local)
    return spectrum


def apply_accumulator(log_spectrum: np.ndarray) -> np.ndarray:
    block_size = len(log_spectrum) - 1
    result = np.full(block_size + 1, -np.inf, dtype=np.float64)
    result[0] = log_spectrum[0]
    for input_weight in range(1, block_size + 1):
        if not math.isfinite(float(log_spectrum[input_weight])):
            continue
        half_down = input_weight // 2
        half_up = (input_weight + 1) // 2
        low = half_up
        high = block_size - half_down
        if low > high:
            continue
        output_weights = np.arange(low, high + 1)
        log_counts = (
            log2_comb_vector(block_size - output_weights, half_down)
            + log2_comb_vector(output_weights - 1, half_up - 1)
        )
        contributions = (
            log_spectrum[input_weight]
            - log2_comb_scalar(block_size, input_weight)
            + log_counts
        )
        result[low : high + 1] = np.logaddexp2(
            result[low : high + 1], contributions
        )
    return result


def cumulative_log2(log_spectrum: np.ndarray) -> np.ndarray:
    result = np.full(len(log_spectrum), -np.inf, dtype=np.float64)
    running = -math.inf
    for weight in range(1, len(log_spectrum)):
        running = float(np.logaddexp2(running, log_spectrum[weight]))
        result[weight] = running
    return result


def first_at_least(values: np.ndarray, threshold: float) -> int | None:
    for weight in range(1, len(values)):
        if values[weight] >= threshold:
            return weight
    return None


def summary(log_spectrum: np.ndarray) -> dict[str, object]:
    cumulative = cumulative_log2(log_spectrum)
    block_size = len(log_spectrum) - 1
    return {
        "first_weight_expected_at_least_1": first_at_least(log_spectrum, 0.0),
        "first_weight_cumulative_expected_at_least_1": first_at_least(cumulative, 0.0),
        "cumulative_log2_expected_nonzero_codewords_through_weight": {
            str(cutoff): float(cumulative[cutoff])
            for cutoff in (31, 37, 63, 95, 127)
            if cutoff <= block_size
        },
    }


def analyze(message_size: int, degree: int, accumulators: int) -> dict[str, object]:
    if message_size <= 0:
        raise ValueError("message size must be positive")
    if degree <= 0 or degree > message_size:
        raise ValueError("degree must be in [1,message size]")
    if accumulators < 0:
        raise ValueError("accumulators must be nonnegative")

    begin = perf_counter()
    spectrum = systematic_expander_spectrum(message_size, degree)
    stage_summaries = [{"stage": "systematic_expander", **summary(spectrum)}]
    stage_seconds = [perf_counter() - begin]
    for index in range(1, accumulators + 1):
        begin = perf_counter()
        spectrum = apply_accumulator(spectrum)
        stage_seconds.append(perf_counter() - begin)
        stage_summaries.append({"stage": f"accumulator_{index}", **summary(spectrum)})

    total_log_mass = float(logsumexp(spectrum * LN2) / LN2)
    if abs(total_log_mass - message_size) > 1e-8:
        raise AssertionError(f"spectrum mass changed: log2 mass={total_log_mass}")
    if abs(float(spectrum[0])) > 1e-12:
        raise AssertionError("zero codeword multiplicity changed")

    block_size = 2 * message_size
    random_nonzero_log_probability = (
        math.log2((1 << message_size) - 1)
        - math.log2((1 << block_size) - 1)
    )
    rows = []
    for weight in range(block_size + 1):
        random_value = (
            0.0
            if weight == 0
            else log2_comb_scalar(block_size, weight) + random_nonzero_log_probability
        )
        value = float(spectrum[weight])
        rows.append(
            {
                "weight": weight,
                "log2_expected_multiplicity": None if not math.isfinite(value) else value,
                "random_linear_log2_expected_multiplicity": random_value,
                "excess_over_random_bits": (
                    None if not math.isfinite(value) else value - random_value
                ),
            }
        )
    return {
        "schema": "riffle-systematic-expander-accumulate-expected-spectrum-v1",
        "ensemble": {
            "message_size": message_size,
            "block_size": block_size,
            "rate": 0.5,
            "expander_degree": degree,
            "accumulators": accumulators,
            "randomness": "independent random check supports and independent uniform interleavers",
        },
        "method": {
            "arithmetic": "binary64 log-domain composition of nonnegative exact combinatorial terms",
            "mass_check": "PASS",
        },
        "timing_seconds": stage_seconds,
        "summary": summary(spectrum),
        "stage_summaries": stage_summaries,
        "spectrum": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--message-size", type=int, default=512)
    parser.add_argument("--degree", type=int, default=8)
    parser.add_argument("--accumulators", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.message_size, args.degree, args.accumulators)
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
