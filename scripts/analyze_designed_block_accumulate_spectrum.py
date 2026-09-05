#!/usr/bin/env python3
"""Spectrum of a direct sum of fixed block codes followed by accumulators."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from flint import fmpz_poly
from scipy.special import logsumexp

from analyze_systematic_expander_accumulate_spectrum import (
    LN2,
    apply_accumulator,
    log2_comb_scalar,
    summary,
)


def log2_positive_integer(value: int) -> float:
    if value <= 0:
        raise ValueError("value must be positive")
    high_bit = value.bit_length() - 1
    shift = max(0, high_bit - 52)
    return math.log2(value >> shift) + shift


def read_spectrum(path: Path, local_length: int) -> list[int]:
    coefficients = [0] * (local_length + 1)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            weight = int(row["weight"])
            count = int(row["count"])
            if weight < 0 or weight > local_length or count < 0:
                raise ValueError(f"invalid spectrum row: {row}")
            coefficients[weight] += count
    return coefficients


def analyze(
    spectrum_path: Path,
    code_name: str,
    local_length: int,
    local_dimension: int,
    blocks: int,
    accumulators: int,
) -> dict[str, object]:
    if local_length <= 0 or local_dimension <= 0 or blocks <= 0:
        raise ValueError("code parameters and block count must be positive")
    if accumulators < 0:
        raise ValueError("accumulators must be nonnegative")

    local = read_spectrum(spectrum_path, local_length)
    if local[0] != 1:
        raise ValueError("a linear injective constituent must have A_0=1")
    if sum(local) != 1 << local_dimension:
        raise ValueError("constituent spectrum mass does not equal 2^dimension")

    begin = perf_counter()
    global_polynomial = fmpz_poly(local) ** blocks
    block_size = local_length * blocks
    dimension = local_dimension * blocks
    spectrum = np.full(block_size + 1, -np.inf, dtype=np.float64)
    for weight in range(block_size + 1):
        coefficient = int(global_polynomial[weight])
        if coefficient:
            spectrum[weight] = log2_positive_integer(coefficient)
    stage_summaries = [{"stage": "direct_sum_constituent", **summary(spectrum)}]
    stage_seconds = [perf_counter() - begin]

    for index in range(1, accumulators + 1):
        begin = perf_counter()
        spectrum = apply_accumulator(spectrum)
        stage_seconds.append(perf_counter() - begin)
        stage_summaries.append({"stage": f"accumulator_{index}", **summary(spectrum)})

    total_log_mass = float(logsumexp(spectrum * LN2) / LN2)
    if abs(total_log_mass - dimension) > 1e-8:
        raise AssertionError(f"spectrum mass changed: log2 mass={total_log_mass}")
    if abs(float(spectrum[0])) > 1e-12:
        raise AssertionError("zero codeword multiplicity changed")

    random_nonzero_log_probability = (
        math.log2((1 << dimension) - 1)
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
        "schema": "riffle-designed-block-accumulate-expected-spectrum-v1",
        "construction": {
            "name": f"Outer DBA-3({code_name})",
            "constituent": code_name,
            "local_length": local_length,
            "local_dimension": local_dimension,
            "blocks": blocks,
            "block_size": block_size,
            "dimension": dimension,
            "accumulators": accumulators,
            "constituent_randomness": "none; one fixed designed code is reused",
            "interleavers": "independent uniform permutations before accumulator stages",
        },
        "method": {
            "constituent_spectrum": str(spectrum_path),
            "direct_sum": "exact integer polynomial power",
            "accumulator_composition": "binary64 log domain",
            "mass_check": "PASS",
        },
        "timing_seconds": stage_seconds,
        "summary": summary(spectrum),
        "stage_summaries": stage_summaries,
        "spectrum": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spectrum", type=Path, required=True)
    parser.add_argument("--code-name", required=True)
    parser.add_argument("--local-length", type=int, required=True)
    parser.add_argument("--local-dimension", type=int, required=True)
    parser.add_argument("--blocks", type=int, required=True)
    parser.add_argument("--accumulators", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(
        args.spectrum,
        args.code_name,
        args.local_length,
        args.local_dimension,
        args.blocks,
        args.accumulators,
    )
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
