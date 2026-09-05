#!/usr/bin/env python3
"""Derive expected spectra for inexpensive parity-breaking outer shears."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_riffle_bitshuffle_splitstate_regular_bulk import (
    LOG2,
    logsumexp,
    modeled_even_floor_spectrum_logs,
    spectrum_bernoulli_envelope_log,
)


DEFAULT_OUTPUT = Path(
    "constructions/"
    "riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/"
    "receipts/parity_break_outer256_d38_expected_spectrum.json"
)


def expected_message_linear_spectrum(source_logs: np.ndarray) -> np.ndarray:
    length = len(source_logs) - 1
    result = np.full(length + 1, -math.inf)
    result[0] = 0.0
    for source_weight in range(1, length + 1):
        source_log = float(source_logs[source_weight])
        if not math.isfinite(source_log):
            continue
        contributions = (
            (source_weight, 0.5),
            (source_weight - 1, 0.5 * source_weight / length),
            (source_weight + 1, 0.5 * (length - source_weight) / length),
        )
        for destination, probability in contributions:
            if not probability:
                continue
            result[destination] = float(
                np.logaddexp(
                    result[destination], source_log + math.log(probability)
                )
            )
    return result


def odd_hypergeometric_probability(
    population: int, marked: int, draws: int
) -> float:
    denominator = math.comb(population, draws)
    numerator = sum(
        math.comb(marked, selected)
        * math.comb(population - marked, draws - selected)
        for selected in range(1, draws + 1, 2)
        if selected <= marked and draws - selected <= population - marked
    )
    return numerator / denominator


def expected_output_subset_spectrum(
    source_logs: np.ndarray, subset_size: int, target_size: int = 1
) -> np.ndarray:
    """Average a disjoint source/target rank-one shear over each shell.

    Setup samples a target set T and a disjoint source set S.  The forward
    map adds parity(x_S) to every coordinate in T.  Target size one is the
    original ParityShear construction.  The map is invertible because S and
    T are disjoint.
    """
    length = len(source_logs) - 1
    if not 1 <= subset_size < length or not 1 <= target_size < length:
        raise ValueError("source and target sizes must lie in [1, outer_bits)")
    if subset_size + target_size > length:
        raise ValueError("disjoint source and target sets do not fit")
    result = np.full(length + 1, -math.inf)
    result[0] = 0.0
    target_denominator = math.comb(length, target_size)
    for source_weight in range(1, length + 1):
        source_log = float(source_logs[source_weight])
        if not math.isfinite(source_log):
            continue
        for target_ones in range(target_size + 1):
            if target_ones > source_weight:
                continue
            target_zeros = target_size - target_ones
            if target_zeros > length - source_weight:
                continue
            target_probability = (
                math.comb(source_weight, target_ones)
                * math.comb(length - source_weight, target_zeros)
                / target_denominator
            )
            remaining_ones = source_weight - target_ones
            odd_probability = odd_hypergeometric_probability(
                length - target_size, remaining_ones, subset_size
            )
            contributions = (
                (source_weight, target_probability * (1.0 - odd_probability)),
                (
                    source_weight + target_size - 2 * target_ones,
                    target_probability * odd_probability,
                ),
            )
            for destination, probability in contributions:
                if probability <= 0.0:
                    continue
                result[destination] = float(
                    np.logaddexp(
                        result[destination], source_log + math.log(probability)
                    )
                )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-bits", type=int, default=256)
    parser.add_argument("--outer-dimension", type=int, default=128)
    parser.add_argument("--source-minimum-distance", type=int, default=38)
    parser.add_argument(
        "--mode",
        choices=("message-linear", "output-subset"),
        default="output-subset",
    )
    parser.add_argument("--subset-size", type=int, default=16)
    parser.add_argument("--target-size", type=int, default=1)
    parser.add_argument("--tail-endpoint-width", type=int, default=0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    source = modeled_even_floor_spectrum_logs(
        args.outer_bits, args.source_minimum_distance
    )
    transformed = (
        expected_message_linear_spectrum(source)
        if args.mode == "message-linear"
        else expected_output_subset_spectrum(
            source, args.subset_size, args.target_size
        )
    )
    mass_log2 = logsumexp(transformed) / LOG2
    envelope_log, envelope_weight = spectrum_bernoulli_envelope_log(
        args.outer_bits,
        transformed,
        0.5,
        args.tail_endpoint_width,
    )
    rows = [
        {
            "weight": weight,
            "log2_expected_multiplicity": value / LOG2,
        }
        for weight, value in enumerate(transformed)
        if math.isfinite(float(value))
    ]
    payload = {
        "schema": "riffle-parity-break-outer-spectrum-v1",
        "candidate": (
            "ParityBreak-MessageLinear"
            if args.mode == "message-linear"
            else (
                f"ParityBreak-OutputSubset-{args.subset_size}"
                f"-TargetSet-{args.target_size}"
            )
        ),
        "parameters": {
            "outer_bits": args.outer_bits,
            "outer_dimension": args.outer_dimension,
            "source_minimum_distance": args.source_minimum_distance,
            "transformed_minimum_distance": args.source_minimum_distance - 1,
            "tail_endpoint_width": args.tail_endpoint_width,
            "mode": args.mode,
            "subset_size": (
                args.subset_size if args.mode == "output-subset" else None
            ),
            "target_size": (
                args.target_size if args.mode == "output-subset" else None
            ),
        },
        "checks": {
            "expected_spectrum_log2_mass": mass_log2,
            "bernoulli_envelope_weight": envelope_weight,
            "bernoulli_envelope_bits_per_active_block": envelope_log / LOG2,
            "excess_envelope_bits_over_dimension": (
                envelope_log / LOG2 - args.outer_dimension
            ),
        },
        "spectrum": rows,
        "scope": (
            "Expected spectrum over an independently sampled parity-breaking "
            "shear for each outer-block instance. The zero word remains zero. "
            "MessageLinear samples a message-linear functional and a target. "
            "OutputSubset samples a target and a uniformly random subset of "
            "the remaining output coordinates, then XORs their parity into "
            "the target."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"mass_log2,{mass_log2:.12f}")
    print(f"minimum_distance,{args.source_minimum_distance - 1}")
    print(f"envelope_weight,{envelope_weight}")
    print(f"envelope_bits,{envelope_log / LOG2:.12f}")
    print(
        "excess_over_dimension_bits,"
        f"{envelope_log / LOG2 - args.outer_dimension:.12f}"
    )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
