#!/usr/bin/env python3
"""Analyze all messages through prefix ranks and a random Toeplitz inner.

The inner is lower-triangular Toeplitz convolution with h[0] = 1 and uniform
independent h[1],...,h[N-1]. It emits zero before the first nonzero input, one
at the activation coordinate, and independent fair bits afterward.

The script evaluates three objects:

* the exact first moment for independent uniform random outer injections;
* the exact first moment for independently sampled Golay--BA-3 outer rows;
* fixed setups with one sampled Golay--BA-3 code reused in every row.

All arithmetic is nearest binary64. Fixed setups use NumPy PCG64 permutations
and are reproducible diagnostics, not outward-rounded certificates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.special import gammaln, logsumexp

from analyze_golay_ba_rm2sub_joint import expected_ba_log_spectrum


WORKSTREAM = Path(__file__).resolve().parent
LN2 = math.log(2.0)
GOLAY_GENERATOR = sum(1 << i for i in (11, 9, 7, 6, 5, 1, 0))
GOLAY_SPECTRUM = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}


def log_two_power_minus_one(bits: int) -> float:
    if bits <= 0:
        return -math.inf
    value = bits * LN2
    if value > 50.0:
        return value + math.log1p(-math.exp(-value))
    return math.log(math.expm1(value))


def log_binomial_lower_tail(trials: int, threshold: int) -> float:
    """Return log Pr[Bin(trials,1/2) <= threshold]."""
    if threshold < 0:
        return -math.inf
    if threshold >= trials:
        return 0.0
    log_largest = (
        math.lgamma(trials + 1)
        - math.lgamma(threshold + 1)
        - math.lgamma(trials - threshold + 1)
        - trials * LN2
    )
    total = 1.0
    ratio_product = 1.0
    for weight in range(threshold, 0, -1):
        ratio_product *= weight / (trials - weight + 1)
        total += ratio_product
        if ratio_product < total * 1.0e-18:
            break
    return log_largest + math.log(total)


def golay_rows() -> list[int]:
    rows = []
    for message_bit in range(12):
        word = GOLAY_GENERATOR << message_bit
        word |= (word.bit_count() & 1) << 23
        rows.append(word)

    spectrum: dict[int, int] = {}
    for message in range(1 << 12):
        word = 0
        for bit, row in enumerate(rows):
            if (message >> bit) & 1:
                word ^= row
        spectrum[word.bit_count()] = spectrum.get(word.bit_count(), 0) + 1
    if spectrum != GOLAY_SPECTRUM:
        raise ArithmeticError(f"unexpected extended Golay spectrum: {spectrum}")
    return rows


def golay_direct_sum_columns(outer_bits: int) -> list[int]:
    if outer_bits % 24:
        raise ValueError("outer bits must be divisible by 24")
    rows = golay_rows()
    columns = []
    for block in range(outer_bits // 24):
        shift = 12 * block
        for coordinate in range(24):
            column = 0
            for message_bit, row in enumerate(rows):
                column |= ((row >> coordinate) & 1) << (shift + message_bit)
            columns.append(column)
    return columns


def apply_accumulator(columns: list[int], permutation: np.ndarray) -> list[int]:
    output = []
    state = 0
    for coordinate in permutation:
        state ^= columns[int(coordinate)]
        output.append(state)
    return output


def matrix_rank(columns: list[int], dimension: int) -> int:
    basis = [0] * dimension
    rank = 0
    for column in columns:
        value = column
        while value:
            pivot = value.bit_length() - 1
            if basis[pivot]:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                rank += 1
                break
    return rank


def sample_ba_columns(outer_bits: int, rng: np.random.Generator) -> list[int]:
    columns = golay_direct_sum_columns(outer_bits)
    for _ in range(2):
        columns = apply_accumulator(columns, rng.permutation(outer_bits))
    dimension = outer_bits // 2
    rank = matrix_rank(columns, dimension)
    if rank != dimension:
        raise ArithmeticError(f"BA generator rank {rank}, expected {dimension}")
    return columns


def sample_rank_increments(
    columns: list[int],
    outer_rows: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return increment[j,i] for row i after its j-th prior coordinate."""
    outer_bits = len(columns)
    dimension = outer_bits // 2
    increments = np.zeros((outer_bits, outer_rows), dtype=np.uint8)
    for row in range(outer_rows):
        basis = [0] * dimension
        rank = 0
        for position, coordinate in enumerate(rng.permutation(outer_bits)):
            if rank == dimension:
                break
            value = columns[int(coordinate)]
            while value:
                pivot = value.bit_length() - 1
                if basis[pivot]:
                    value ^= basis[pivot]
                else:
                    basis[pivot] = value
                    rank += 1
                    increments[position, row] = 1
                    break
        if rank != dimension:
            raise ArithmeticError(f"row {row} prefix order ended at rank {rank}")
    return increments


def fixed_prefix_dimensions(
    increments: np.ndarray,
    dimension: int,
    rng: np.random.Generator,
) -> np.ndarray:
    outer_bits, outer_rows = increments.shape
    output_bits = outer_bits * outer_rows
    prefix_dimension = np.empty(output_bits + 1, dtype=np.int32)
    current = dimension * outer_rows
    prefix_dimension[0] = current
    position = 0
    for region in range(outer_bits):
        for row in rng.permutation(outer_rows):
            current -= int(increments[region, int(row)])
            position += 1
            prefix_dimension[position] = current
    if current != 0:
        raise ArithmeticError(f"full routed generator has kernel dimension {current}")
    return prefix_dimension


def log_nonzero_counts(log_total_counts: np.ndarray) -> np.ndarray:
    result = np.full_like(log_total_counts, -math.inf)
    positive = log_total_counts > 0.0
    large = positive & (log_total_counts >= 50.0)
    result[large] = log_total_counts[large]
    small = positive & ~large
    result[small] = np.log(np.expm1(log_total_counts[small]))
    return result


def inner_first_moment_from_log_prefix_counts(
    log_nonzero_prefix: np.ndarray,
    output_bits: int,
    distance: int,
) -> dict[str, object]:
    """Evaluate the exact summation-by-parts identity for the ideal inner."""
    maximum_prefix = output_bits - distance
    if len(log_nonzero_prefix) < maximum_prefix + 1:
        raise ValueError("prefix-count array is too short")

    log_full_tail = log_binomial_lower_tail(output_bits - 1, distance - 1)
    base_log = float(log_nonzero_prefix[0]) + log_full_tail

    prefix = np.arange(1, maximum_prefix + 1, dtype=np.int64)
    suffix = output_bits - prefix
    log_delta = (
        gammaln(suffix)
        - gammaln(distance)
        - gammaln(suffix - distance + 1)
        - suffix * LN2
    )
    terms = log_nonzero_prefix[1 : maximum_prefix + 1] + log_delta
    aggregate_log = float(np.logaddexp(base_log, logsumexp(terms)))
    dominant_index = int(np.argmax(terms))
    dominant_prefix = dominant_index + 1

    return {
        "expected_bad_log2": aggregate_log / LN2,
        "margin_bits": -aggregate_log / LN2,
        "full_length_term_log2": base_log / LN2,
        "full_length_term_margin_bits": -base_log / LN2,
        "dominant_nonzero_prefix": dominant_prefix,
        "dominant_suffix_length": output_bits - dominant_prefix,
        "dominant_nonzero_prefix_count_log2": (
            float(log_nonzero_prefix[dominant_prefix]) / LN2
        ),
        "dominant_term_log2": float(terms[dominant_index]) / LN2,
    }


def fixed_setup_result(
    seed: int,
    outer_bits: int,
    outer_rows: int,
    distance: int,
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    columns = sample_ba_columns(outer_bits, rng)
    increments = sample_rank_increments(columns, outer_rows, rng)
    dimension = outer_bits // 2
    prefix_dimension = fixed_prefix_dimensions(increments, dimension, rng)
    log_prefix = np.asarray(
        [log_two_power_minus_one(int(value)) for value in prefix_dimension],
        dtype=np.float64,
    )
    result = inner_first_moment_from_log_prefix_counts(
        log_prefix,
        outer_bits * outer_rows,
        distance,
    )
    ideal_dimension = np.maximum(
        dimension * outer_rows - np.arange(len(prefix_dimension)),
        0,
    )
    deficit = prefix_dimension.astype(np.int64) - ideal_dimension
    nonzero_deficit = np.flatnonzero(deficit)
    result.update(
        {
            "seed": seed,
            "prefix_dimension_sha256": hashlib.sha256(
                prefix_dimension.tobytes()
            ).hexdigest(),
            "earliest_prefix_rank_deficit": (
                int(nonzero_deficit[0]) if len(nonzero_deficit) else None
            ),
            "maximum_prefix_rank_deficit": int(np.max(deficit)),
            "prefix_positions_with_rank_deficit": int(np.count_nonzero(deficit)),
            "sum_prefix_rank_deficit": int(np.sum(deficit, dtype=np.int64)),
            "dominant_region": int(result["dominant_nonzero_prefix"]) // outer_rows,
            "dominant_position_in_region": (
                int(result["dominant_nonzero_prefix"]) % outer_rows
            ),
            "rank_increment_counts_by_region_sha256": hashlib.sha256(
                np.sum(increments, axis=1, dtype=np.int32).tobytes()
            ).hexdigest(),
        }
    )
    return result


def random_injection_log_local_prefix_counts(
    outer_bits: int,
    dimension: int,
) -> np.ndarray:
    result = np.empty(outer_bits + 1, dtype=np.float64)
    log_nonzero_messages = log_two_power_minus_one(dimension)
    log_nonzero_words = log_two_power_minus_one(outer_bits)
    for constraints in range(outer_bits + 1):
        if constraints == 0:
            result[constraints] = dimension * LN2
        elif constraints == outer_bits:
            result[constraints] = 0.0
        else:
            log_surviving_nonzero = (
                log_nonzero_messages
                + log_two_power_minus_one(outer_bits - constraints)
                - log_nonzero_words
            )
            result[constraints] = float(np.logaddexp(0.0, log_surviving_nonzero))
    return result


def ba_log_local_prefix_counts(outer_bits: int) -> np.ndarray:
    log_spectrum = expected_ba_log_spectrum(outer_bits)
    result = np.empty(outer_bits + 1, dtype=np.float64)
    for constraints in range(outer_bits + 1):
        denominator = (
            gammaln(outer_bits + 1)
            - gammaln(constraints + 1)
            - gammaln(outer_bits - constraints + 1)
        )
        terms = []
        for weight in range(outer_bits - constraints + 1):
            if not math.isfinite(float(log_spectrum[weight])):
                continue
            log_avoid = (
                gammaln(outer_bits - weight + 1)
                - gammaln(constraints + 1)
                - gammaln(outer_bits - weight - constraints + 1)
                - denominator
            )
            terms.append(float(log_spectrum[weight]) + float(log_avoid))
        result[constraints] = float(logsumexp(terms))
    result[0] = (outer_bits // 2) * LN2
    result[-1] = 0.0
    return result


def independent_rows_result(
    name: str,
    log_local_prefix: np.ndarray,
    outer_rows: int,
    distance: int,
) -> dict[str, object]:
    outer_bits = len(log_local_prefix) - 1
    output_bits = outer_bits * outer_rows
    maximum_prefix = output_bits - distance
    prefix = np.arange(maximum_prefix + 1, dtype=np.int64)
    full_regions = prefix // outer_rows
    partial = prefix % outer_rows
    log_total = (
        (outer_rows - partial) * log_local_prefix[full_regions]
        + partial * log_local_prefix[full_regions + (partial > 0)]
    )
    log_nonzero = log_nonzero_counts(log_total)
    result = inner_first_moment_from_log_prefix_counts(
        log_nonzero,
        output_bits,
        distance,
    )
    dominant_prefix = int(result["dominant_nonzero_prefix"])
    result.update(
        {
            "name": name,
            "dominant_region": dominant_prefix // outer_rows,
            "dominant_position_in_region": dominant_prefix % outer_rows,
            "dominant_prefix_fraction": dominant_prefix / output_bits,
            "dominant_suffix_fraction": (
                output_bits - dominant_prefix
            ) / output_bits,
        }
    )
    return result


def ideal_prefix_result(
    outer_bits: int,
    outer_rows: int,
    distance: int,
) -> dict[str, object]:
    output_bits = outer_bits * outer_rows
    message_bits = (outer_bits // 2) * outer_rows
    maximum_prefix = output_bits - distance
    dimensions = np.maximum(
        message_bits - np.arange(maximum_prefix + 1, dtype=np.int64),
        0,
    )
    log_prefix = np.asarray(
        [log_two_power_minus_one(int(value)) for value in dimensions],
        dtype=np.float64,
    )
    result = inner_first_moment_from_log_prefix_counts(
        log_prefix,
        output_bits,
        distance,
    )
    result["name"] = "ideal full-rank prefix benchmark"
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-bits", type=int, default=720)
    parser.add_argument("--outer-rows", type=int, default=2944)
    parser.add_argument("--relative-distance", type=float, default=0.11)
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=[20260902, 20260903, 20260904, 20260905],
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=WORKSTREAM / "ba3_B720_ideal_causal_g1_prefix_diagnostic.json",
    )
    parser.add_argument(
        "--skip-fixed",
        action="store_true",
        help="evaluate only the independent-row analytic baselines",
    )
    args = parser.parse_args()

    if args.outer_bits % 24:
        raise ValueError("outer bits must be divisible by 24")
    output_bits = args.outer_bits * args.outer_rows
    parent_message_bits = (args.outer_bits // 2) * args.outer_rows
    distance = math.floor(args.relative_distance * output_bits)

    random_local = random_injection_log_local_prefix_counts(
        args.outer_bits,
        args.outer_bits // 2,
    )
    ba_local = ba_log_local_prefix_counts(args.outer_bits)
    random_result = independent_rows_result(
        "independent uniform random injections",
        random_local,
        args.outer_rows,
        distance,
    )
    independent_ba_result = independent_rows_result(
        "independently sampled Golay--BA-3 rows",
        ba_local,
        args.outer_rows,
        distance,
    )
    benchmark_result = ideal_prefix_result(
        args.outer_bits,
        args.outer_rows,
        distance,
    )

    fixed_results = []
    if not args.skip_fixed:
        for seed in args.seeds:
            print(f"evaluating fixed repeated-BA setup seed {seed}", flush=True)
            fixed_results.append(
                fixed_setup_result(
                    seed,
                    args.outer_bits,
                    args.outer_rows,
                    distance,
                )
            )

    result = {
        "schema": "ba3-b720-ideal-causal-g1-prefix-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "outer_bits": args.outer_bits,
            "outer_dimension": args.outer_bits // 2,
            "outer_rows": args.outer_rows,
            "parent_message_bits": parent_message_bits,
            "target_message_bits": 2**20,
            "zero_shortened_input_bits": parent_message_bits - 2**20,
            "output_bits": output_bits,
            "bad_output_weight_at_most": distance,
            "relative_distance_target": args.relative_distance,
            "packet_bits": 1,
            "inner": "invertible random lower-triangular Toeplitz convolution",
        },
        "inner_law": [
            "Sample h[1],...,h[N-1] independently and uniformly, set h[0]=1, and compute y[t]=sum_{i=0}^t h[i]u[t-i].",
            "The output is zero before the first nonzero routed input bit.",
            "The activation output equals one.",
            "For each fixed nonzero input difference, every later output bit is independent Bernoulli one half because the suffix map from fresh h coefficients is triangular.",
        ],
        "probability_spaces": {
            "independent_uniform_outer": [
                "one independent uniform linear injection for each outer row",
                "independent local coordinate and region permutations",
                "one random Toeplitz convolution shared by all messages",
            ],
            "independent_golay_ba3_outer": [
                "two independent accumulator interleavers for each outer row",
                "independent local coordinate and region permutations",
                "one random Toeplitz convolution shared by all messages",
            ],
            "fixed_repeated_golay_ba3": [
                "the recorded PCG64 seed deterministically fixes two accumulator interleavers, one repeated BA generator, and all routing permutations",
                "the only remaining randomness in the displayed conditional first moment is the one Toeplitz convolution shared by all messages",
            ],
        },
        "independent_row_ensemble_results": {
            "uniform_random_outer": random_result,
            "golay_ba3_outer": independent_ba_result,
        },
        "ideal_prefix_benchmark": benchmark_result,
        "fixed_repeated_ba_setups": fixed_results,
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "rng": "numpy.random.Generator(PCG64)",
        },
        "interpretation": [
            "The independent-row calculations are exact in form because the local prefix-count expectations factor across independently sampled rows.",
            "The fixed repeated-BA calculations condition on one sampled BA generator and all sampled routing permutations; they use actual prefix ranks and do not multiply BA spectrum expectations.",
            "The fixed-setup inner first moment covers all nonzero parent messages at once; it is not an occupation truncation.",
            "The repeated-BA ensemble average cannot be reconstructed from the expected BA spectrum because it requires high moments of local prefix counts.",
        ],
        "limitations": [
            "All arithmetic is nearest binary64.",
            "The PCG64 fixed setups are reproducible diagnostics, not frozen cryptographic setup specifications.",
            "Direct Toeplitz convolution has unbounded memory; fast polynomial multiplication is not the current linear-time RM2Sub-S19 implementation.",
            "No probability bound on obtaining a fixed setup with the displayed prefix profile is proved.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "random_outer_margin_bits": random_result["margin_bits"],
                "independent_ba_margin_bits": independent_ba_result["margin_bits"],
                "fixed_repeated_ba_margins": [
                    row["margin_bits"] for row in fixed_results
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
