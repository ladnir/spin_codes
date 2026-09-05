#!/usr/bin/env python3
"""Local constituent audit for Riffle LDPCSplitState t=256 s=64.

The audit has two independent parts.  First, it computes the exact ensemble
moment of four independently interleaved length-64 accumulators applied to a
uniform nonzero state.  Second, it samples sparse 64-by-256 compression
matrices with unique fixed-weight columns and counts kernel words through
weight four exactly.

These calculations do not constitute an end-to-end distance certificate.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path

import numpy as np

from explore_riffle_packet4_longacc_fieldchecksum import (
    accumulator_weight_counts,
)


DEFAULT_OUTPUT = Path(
    "constructions/riffle_ldpcsplitstate_g4_t256_s64/"
    "receipts/local_structure_audit.json"
)


def gf2_rank(columns: list[int], rows: int) -> int:
    basis = [0] * rows
    rank = 0
    for value in columns:
        x = value
        while x:
            pivot = x.bit_length() - 1
            if basis[pivot]:
                x ^= basis[pivot]
            else:
                basis[pivot] = x
                rank += 1
                break
    return rank


def sample_unique_columns(
    *, rows: int, columns: int, column_weight: int, rng: random.Random
) -> list[int]:
    result: set[int] = set()
    while len(result) < columns:
        positions = rng.sample(range(rows), column_weight)
        mask = sum(1 << position for position in positions)
        result.add(mask)
    return list(result)


def sample_unique_mixed_columns(
    *, rows: int, degrees: list[int], rng: random.Random
) -> list[int]:
    result: set[int] = set()
    columns = []
    shuffled_degrees = degrees.copy()
    rng.shuffle(shuffled_degrees)
    for degree in shuffled_degrees:
        while True:
            positions = rng.sample(range(rows), degree)
            mask = sum(1 << position for position in positions)
            if mask not in result:
                result.add(mask)
                columns.append(mask)
                break
    return columns


def low_kernel_counts(columns: list[int]) -> dict[str, int]:
    """Count distinct kernel supports of weights two through four."""
    column_set = set(columns)
    column_index = {value: index for index, value in enumerate(columns)}
    weight_two = len(columns) - len(column_set)

    weight_three_hits: set[tuple[int, int, int]] = set()
    pairs_by_syndrome: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for left in range(len(columns)):
        for right in range(left + 1, len(columns)):
            syndrome = columns[left] ^ columns[right]
            if syndrome in column_set:
                third = column_index[syndrome]
                if third != left and third != right:
                    weight_three_hits.add(tuple(sorted((left, right, third))))
            pairs_by_syndrome[syndrome].append((left, right))

    weight_four_hits: set[tuple[int, int, int, int]] = set()
    for pairs in pairs_by_syndrome.values():
        for first in range(len(pairs)):
            a, b = pairs[first]
            for second in range(first + 1, len(pairs)):
                c, d = pairs[second]
                if len({a, b, c, d}) == 4:
                    weight_four_hits.add(tuple(sorted((a, b, c, d))))
    return {
        "weight_2": weight_two,
        "weight_3": len(weight_three_hits),
        "weight_4": len(weight_four_hits),
    }


def sparse_compression_scan(
    *, rows: int, columns: int, degrees: tuple[int, ...], samples: int, seed: int
) -> list[dict[str, object]]:
    rng = random.Random(seed)
    results = []
    for degree in degrees:
        best: dict[str, object] | None = None
        full_rank_samples = 0
        for sample_index in range(samples):
            candidate = sample_unique_columns(
                rows=rows,
                columns=columns,
                column_weight=degree,
                rng=rng,
            )
            rank = gf2_rank(candidate, rows)
            counts = low_kernel_counts(candidate)
            if rank == rows:
                full_rank_samples += 1
            score = (
                int(rank != rows),
                counts["weight_2"],
                counts["weight_3"],
                counts["weight_4"],
            )
            if best is None or score < tuple(best["score"]):
                best = {
                    "sample_index": sample_index,
                    "rank": rank,
                    "counts": counts,
                    "score": list(score),
                    "columns_hex": [f"{value:016x}" for value in candidate],
                }
        assert best is not None
        results.append(
            {
                "column_weight": degree,
                "samples": samples,
                "full_rank_samples": full_rank_samples,
                "best": best,
                "xor_contributions_per_output": float(degree),
            }
        )
    return results


def mixed_compression_scan(
    *, rows: int, columns: int, odd_fractions: tuple[float, ...], samples: int, seed: int
) -> list[dict[str, object]]:
    rng = random.Random(seed ^ 0x233233)
    results = []
    for odd_fraction in odd_fractions:
        degree_three = round(columns * odd_fraction)
        degrees = [3] * degree_three + [2] * (columns - degree_three)
        best: dict[str, object] | None = None
        full_rank_samples = 0
        for sample_index in range(samples):
            candidate = sample_unique_mixed_columns(
                rows=rows, degrees=degrees, rng=rng
            )
            rank = gf2_rank(candidate, rows)
            counts = low_kernel_counts(candidate)
            if rank == rows:
                full_rank_samples += 1
            score = (
                int(rank != rows),
                counts["weight_2"],
                counts["weight_3"],
                counts["weight_4"],
            )
            if best is None or score < tuple(best["score"]):
                best = {
                    "sample_index": sample_index,
                    "rank": rank,
                    "counts": counts,
                    "score": list(score),
                    "columns_hex": [f"{value:016x}" for value in candidate],
                }
        assert best is not None
        results.append(
            {
                "degree_three_fraction": odd_fraction,
                "average_column_weight": sum(degrees) / columns,
                "samples": samples,
                "full_rank_samples": full_rank_samples,
                "best": best,
            }
        )
    return results


def accumulator_state_audit(
    *, state_bits: int, lanes: int, target_distance: int, epochs: int
) -> dict[str, object]:
    distributions = [
        accumulator_weight_counts(state_bits, weight)
        for weight in range(state_bits + 1)
    ]
    denominators = [math.comb(state_bits, weight) for weight in range(state_bits + 1)]

    local_threshold = math.floor(0.09 * state_bits * lanes)
    local_low_probability = 0.0
    # Exact fourfold convolution is used only for the low local threshold.
    for state_weight in range(1, state_bits + 1):
        total = {0: 1}
        for _ in range(lanes):
            following: dict[int, int] = defaultdict(int)
            for left_weight, left_count in total.items():
                for right_weight, right_count in distributions[state_weight].items():
                    following[left_weight + right_weight] += left_count * right_count
            total = dict(following)
        low = sum(count for weight, count in total.items() if weight <= local_threshold)
        probability_of_state_weight = (
            math.comb(state_bits, state_weight) / ((1 << state_bits) - 1)
        )
        local_low_probability += probability_of_state_weight * (
            low / denominators[state_weight] ** lanes
        )

    grid = (
        -10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
        -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0,
        -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5,
        0.75, 1.0,
    )
    best_log2 = math.inf
    best_tilt = math.nan
    rows = []
    for log_surprisal in grid:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        live_moment = 0.0
        for state_weight in range(1, state_bits + 1):
            lane_moment = sum(
                count * z**weight
                for weight, count in distributions[state_weight].items()
            ) / denominators[state_weight]
            live_moment += (
                math.comb(state_bits, state_weight)
                * lane_moment**lanes
                / ((1 << state_bits) - 1)
            )
        log2_bound = (
            epochs * math.log2(live_moment)
            + target_distance * surprisal / math.log(2.0)
        )
        rows.append(
            {
                "log_surprisal": log_surprisal,
                "live_moment_log2": math.log2(live_moment),
                "repeated_target_log2_bound": log2_bound,
            }
        )
        if log2_bound < best_log2:
            best_log2 = log2_bound
            best_tilt = log_surprisal
    return {
        "state_bits": state_bits,
        "lanes": lanes,
        "output_bits_per_epoch": state_bits * lanes,
        "local_distance_threshold": local_threshold,
        "local_low_weight_probability": local_low_probability,
        "local_low_weight_bits": -math.log2(local_low_probability),
        "repeated_live_inner_log2_bound": best_log2,
        "repeated_live_inner_suppression_bits": -best_log2,
        "best_log_surprisal": best_tilt,
        "tilt_rows": rows,
        "scope": (
            "ensemble average over four independent state interleavers; "
            "zero-state, affine-coset, and termination transitions omitted"
        ),
    }


def odd_intersection_probability(
    *, population: int, selected: int, marked: int
) -> float:
    numerator = 0
    for intersection in range(1, selected + 1, 2):
        if intersection <= marked and selected - intersection <= population - marked:
            numerator += math.comb(marked, intersection) * math.comb(
                population - marked, selected - intersection
            )
    return numerator / math.comb(population, selected)


def repeated_live_bound(
    *,
    state_bits: int,
    output_bits: int,
    epochs: int,
    target_distance: int,
    one_bit_probabilities: list[float],
) -> dict[str, object]:
    grid = (
        -10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
        -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0,
        -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5,
        0.75, 1.0,
    )
    best_log2 = math.inf
    best_tilt = math.nan
    rows = []
    denominator = (1 << state_bits) - 1
    for log_surprisal in grid:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        live_moment = 0.0
        for state_weight in range(1, state_bits + 1):
            probability_one = one_bit_probabilities[state_weight]
            output_moment = (
                1.0 - probability_one + probability_one * z
            ) ** output_bits
            live_moment += (
                math.comb(state_bits, state_weight)
                * output_moment
                / denominator
            )
        log2_bound = (
            epochs * math.log2(live_moment)
            + target_distance * surprisal / math.log(2.0)
        )
        rows.append(
            {
                "log_surprisal": log_surprisal,
                "live_moment_log2": math.log2(live_moment),
                "repeated_target_log2_bound": log2_bound,
            }
        )
        if log2_bound < best_log2:
            best_log2 = log2_bound
            best_tilt = log_surprisal
    return {
        "repeated_live_inner_log2_bound": best_log2,
        "repeated_live_inner_suppression_bits": -best_log2,
        "best_log_surprisal": best_tilt,
        "tilt_rows": rows,
    }


def sparse_expansion_audit(
    *,
    state_bits: int,
    output_bits: int,
    row_degrees: tuple[int, ...],
    epochs: int,
    target_distance: int,
) -> list[dict[str, object]]:
    local_threshold = math.floor(0.09 * output_bits)
    denominator = (1 << state_bits) - 1
    rows = []
    for degree in row_degrees:
        probabilities = [
            odd_intersection_probability(
                population=state_bits,
                selected=degree,
                marked=state_weight,
            )
            for state_weight in range(state_bits + 1)
        ]
        local_low_probability = 0.0
        for state_weight in range(1, state_bits + 1):
            p = probabilities[state_weight]
            low_probability = sum(
                math.comb(output_bits, weight)
                * p**weight
                * (1.0 - p) ** (output_bits - weight)
                for weight in range(local_threshold + 1)
            )
            local_low_probability += (
                math.comb(state_bits, state_weight)
                * low_probability
                / denominator
            )
        repeated = repeated_live_bound(
            state_bits=state_bits,
            output_bits=output_bits,
            epochs=epochs,
            target_distance=target_distance,
            one_bit_probabilities=probabilities,
        )
        rows.append(
            {
                "row_weight": degree,
                "xor_per_output": degree - 1,
                "expected_basis_codeword_weight": output_bits * degree / state_bits,
                "local_low_weight_bits": -math.log2(local_low_probability),
                **repeated,
                "scope": (
                    "ensemble average with each output row sampled "
                    "independently and uniformly at the stated weight; "
                    "the constraint BA=0 is omitted"
                ),
            }
        )
    return rows


def dense_random_expansion_reference(
    *, output_bits: int, epochs: int, target_distance: int
) -> dict[str, float]:
    best_log2 = math.inf
    best_tilt = math.nan
    for log_surprisal in (
        -10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
        -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0,
        -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5,
        0.75, 1.0,
    ):
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        moment = ((1.0 + z) ** output_bits - 1.0) / (
            math.ldexp(1.0, output_bits) - 1.0
        )
        candidate = (
            epochs * math.log2(moment)
            + target_distance * surprisal / math.log(2.0)
        )
        if candidate < best_log2:
            best_log2 = candidate
            best_tilt = log_surprisal
    return {
        "repeated_live_inner_log2_bound": best_log2,
        "repeated_live_inner_suppression_bits": -best_log2,
        "best_log_surprisal": best_tilt,
    }


def accumulator_moments_all_input_weights(length: int, z: float) -> np.ndarray:
    """Return sum_y z^wt(y) for every accumulator input weight."""
    current = np.zeros((length + 1, 2), dtype=np.float64)
    current[0, 0] = 1.0
    maximum_used = 0
    for _ in range(length):
        following = np.zeros_like(current)
        # Input zero preserves the accumulator state.
        following[: maximum_used + 1, 0] += current[: maximum_used + 1, 0]
        following[: maximum_used + 1, 1] += z * current[: maximum_used + 1, 1]
        # Input one toggles the accumulator state.
        following[1 : maximum_used + 2, 1] += z * current[: maximum_used + 1, 0]
        following[1 : maximum_used + 2, 0] += current[: maximum_used + 1, 1]
        current = following
        maximum_used += 1
    return np.sum(current, axis=1)


def repeat_accumulate_expansion_audit(
    *,
    state_bits: int,
    repetition: int,
    epochs: int,
    target_distance: int,
) -> dict[str, object]:
    output_bits = state_bits * repetition
    denominator = (1 << state_bits) - 1
    grid = (
        -10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
        -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0,
        -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5,
        0.75, 1.0,
    )
    best_log2 = math.inf
    best_tilt = math.nan
    rows = []
    for log_surprisal in grid:
        surprisal = math.exp(log_surprisal)
        z = math.exp(-surprisal)
        unnormalized = accumulator_moments_all_input_weights(output_bits, z)
        live_moment = 0.0
        for state_weight in range(1, state_bits + 1):
            repeated_weight = repetition * state_weight
            conditional_moment = unnormalized[repeated_weight] / math.comb(
                output_bits, repeated_weight
            )
            live_moment += (
                math.comb(state_bits, state_weight)
                * conditional_moment
                / denominator
            )
        log2_bound = (
            epochs * math.log2(live_moment)
            + target_distance * surprisal / math.log(2.0)
        )
        rows.append(
            {
                "log_surprisal": log_surprisal,
                "live_moment_log2": math.log2(live_moment),
                "repeated_target_log2_bound": log2_bound,
            }
        )
        if log2_bound < best_log2:
            best_log2 = log2_bound
            best_tilt = log_surprisal
    return {
        "state_bits": state_bits,
        "repetition": repetition,
        "output_bits": output_bits,
        "accumulator_xor_per_output": 1.0,
        "repeated_live_inner_log2_bound": best_log2,
        "repeated_live_inner_suppression_bits": -best_log2,
        "best_log_surprisal": best_tilt,
        "tilt_rows": rows,
        "scope": (
            "ensemble average over a fresh uniform interleaver of the four "
            "state copies; affine input translates and fixed-interleaver "
            "concentration omitted"
        ),
    }


def accumulator_weight_transition(length: int) -> np.ndarray:
    """Return P[input weight, output weight] for a random input slice."""
    current = np.zeros((length + 1, 2, length + 1), dtype=np.float64)
    current[0, 0, 0] = 1.0
    maximum_used = 0
    maximum_weight = 0
    for _ in range(length):
        following = np.zeros_like(current)
        used = slice(0, maximum_used + 1)
        weights = slice(0, maximum_weight + 1)
        # Input zero, old state zero.
        following[used, 0, weights] += current[used, 0, weights]
        # Input zero, old state one: emit one.
        following[0 : maximum_used + 1, 1, 1 : maximum_weight + 2] += current[
            used, 1, weights
        ]
        # Input one, old state zero: toggle to one and emit one.
        following[1 : maximum_used + 2, 1, 1 : maximum_weight + 2] += current[
            used, 0, weights
        ]
        # Input one, old state one: toggle to zero and emit zero.
        following[1 : maximum_used + 2, 0, 0 : maximum_weight + 1] += current[
            used, 1, weights
        ]
        current = following
        maximum_used += 1
        maximum_weight += 1
    counts = np.sum(current, axis=1)
    for input_weight in range(length + 1):
        counts[input_weight] /= math.comb(length, input_weight)
    return counts


def multi_accumulate_expansion_audit(
    *,
    state_bits: int,
    repetition: int,
    accumulator_depths: tuple[int, ...],
    epochs: int,
    target_distance: int,
) -> list[dict[str, object]]:
    output_bits = state_bits * repetition
    transition = accumulator_weight_transition(output_bits)
    denominator = (1 << state_bits) - 1
    local_threshold = math.floor(0.09 * output_bits)
    state_input = np.zeros(output_bits + 1, dtype=np.float64)
    for state_weight in range(1, state_bits + 1):
        state_input[repetition * state_weight] = (
            math.comb(state_bits, state_weight) / denominator
        )

    distribution = state_input
    rows = []
    for depth in range(1, max(accumulator_depths) + 1):
        distribution = distribution @ transition
        if depth not in accumulator_depths:
            continue
        local_low_probability = float(
            np.sum(distribution[: local_threshold + 1])
        )
        best_log2 = math.inf
        best_tilt = math.nan
        for log_surprisal in (
            -10.0, -9.0, -8.5, -8.0, -7.5, -7.0, -6.5, -6.0,
            -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0,
            -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5,
            0.75, 1.0,
        ):
            surprisal = math.exp(log_surprisal)
            z = math.exp(-surprisal)
            powers = z ** np.arange(output_bits + 1)
            live_moment = float(distribution @ powers)
            candidate = (
                epochs * math.log2(live_moment)
                + target_distance * surprisal / math.log(2.0)
            )
            if candidate < best_log2:
                best_log2 = candidate
                best_tilt = log_surprisal
        rows.append(
            {
                "accumulator_depth": depth,
                "accumulator_xor_per_output": float(depth),
                "local_low_weight_bits": -math.log2(local_low_probability),
                "repeated_live_inner_log2_bound": best_log2,
                "repeated_live_inner_suppression_bits": -best_log2,
                "best_log_surprisal": best_tilt,
                "scope": (
                    "ensemble average over an independent uniform "
                    "interleaver before every accumulator"
                ),
            }
        )
    return rows


def evaluate(args: argparse.Namespace) -> dict[str, object]:
    target_distance = math.floor(
        args.relative_distance * 2 * args.message_bits
    )
    epochs = (2 * args.message_bits) // args.step_bits
    return {
        "schema": "riffle-ldpcsplitstate-local-v1",
        "candidate": "Riffle LDPCSplitState g=4 t=256 s=64",
        "parameters": {
            "state_bits": args.state_bits,
            "step_bits": args.step_bits,
            "message_bits": args.message_bits,
            "relative_distance": args.relative_distance,
            "compression_samples": args.compression_samples,
            "seed": args.seed,
        },
        "accumulator_state_expansion": accumulator_state_audit(
            state_bits=args.state_bits,
            lanes=args.step_bits // args.state_bits,
            target_distance=target_distance,
            epochs=epochs,
        ),
        "sparse_state_expansion": sparse_expansion_audit(
            state_bits=args.state_bits,
            output_bits=args.step_bits,
            row_degrees=tuple(args.expansion_row_degrees),
            epochs=epochs,
            target_distance=target_distance,
        ),
        "dense_random_expansion_reference": dense_random_expansion_reference(
            output_bits=args.step_bits,
            epochs=epochs,
            target_distance=target_distance,
        ),
        "repeat_accumulate_state_expansion": repeat_accumulate_expansion_audit(
            state_bits=args.state_bits,
            repetition=args.step_bits // args.state_bits,
            epochs=epochs,
            target_distance=target_distance,
        ),
        "multi_accumulate_state_expansion": multi_accumulate_expansion_audit(
            state_bits=args.state_bits,
            repetition=args.step_bits // args.state_bits,
            accumulator_depths=(1, 2, 3),
            epochs=epochs,
            target_distance=target_distance,
        ),
        "sparse_compression": sparse_compression_scan(
            rows=args.state_bits,
            columns=args.step_bits,
            degrees=tuple(args.compression_degrees),
            samples=args.compression_samples,
            seed=args.seed,
        ),
        "mixed_sparse_compression": mixed_compression_scan(
            rows=args.state_bits,
            columns=args.step_bits,
            odd_fractions=(0.125, 0.25, 0.5),
            samples=args.compression_samples,
            seed=args.seed,
        ),
        "scope": (
            "local constituent diagnostics only; no outer spectrum, packet "
            "placement, affine-coset maximum, or end-to-end union bound"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-bits", type=int, default=64)
    parser.add_argument("--step-bits", type=int, default=256)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--relative-distance", type=float, default=0.09)
    parser.add_argument(
        "--compression-degrees", type=int, nargs="+", default=(3, 5)
    )
    parser.add_argument(
        "--expansion-row-degrees",
        type=int,
        nargs="+",
        default=(2, 3, 4, 5, 6, 8, 10),
    )
    parser.add_argument("--compression-samples", type=int, default=16)
    parser.add_argument("--seed", type=int, default=0x4C445043)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = evaluate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    expansion = payload["accumulator_state_expansion"]
    print(
        "accumulator_live_suppression_bits,"
        f"{expansion['repeated_live_inner_suppression_bits']:.6f}"
    )
    print(
        "accumulator_local_low_weight_bits,"
        f"{expansion['local_low_weight_bits']:.6f}"
    )
    for row in payload["sparse_compression"]:
        counts = row["best"]["counts"]
        print(
            "compression_degree,"
            f"{row['column_weight']},rank,{row['best']['rank']},"
            f"w3,{counts['weight_3']},w4,{counts['weight_4']}"
        )
    for row in payload["mixed_sparse_compression"]:
        counts = row["best"]["counts"]
        print(
            "mixed_compression_average_degree,"
            f"{row['average_column_weight']:.3f},rank,{row['best']['rank']},"
            f"w3,{counts['weight_3']},w4,{counts['weight_4']}"
        )
    for row in payload["sparse_state_expansion"]:
        print(
            "expansion_row_weight,"
            f"{row['row_weight']},xor_per_output,{row['xor_per_output']},"
            "live_suppression_bits,"
            f"{row['repeated_live_inner_suppression_bits']:.6f},"
            f"local_low_bits,{row['local_low_weight_bits']:.6f}"
        )
    dense = payload["dense_random_expansion_reference"]
    print(
        "dense_random_live_suppression_bits,"
        f"{dense['repeated_live_inner_suppression_bits']:.6f}"
    )
    repeat_accumulate = payload["repeat_accumulate_state_expansion"]
    print(
        "repeat_accumulate_live_suppression_bits,"
        f"{repeat_accumulate['repeated_live_inner_suppression_bits']:.6f}"
    )
    for row in payload["multi_accumulate_state_expansion"]:
        print(
            "multi_accumulate_depth,"
            f"{row['accumulator_depth']},live_suppression_bits,"
            f"{row['repeated_live_inner_suppression_bits']:.6f},"
            f"local_low_bits,{row['local_low_weight_bits']:.6f}"
        )
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
