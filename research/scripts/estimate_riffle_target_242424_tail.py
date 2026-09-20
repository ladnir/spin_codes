#!/usr/bin/env python3
"""Estimate the full-size (24,24,24) shell with exact gap integration."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import log_ndtr, logsumexp

from analyze_riffle_parity_ladder_zero_accumulator import (
    GAP_PERIOD,
    GapCounter,
    compositions,
)
from analyze_riffle_shiftalpha64_lowblock_kernel import packetize, sample_weight_slice


PACKET_WIDTH = 4
PART_COUNT = 3
PART_WEIGHT = 24
PACKET_POSITIONS = 524_352
BINARY_LENGTH = PACKET_WIDTH * PACKET_POSITIONS
EXACT_OUTER_WORDS = 60_462_011
DEFAULT_RELATIVE_DISTANCES = (
    0.08,
    0.09,
    0.10,
    0.11,
    0.12,
    0.13,
    0.14,
    0.15,
    0.16,
    0.17,
    0.18,
    0.19,
    0.20,
    0.21,
    0.22,
)


def log_integer(value: int) -> float:
    if value <= 0:
        return -math.inf
    top_bit = value.bit_length() - 1
    return top_bit * math.log(2.0) + math.log(value / (1 << top_bit))


class FlexibleGapCounter(GapCounter):
    def __init__(self) -> None:
        self._q_cache: dict[tuple[tuple[int, ...], int], int] = {}
        self._interpolation_cache: dict[
            tuple[tuple[int, ...], int], tuple[int, ...]
        ] = {}

    @staticmethod
    @lru_cache(maxsize=None)
    def derivative_compositions(zero_visits: int) -> tuple[tuple[int, ...], ...]:
        return tuple(
            counts
            for degree in range(zero_visits + 1)
            for counts in compositions(degree, PACKET_WIDTH)
        )

    def positive_gap_count(self, exponents: tuple[int, ...], maximum_cost: int) -> int:
        if maximum_cost < 0:
            return 0
        if not any(exponents):
            return 1
        degree = sum(exponents)
        residue = maximum_cost % GAP_PERIOD
        key = (exponents, residue)
        differences = self._interpolation_cache.get(key)
        if differences is None:
            interpolation_stop = residue + GAP_PERIOD * degree
            coefficients = [0] * (interpolation_stop + 1)
            coefficients[0] = 1
            for weight, multiplicity in enumerate(exponents, 1):
                for _ in range(multiplicity):
                    for index in range(weight, interpolation_stop + 1):
                        coefficients[index] += coefficients[index - weight]
            for index in range(1, interpolation_stop + 1):
                coefficients[index] += coefficients[index - 1]
            values = [
                coefficients[residue + GAP_PERIOD * sample]
                for sample in range(degree + 1)
            ]
            first_differences = []
            while values:
                first_differences.append(values[0])
                values = [
                    values[index + 1] - values[index]
                    for index in range(len(values) - 1)
                ]
            differences = tuple(first_differences)
            self._interpolation_cache[key] = differences
        quotient = maximum_cost // GAP_PERIOD
        return sum(
            difference * math.comb(quotient, order)
            for order, difference in enumerate(differences)
        )

    def all_gap_count(
        self,
        visits: tuple[int, ...],
        inactive_packets: int,
        maximum_extra_weight: int,
    ) -> int:
        zero_visits = visits[0]
        positive_visits = visits[1:]
        result = 0
        for derivative_counts in self.derivative_compositions(zero_visits):
            derivative_degree = sum(derivative_counts)
            coefficient = (-1) ** derivative_degree * math.comb(
                inactive_packets + zero_visits - derivative_degree,
                zero_visits - derivative_degree,
            )
            for multiplicity, derivative_count in zip(
                positive_visits, derivative_counts
            ):
                if not derivative_count:
                    continue
                if not multiplicity:
                    coefficient = 0
                    break
                coefficient *= math.comb(
                    multiplicity + derivative_count - 1,
                    derivative_count,
                )
            if not coefficient:
                continue
            shift = sum(
                weight * derivative_counts[weight - 1]
                for weight in range(1, PACKET_WIDTH + 1)
            )
            exponents = tuple(
                multiplicity + derivative_count
                for multiplicity, derivative_count in zip(
                    positive_visits, derivative_counts
                )
            )
            result += coefficient * self.positive_gap_count(
                exponents,
                maximum_extra_weight - shift,
            )
        return result


def sample_visit_histograms(samples: int, seed: int) -> Counter[tuple[int, ...]]:
    rng = np.random.default_rng(seed)
    result: Counter[tuple[int, ...]] = Counter()
    for _ in range(samples):
        packets = np.concatenate(
            [
                packetize(sample_weight_slice(rng, PART_WEIGHT))
                for _ in range(PART_COUNT)
            ]
        )
        active = packets[packets != 0]
        rng.shuffle(active)
        state = 0
        visits = [0] * (PACKET_WIDTH + 1)
        for packet in active:
            state ^= int(packet)
            visits[state.bit_count()] += 1
        result[tuple(visits)] += 1
    return result


def estimate_distance(
    grouped: Counter[tuple[int, ...]],
    samples: int,
    distance: int,
    counter: FlexibleGapCounter,
    gap_method: str,
) -> dict[str, object]:
    log_terms = []
    log_square_terms = []
    nonzero_sample_paths = 0
    maximum_conditional_log = -math.inf
    for visits, frequency in grouped.items():
        active_packets = sum(visits)
        active_output_weight = sum(
            weight * visits[weight] for weight in range(1, PACKET_WIDTH + 1)
        )
        if gap_method == "exact":
            numerator = counter.all_gap_count(
                visits,
                PACKET_POSITIONS - active_packets,
                distance - active_output_weight,
            )
            if not numerator:
                continue
            log_probability = log_integer(numerator) - (
                math.lgamma(PACKET_POSITIONS + 1)
                - math.lgamma(active_packets + 1)
                - math.lgamma(PACKET_POSITIONS - active_packets + 1)
            )
        else:
            log_probability = saddle_gap_log_probability(
                visits,
                distance,
                active_output_weight,
            )
            if not math.isfinite(log_probability):
                continue
        nonzero_sample_paths += frequency
        maximum_conditional_log = max(maximum_conditional_log, log_probability)
        log_terms.append(math.log(frequency) + log_probability)
        log_square_terms.append(math.log(frequency) + 2.0 * log_probability)

    log_mean = float(logsumexp(log_terms)) - math.log(samples)
    log_second_moment = float(logsumexp(log_square_terms)) - math.log(samples)
    mean = math.exp(log_mean)
    second_moment = math.exp(log_second_moment)
    variance = max(0.0, second_moment - mean * mean)
    relative_standard_error = (
        math.sqrt(variance / samples) / mean if mean > 0.0 else math.inf
    )
    expected_log2 = (log_mean + math.log(EXACT_OUTER_WORDS)) / math.log(2.0)
    return {
        "distance": distance,
        "relative_binary_weight": distance / BINARY_LENGTH,
        "estimated_per_word_log2_probability": log_mean / math.log(2.0),
        "estimated_shell_expected_count_log2": expected_log2,
        "estimated_shell_expected_count": math.exp(
            log_mean + math.log(EXACT_OUTER_WORDS)
        ),
        "relative_standard_error": relative_standard_error,
        "sample_paths_with_nonzero_conditional_tail": nonzero_sample_paths,
        "maximum_sampled_conditional_log2_probability": (
            maximum_conditional_log / math.log(2.0)
        ),
    }


def saddle_gap_log_probability(
    visits: tuple[int, ...],
    distance: int,
    active_output_weight: int,
) -> float:
    """Lugannani--Rice approximation for the Dirichlet spacing tail."""
    active_packets = sum(visits)
    inactive_packets = PACKET_POSITIONS - active_packets
    threshold = (distance - active_output_weight) / inactive_packets
    if threshold <= 0.0:
        return -math.inf
    if threshold >= PACKET_WIDTH:
        return 0.0

    alpha = np.asarray(visits, dtype=np.float64)
    alpha[0] += 1.0  # the gap before the first active packet
    coefficients = np.arange(PACKET_WIDTH + 1, dtype=np.float64) - threshold

    def derivative(tilt: float) -> float:
        return float(
            np.sum(alpha * coefficients / (1.0 - tilt * coefficients))
        )

    lower = -1.0 / threshold + 1e-12
    tilt = brentq(derivative, lower, -1e-15, xtol=1e-14, rtol=1e-14)
    denominators = 1.0 - tilt * coefficients
    cumulant = -float(np.sum(alpha * np.log(denominators)))
    second = float(
        np.sum(alpha * coefficients * coefficients / (denominators * denominators))
    )
    w = -math.sqrt(max(0.0, -2.0 * cumulant))
    u = tilt * math.sqrt(second)
    if not w < 0.0 or not u < 0.0:
        raise RuntimeError("invalid lower-tail saddle")
    log_normal = float(log_ndtr(w))
    mills_ratio = math.exp(-0.5 * w * w - 0.5 * math.log(2.0 * math.pi) - log_normal)
    correction = 1.0 + mills_ratio * (1.0 / w - 1.0 / u)
    if correction <= 0.0:
        return log_normal
    return log_normal + math.log(correction)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=20_260_825)
    parser.add_argument(
        "--gap-method",
        choices=("saddle", "exact"),
        default="saddle",
    )
    parser.add_argument(
        "--relative-distances",
        type=float,
        nargs="+",
        default=list(DEFAULT_RELATIVE_DISTANCES),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/"
            "receipts/target_occ3_242424_exact_gap_estimate.json"
        ),
    )
    args = parser.parse_args()
    if args.samples <= 1:
        raise ValueError("samples must exceed one")
    if any(not 0.0 < value < 0.25 for value in args.relative_distances):
        raise ValueError("relative distances must lie in (0, 0.25)")

    grouped = sample_visit_histograms(args.samples, args.seed)
    distances = sorted(
        {int(math.floor(value * BINARY_LENGTH)) for value in args.relative_distances}
    )
    counter = FlexibleGapCounter()
    rows = [
        estimate_distance(grouped, args.samples, distance, counter, args.gap_method)
        for distance in distances
    ]
    crossing = next(
        (row for row in rows if row["estimated_shell_expected_count_log2"] >= 0.0),
        None,
    )
    payload = {
        "schema": "riffle-target-242424-exact-gap-estimate-v1",
        "evidence_label": "MONTE_CARLO_ACTIVE_PATHS_WITH_EXACT_GAP_INTEGRATION",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "profile": [24, 24, 24],
        "exact_outer_words": EXACT_OUTER_WORDS,
        "packet_positions": PACKET_POSITIONS,
        "binary_length": BINARY_LENGTH,
        "samples": args.samples,
        "seed": args.seed,
        "gap_method": args.gap_method,
        "distinct_visit_histograms": len(grouped),
        "rows": rows,
        "first_grid_crossing": crossing,
        "gap_cache_entries": len(counter._q_cache),
        "gap_interpolation_cache_entries": len(counter._interpolation_cache),
        "scope": (
            "The outer shell count is exact. Monte Carlo samples the "
            "independent within-block packet law and active-packet order. "
            + (
                "Each conditional gap probability is exact. "
                if args.gap_method == "exact"
                else "Conditional gap probabilities use the continuous-spacing "
                "Lugannani--Rice saddle approximation. "
            )
            + "The result estimates only the (24,24,24) shell, not the "
            "complete code enumerator."
        ),
    }
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "samples": args.samples,
                "distinct_visit_histograms": len(grouped),
                "rows": [
                    {
                        "relative_weight": row["relative_binary_weight"],
                        "shell_log2": row["estimated_shell_expected_count_log2"],
                        "relative_standard_error": row[
                            "relative_standard_error"
                        ],
                    }
                    for row in rows
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
