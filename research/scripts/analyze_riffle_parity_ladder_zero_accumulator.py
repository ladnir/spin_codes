#!/usr/bin/env python3
"""Compose the zero-parity weight-22 shell with the packet accumulator.

The calculation is exact.  It enumerates active four-bit state paths by their
state-weight visit counts.  A root-of-unity quasipolynomial evaluates all
placements of the inactive packets without iterating over the global length.
"""

from __future__ import annotations

import argparse
import functools
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = (
    REPOSITORY_ROOT
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "receipts"
    / "goal17_zero_parity_accumulator_exact.json"
)

PACKET_WIDTH = 4
PACKETS_PER_BLOCK = 32
BCH_LENGTH = 128
BCH_WEIGHT = 22
BCH_WEIGHT22_WORDS = 243_840
DATA_BLOCKS = 1 << 14
PACKET_POSITIONS = DATA_BLOCKS * PACKETS_PER_BLOCK
MINIMUM_SHELL_WORDS = DATA_BLOCKS * BCH_WEIGHT22_WORDS
GAP_PERIOD = math.lcm(*range(1, PACKET_WIDTH + 1))


def transition_multiplicity(old_weight: int, new_weight: int, packet_weight: int) -> int:
    numerator = old_weight + new_weight - packet_weight
    if numerator & 1:
        return 0
    overlap = numerator // 2
    added = new_weight - overlap
    if not 0 <= overlap <= old_weight or not 0 <= added <= PACKET_WIDTH - old_weight:
        return 0
    return math.comb(old_weight, overlap) * math.comb(
        PACKET_WIDTH - old_weight, added
    )


def compositions(total: int, parts: int, prefix: tuple[int, ...] = ()):
    if parts == 1:
        yield prefix + (total,)
        return
    for first in range(total + 1):
        yield from compositions(total - first, parts - 1, prefix + (first,))


def active_visit_counts() -> list[tuple[int, dict[tuple[int, ...], int]]]:
    """Count ordered nonzero packet paths of total binary weight 22."""
    current: dict[tuple[int, int, tuple[int, ...]], int] = {
        (0, 0, (0, 0, 0, 0, 0)): 1
    }
    result = []
    for active_packets in range(1, BCH_WEIGHT + 1):
        following: defaultdict[tuple[int, int, tuple[int, ...]], int] = defaultdict(int)
        for (input_weight, old_weight, visits), ways in current.items():
            for packet_weight in range(1, PACKET_WIDTH + 1):
                if input_weight + packet_weight > BCH_WEIGHT:
                    continue
                for new_weight in range(PACKET_WIDTH + 1):
                    multiplicity = transition_multiplicity(
                        old_weight, new_weight, packet_weight
                    )
                    if not multiplicity:
                        continue
                    next_visits = list(visits)
                    next_visits[new_weight] += 1
                    following[
                        (
                            input_weight + packet_weight,
                            new_weight,
                            tuple(next_visits),
                        )
                    ] += ways * multiplicity
        current = following
        grouped: defaultdict[tuple[int, ...], int] = defaultdict(int)
        for (input_weight, _state_weight, visits), ways in current.items():
            if input_weight == BCH_WEIGHT:
                grouped[visits] += ways
        if grouped:
            result.append((active_packets, dict(grouped)))
    return result


def log2_integer(value: int) -> float:
    if value <= 0:
        return -math.inf
    top_bit = value.bit_length() - 1
    return top_bit + math.log2(value / (1 << top_bit))


def log2_fraction(value: Fraction) -> float:
    if value <= 0:
        return -math.inf
    return log2_integer(value.numerator) - log2_integer(value.denominator)


def rational_payload(value: Fraction) -> dict[str, int | str | float]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": format(float(value), ".17g"),
        "log2": log2_fraction(value),
    }


class GapCounter:
    def __init__(self) -> None:
        self._q_cache: dict[tuple[tuple[int, ...], int], int] = {}
        self._derivative_compositions = {
            zero_visits: tuple(
                q
                for degree in range(zero_visits + 1)
                for q in compositions(degree, PACKET_WIDTH)
            )
            for zero_visits in range(BCH_WEIGHT // 2 + 1)
        }

    def positive_gap_count(self, exponents: tuple[int, ...], maximum_cost: int) -> int:
        """Return [z^<=R] product_b (1-z^b)^(-exponents[b])."""
        if maximum_cost < 0:
            return 0
        if not any(exponents):
            return 1
        key = (exponents, maximum_cost)
        cached = self._q_cache.get(key)
        if cached is not None:
            return cached

        degree = sum(exponents)
        residue = maximum_cost % GAP_PERIOD
        quotient = maximum_cost // GAP_PERIOD
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
        result = sum(
            difference * math.comb(quotient, order)
            for order, difference in enumerate(first_differences)
        )
        self._q_cache[key] = result
        return result

    def all_gap_count(
        self,
        visits: tuple[int, ...],
        inactive_packets: int,
        maximum_extra_weight: int,
    ) -> int:
        """Count all gap placements below one extra-output threshold."""
        zero_visits = visits[0]
        positive_visits = visits[1:]
        result = 0
        for derivative_counts in self._derivative_compositions[zero_visits]:
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
                exponents, maximum_extra_weight - shift
            )
        return result


def exact_shell_probability(
    distance: int,
    active_groups: list[tuple[int, dict[tuple[int, ...], int]]],
) -> tuple[Fraction, list[dict[str, object]], int]:
    counter = GapCounter()
    probability = Fraction(0)
    contributions = []
    slice_size = math.comb(BCH_LENGTH, BCH_WEIGHT)
    for active_packets, groups in active_groups:
        inactive_packets = PACKET_POSITIONS - active_packets
        bad_sequences = 0
        for visits, active_paths in groups.items():
            active_output_weight = sum(
                state_weight * visits[state_weight]
                for state_weight in range(1, PACKET_WIDTH + 1)
            )
            gap_count = counter.all_gap_count(
                visits,
                inactive_packets,
                distance - active_output_weight,
            )
            bad_sequences += active_paths * gap_count
        contribution = Fraction(
            math.comb(PACKETS_PER_BLOCK, active_packets) * bad_sequences,
            slice_size * math.comb(PACKET_POSITIONS, active_packets),
        )
        probability += contribution
        contributions.append(
            {
                "active_packets": active_packets,
                "visit_histograms": len(groups),
                "bad_full_input_sequences": bad_sequences,
                "probability": rational_payload(contribution),
            }
        )
    return probability, contributions, len(counter._q_cache)


def direct_positive_gap_count(exponents: tuple[int, ...], maximum_cost: int) -> int:
    coefficients = [0] * (maximum_cost + 1)
    coefficients[0] = 1
    for weight, multiplicity in enumerate(exponents, 1):
        for _ in range(multiplicity):
            for index in range(weight, maximum_cost + 1):
                coefficients[index] += coefficients[index - weight]
    return sum(coefficients)


def audit_quasipolynomial() -> int:
    counter = GapCounter()
    comparisons = 0
    for degree in range(1, 7):
        for exponents in compositions(degree, PACKET_WIDTH):
            for maximum_cost in range(73):
                exact = direct_positive_gap_count(exponents, maximum_cost)
                observed = counter.positive_gap_count(exponents, maximum_cost)
                if observed != exact:
                    raise RuntimeError("gap quasipolynomial failed a direct coefficient check")
                comparisons += 1
    return comparisons


def audit_all_gap_formula() -> int:
    counter = GapCounter()
    comparisons = 0
    for active_packets in range(1, 6):
        for visits in compositions(active_packets, PACKET_WIDTH + 1):
            inactive_packets = 8
            for maximum_extra_weight in range(inactive_packets):
                direct = 0
                state_weights = tuple(
                    weight
                    for weight, multiplicity in enumerate(visits)
                    for _ in range(multiplicity)
                )
                for gaps in compositions(inactive_packets, active_packets + 1):
                    cost = sum(
                        state_weights[index] * gaps[index + 1]
                        for index in range(active_packets)
                    )
                    direct += int(cost <= maximum_extra_weight)
                observed = counter.all_gap_count(
                    visits, inactive_packets, maximum_extra_weight
                )
                if observed != direct:
                    raise RuntimeError("all-gap formula failed exhaustive enumeration")
                comparisons += 1
    return comparisons


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--distances",
        type=int,
        nargs="+",
        default=[76_000, 77_000, 188_743],
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if any(distance < BCH_WEIGHT // 2 for distance in args.distances):
        raise ValueError("distance lies below the deterministic shell floor")
    if any(distance >= PACKET_POSITIONS for distance in args.distances):
        raise ValueError("exact gap chamber requires distance < packet positions")

    quasipolynomial_comparisons = audit_quasipolynomial()
    all_gap_comparisons = audit_all_gap_formula()
    active_groups = active_visit_counts()
    rows = []
    for distance in args.distances:
        probability, contributions, gap_cache_entries = exact_shell_probability(
            distance, active_groups
        )
        expected_bad_words = probability * MINIMUM_SHELL_WORDS
        rows.append(
            {
                "bad_output_weight_inclusive": distance,
                "relative_binary_output_weight": distance
                / (PACKET_WIDTH * PACKET_POSITIONS),
                "bad_probability_per_minimum_shell_word": rational_payload(probability),
                "expected_bad_minimum_shell_words": rational_payload(
                    expected_bad_words
                ),
                "first_moment_margin_bits": -log2_fraction(expected_bad_words),
                "passes_minimum_shell_first_moment": expected_bad_words < 1,
                "active_support_contributions": contributions,
                "gap_quasipolynomial_cache_entries": gap_cache_entries,
            }
        )

    payload = {
        "schema": "riffle-parity-ladder-zero-accumulator-exact-v1",
        "evidence_label": "EXACT_MINIMUM_SHELL_EXPECTED_ENUMERATOR",
        "stage": "zero parity",
        "parameters": {
            "data_blocks": DATA_BLOCKS,
            "packet_positions": PACKET_POSITIONS,
            "binary_length": PACKET_WIDTH * PACKET_POSITIONS,
            "minimum_bch_weight": BCH_WEIGHT,
            "minimum_shell_outer_words": MINIMUM_SHELL_WORDS,
            "minimum_shell_outer_words_log2": math.log2(MINIMUM_SHELL_WORDS),
        },
        "distance_rows": rows,
        "validation": {
            "direct_positive_gap_coefficient_comparisons": quasipolynomial_comparisons,
            "exhaustive_all_gap_comparisons": all_gap_comparisons,
            "quasipolynomial_period": GAP_PERIOD,
            "all_checks": "PASS",
        },
        "scope": (
            "The result exactly composes every weight-22 occupation-one word "
            "with the accumulator. It is a lower contribution to the complete "
            "zero-parity expected enumerator because higher BCH weights and "
            "higher occupations are not included. A positive shell margin does "
            "not by itself prove a distance for the complete code."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "distance_rows": [
                    {
                        "distance": row["bad_output_weight_inclusive"],
                        "margin_bits": row["first_moment_margin_bits"],
                    }
                    for row in rows
                ],
                "status": "PASS",
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
