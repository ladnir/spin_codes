#!/usr/bin/env python3
"""Audit the ordered-state content lost by the coefficient histogram."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import accumulate, binary_rank, build_apply
from audit_riffle_dp_2lap_g4_component_endpoint_dimension30 import observations
from probe_riffle_dp_g4_component_mixing import transpose_columns


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "explorations" / "riffle_dp_2lap_g4_histogram_transition.json"
WINDOW = 24
SLOTS = 16
FIVE_NODE_LENGTH = 80
LOW_FIVE_NODE_THRESHOLD = 24

# These states have the same current nibble histogram and different successor
# histograms.  The audit verifies the witness directly; it does not rely on a
# probabilistic search.
HISTOGRAM_COLLISION_LEFT = 0x7373D52AF3615F19
HISTOGRAM_COLLISION_RIGHT = 0x953F133D5F76721A


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def nibble_histogram(state: int) -> tuple[int, ...]:
    result = [0] * 16
    for slot in range(SLOTS):
        result[(state >> (4 * slot)) & 15] += 1
    return tuple(result)


def value_weight(histogram: tuple[int, ...], value: int) -> int:
    return sum(
        count
        for coefficient, count in enumerate(histogram)
        if (coefficient & value).bit_count() & 1
    )


def nullspace(equations: tuple[int, ...], width: int) -> tuple[int, ...]:
    rows = list(equations)
    pivots = []
    rank = 0
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(rows)) if (rows[index] >> column) & 1),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for index in range(len(rows)):
            if index != rank and ((rows[index] >> column) & 1):
                rows[index] ^= rows[rank]
        pivots.append(column)
        rank += 1
    free_columns = tuple(column for column in range(width) if column not in pivots)
    basis = []
    for free_column in free_columns:
        value = 1 << free_column
        for row_index, pivot_column in enumerate(pivots):
            if (rows[row_index] >> free_column) & 1:
                value |= 1 << pivot_column
        basis.append(value)
    return tuple(basis)


def dual_weight_spectrum(columns: tuple[int, ...]) -> Counter[int]:
    equations = tuple(
        sum((((column >> bit) & 1) << index) for index, column in enumerate(columns))
        for bit in range(64)
    )
    basis = nullspace(equations, len(columns))
    if len(basis) != len(columns) - binary_rank(columns):
        raise RuntimeError("histogram transition: dual dimension mismatch")
    words = [0]
    for vector in basis:
        words += [word ^ vector for word in words]
    return Counter(word.bit_count() for word in words)


def krawtchouk(length: int, output_weight: int, input_weight: int) -> int:
    lower = max(0, output_weight - (length - input_weight))
    upper = min(output_weight, input_weight)
    return sum(
        (-1) ** overlap
        * math.comb(input_weight, overlap)
        * math.comb(length - input_weight, output_weight - overlap)
        for overlap in range(lower, upper + 1)
    )


def primal_weight_spectrum(
    length: int, dual_spectrum: Counter[int]
) -> tuple[int, ...]:
    dual_size = sum(dual_spectrum.values())
    result = []
    for weight in range(length + 1):
        numerator = sum(
            count * krawtchouk(length, weight, dual_weight)
            for dual_weight, count in dual_spectrum.items()
        )
        if numerator % dual_size:
            raise RuntimeError("histogram transition: nonintegral MacWilliams coefficient")
        result.append(numerator // dual_size)
    return tuple(result)


def coefficient_from_observations(
    basis_observations: dict[int, tuple[int, ...]], character: int, coordinate: int
) -> int:
    return sum(
        (((character & basis_observations[value][coordinate]).bit_count() & 1) << bit)
        for bit, value in enumerate((1, 2, 4, 8))
    )


def main() -> None:
    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    transpose_step = build_apply(
        transpose_columns(tuple(step(1 << bit) for bit in range(64)))
    )

    # Verify on a basis of characters that the coefficient word at node t and
    # slot s is nibble s of (T^(t+1) character).  Linearity then proves the
    # identity for all 64-bit characters.
    basis_observations = {
        value: observations(step, value, WINDOW) for value in (1, 2, 4, 8)
    }
    for character_bit in range(64):
        character = 1 << character_bit
        transpose_state = character
        for node in range(WINDOW):
            transpose_state = transpose_step(transpose_state)
            for slot in range(SLOTS):
                coordinate = slot * WINDOW + node
                coefficient = coefficient_from_observations(
                    basis_observations, character, coordinate
                )
                expected = (transpose_state >> (4 * slot)) & 15
                if coefficient != expected:
                    raise RuntimeError("histogram transition: nibble identity failed")

    left_histogram = nibble_histogram(HISTOGRAM_COLLISION_LEFT)
    right_histogram = nibble_histogram(HISTOGRAM_COLLISION_RIGHT)
    left_successor_histogram = nibble_histogram(
        transpose_step(HISTOGRAM_COLLISION_LEFT)
    )
    right_successor_histogram = nibble_histogram(
        transpose_step(HISTOGRAM_COLLISION_RIGHT)
    )
    if left_histogram != right_histogram:
        raise RuntimeError("histogram transition: collision witness changed")
    if left_successor_histogram == right_successor_histogram:
        raise RuntimeError("histogram transition: successor witness changed")

    five_node_rows = []
    for packet_value in range(1, 16):
        prefix_ranks = tuple(
            binary_rank(observations(step, packet_value, nodes))
            for nodes in range(1, 6)
        )
        columns = observations(step, packet_value, 5)
        dual_spectrum = dual_weight_spectrum(columns)
        spectrum = primal_weight_spectrum(FIVE_NODE_LENGTH, dual_spectrum)
        if spectrum[0] != 1 or any(count < 0 for count in spectrum):
            raise RuntimeError("histogram transition: invalid primal spectrum")
        if sum(spectrum) != 1 << prefix_ranks[-1]:
            raise RuntimeError("histogram transition: primal spectrum mass mismatch")
        minimum_weight = next(weight for weight in range(1, 81) if spectrum[weight])
        maximum_weight = max(weight for weight in range(81) if spectrum[weight])
        minimum_complement_weight = FIVE_NODE_LENGTH - maximum_weight
        two_sided_distance = min(minimum_weight, minimum_complement_weight)
        low_weight_scalar_sequence = tuple(
            minimum_weight if node in (4, 9, 14, 19) else 0
            for node in range(WINDOW)
        )
        low_complement_scalar_sequence = tuple(
            minimum_complement_weight if node in (4, 9, 14, 19) else 0
            for node in range(WINDOW)
        )
        if any(
            sum(low_weight_scalar_sequence[start : start + 5]) < minimum_weight
            for start in range(WINDOW - 4)
        ):
            raise RuntimeError("histogram transition: scalar weight witness failed")
        if any(
            sum(low_complement_scalar_sequence[start : start + 5])
            < minimum_complement_weight
            for start in range(WINDOW - 4)
        ):
            raise RuntimeError("histogram transition: scalar complement witness failed")
        low_candidate_count = sum(spectrum[: LOW_FIVE_NODE_THRESHOLD + 1]) + sum(
            spectrum[FIVE_NODE_LENGTH - LOW_FIVE_NODE_THRESHOLD :]
        )
        value15_candidate_count = sum(spectrum[:18]) + sum(spectrum[63:])
        five_node_rows.append(
            {
                "packet_value": packet_value,
                "prefix_ranks_1_through_5_nodes": list(prefix_ranks),
                "dual_dimension": int(math.log2(sum(dual_spectrum.values()))),
                "dual_minimum_weight": min(
                    weight for weight, count in dual_spectrum.items() if weight and count
                ),
                "minimum_weight": minimum_weight,
                "maximum_weight": maximum_weight,
                "minimum_complement_weight": minimum_complement_weight,
                "two_sided_distance": two_sided_distance,
                "scalar_24node_weight_bound": 4 * minimum_weight,
                "scalar_24node_complement_bound": 4 * minimum_complement_weight,
                "scalar_24node_two_sided_bound": 4 * two_sided_distance,
                "five_node_two_sided_count_at_most_24": low_candidate_count,
                "five_node_two_sided_count_at_most_24_log2": math.log2(
                    low_candidate_count
                ),
                "five_node_two_sided_count_at_most_17": value15_candidate_count,
                "five_node_two_sided_count_at_most_17_log2": math.log2(
                    value15_candidate_count
                ),
                "primal_counts_weights_0_through_24": list(
                    spectrum[: LOW_FIVE_NODE_THRESHOLD + 1]
                ),
                "primal_counts_weights_56_through_80": list(
                    spectrum[FIVE_NODE_LENGTH - LOW_FIVE_NODE_THRESHOLD :]
                ),
            }
        )

    payload = {
        "schema": "riffle-dp-2lap-g4-histogram-transition-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_LINEAR_IDENTITIES_AND_EXACT_MACWILLIAMS_SPECTRA",
        "source_sha256": digest(Path(__file__).resolve()),
        "coefficient_nibble_identity": {
            "result": "PASS_EXACT",
            "statement": (
                "At node t and slot s, the four-bit coefficient word equals "
                "nibble s of T^(t+1) chi. The 24-node histogram counts all "
                "nibbles in T chi through T^24 chi."
            ),
            "verification": (
                "Checked every slot and node on all 64 basis characters; "
                "linearity extends the identity to every character."
            ),
        },
        "histogram_transition_counterexample": {
            "result": "REFUTED_EXACT",
            "claim_refuted": (
                "The current 16-bin nibble histogram determines the next histogram."
            ),
            "left_state_hex": hex(HISTOGRAM_COLLISION_LEFT),
            "right_state_hex": hex(HISTOGRAM_COLLISION_RIGHT),
            "common_current_histogram": list(left_histogram),
            "left_successor_histogram": list(left_successor_histogram),
            "right_successor_histogram": list(right_successor_histogram),
            "left_successor_packet_weights_1_through_15": [
                value_weight(left_successor_histogram, value)
                for value in range(1, 16)
            ],
            "right_successor_packet_weights_1_through_15": [
                value_weight(right_successor_histogram, value)
                for value in range(1, 16)
            ],
        },
        "five_node_rows": five_node_rows,
        "summary": {
            "one_two_three_node_ranks": [16, 32, 48],
            "first_uniform_redundancy_node_count": 5,
            "five_node_code_parameters": "[80,64] for every packet value",
            "five_node_two_sided_distance_range": [
                min(row["two_sided_distance"] for row in five_node_rows),
                max(row["two_sided_distance"] for row in five_node_rows),
            ],
            "scalar_24node_two_sided_bound_range": [
                min(row["scalar_24node_two_sided_bound"] for row in five_node_rows),
                max(row["scalar_24node_two_sided_bound"] for row in five_node_rows),
            ],
            "five_node_low_candidate_log2_range_at_two_sided_weight_24": [
                min(
                    row["five_node_two_sided_count_at_most_24_log2"]
                    for row in five_node_rows
                ),
                max(
                    row["five_node_two_sided_count_at_most_24_log2"]
                    for row in five_node_rows
                ),
            ],
        },
        "result": "HISTOGRAM_ONLY_TRANSFER_REFUTED_AND_FIVE_NODE_RELAXATION_TOO_WEAK",
        "scope_limitation": (
            "The identities, counterexample, ranks, and spectra are exact. "
            "The audit does not rule out a stronger ordered-state quotient or "
            "a different algebraic proof of the 24-node distance."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("coefficient_nibble_identity=PASS_EXACT")
    print("histogram_markov_property=REFUTED_EXACT")
    print("five_node_spectra=PASS_EXACT")
    print(
        "five_node_two_sided_distance_range="
        f"{payload['summary']['five_node_two_sided_distance_range']}"
    )
    print(
        "scalar_24node_bound_range="
        f"{payload['summary']['scalar_24node_two_sided_bound_range']}"
    )
    print(f"output={OUTPUT}")
    print(f"status={payload['result']}")


if __name__ == "__main__":
    main()
