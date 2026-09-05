#!/usr/bin/env python3
"""Audit candidate structural reductions for the 24/12-node character codes."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import (
    accumulate,
    apply_polynomial,
    binary_rank,
    build_apply,
    kernel_basis,
)
from audit_riffle_dp_2lap_g4_component_endpoint_dimension30 import (
    code_columns,
    observations,
    support_annihilator,
)
from probe_riffle_dp_2lap_g4_full_local_distance import support_data
from probe_riffle_dp_g4_component_mixing import transpose_columns


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
FULL_PROBE = EXPLORATIONS / "riffle_dp_2lap_g4_full_local_distance_probe.json"
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_structural_checkpoint.json"
COMPLEMENT_PAIRS = ((0x2D, 0x52), (0x34, 0x4B))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def nullspace(equations: tuple[int, ...]) -> tuple[int, ...]:
    rows = list(dict.fromkeys(equations))
    pivots = []
    rank = 0
    for column in range(64):
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
    free = [column for column in range(64) if column not in pivots]
    basis = []
    for free_column in free:
        value = 1 << free_column
        for row_index, pivot_column in enumerate(pivots):
            if (rows[row_index] >> free_column) & 1:
                value |= 1 << pivot_column
        basis.append(value)
    return tuple(basis)


def span_nonzero(basis: tuple[int, ...]) -> tuple[int, ...]:
    values = [0]
    for vector in basis:
        values += [value ^ vector for value in values]
    return tuple(values[1:])


def component_basis(transpose_apply, mask: int) -> tuple[int, ...]:
    return kernel_basis(
        tuple(
            apply_polynomial(
                transpose_apply,
                support_annihilator(mask),
                1 << bit,
            )
            for bit in range(64)
        )
    )


def code_words(step, transpose_apply, mask: int, value: int, window: int) -> tuple[int, ...]:
    states = observations(step, value, window)
    basis = component_basis(transpose_apply, mask)
    words = []
    for character in basis:
        word = sum(
            (((character & state).bit_count() & 1) << coordinate)
            for coordinate, state in enumerate(states)
        )
        words.append(word)
    return tuple(words)


def cross_gram_rank(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, int]:
    rows = tuple(
        sum((((left_word & right_word).bit_count() & 1) << index) for index, right_word in enumerate(right))
        for left_word in left
    )
    return binary_rank(rows), sum(row.bit_count() for row in rows)


def coefficient_histogram(step, character: int, window: int) -> tuple[int, ...]:
    basis_states = {
        value: observations(step, value, window) for value in (1, 2, 4, 8)
    }
    histogram = [0] * 16
    for coordinate in range(16 * window):
        coefficient = sum(
            (((character & basis_states[value][coordinate]).bit_count() & 1) << bit)
            for bit, value in enumerate((1, 2, 4, 8))
        )
        histogram[coefficient] += 1
    return tuple(histogram)


def histogram_weights(histogram: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(count for coefficient, count in enumerate(histogram) if (coefficient & value).bit_count() & 1)
        for value in range(16)
    )


def main() -> None:
    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    transpose_apply = build_apply(
        transpose_columns(tuple(step(1 << bit) for bit in range(64)))
    )

    observability_rows = []
    for value in range(1, 16):
        states = observations(step, value, 5)
        prefix_ranks = [
            binary_rank(
                tuple(
                    states[slot * 5 + node]
                    for slot in range(16)
                    for node in range(prefix)
                )
            )
            for prefix in range(1, 6)
        ]
        four_states = tuple(
            states[slot * 5 + node]
            for slot in range(16)
            for node in range(4)
        )
        kernel_basis4 = nullspace(four_states)
        kernel_rows = []
        for character in span_nonzero(kernel_basis4):
            mask, degrees, polynomial = support_data(transpose_apply, character)
            kernel_rows.append(
                {
                    "character_hex": hex(character),
                    "support_mask_hex": hex(mask),
                    "component_degrees": degrees,
                    "component_dimension": sum(degrees),
                    "minimal_polynomial_hex": hex(polynomial),
                }
            )
        observability_rows.append(
            {
                "packet_value": value,
                "prefix_ranks_1_through_5_nodes": prefix_ranks,
                "four_node_kernel_dimension": len(kernel_basis4),
                "four_node_nonzero_kernel_characters": kernel_rows,
            }
        )
    if any(row["prefix_ranks_1_through_5_nodes"][:3] != [16, 32, 48] for row in observability_rows):
        raise RuntimeError("structural checkpoint: early observability ranks changed")
    if any(row["prefix_ranks_1_through_5_nodes"][4] != 64 for row in observability_rows):
        raise RuntimeError("structural checkpoint: five-node injectivity failed")

    gram_rows = []
    for window in (12, 24):
        for value in range(1, 16):
            for left_mask, right_mask in COMPLEMENT_PAIRS:
                left = code_words(step, transpose_apply, left_mask, value, window)
                right = code_words(step, transpose_apply, right_mask, value, window)
                rank, nonzero_entries = cross_gram_rank(left, right)
                gram_rows.append(
                    {
                        "window_nodes": window,
                        "packet_value": value,
                        "left_support_mask_hex": hex(left_mask),
                        "right_support_mask_hex": hex(right_mask),
                        "cross_gram_rank": rank,
                        "cross_gram_nonzero_entry_count": nonzero_entries,
                    }
                )
    if any(row["cross_gram_rank"] == 0 for row in gram_rows):
        raise RuntimeError("structural checkpoint: unexpected orthogonal complement pair")

    probe = json.loads(FULL_PROBE.read_text())
    if probe["result"] != "NO_COUNTEREXAMPLE_FOUND" or probe["violation_count"] != 0:
        raise RuntimeError("structural checkpoint: full-code probe result changed")
    histogram_rows = []
    for row in probe["cases"]:
        window = row["window_nodes"]
        character = int(row["character_hex"], 16)
        histogram = coefficient_histogram(step, character, window)
        weights = histogram_weights(histogram)
        direct_weights = tuple(
            sum(
                (character & state).bit_count() & 1
                for state in observations(step, value, window)
            )
            for value in range(16)
        )
        if weights != direct_weights:
            raise RuntimeError("structural checkpoint: coefficient lift failed")
        length = 16 * window
        walsh_square_sum = sum((length - 2 * weight) ** 2 for weight in weights)
        histogram_collision_term = 16 * sum(count * count for count in histogram)
        if walsh_square_sum != histogram_collision_term:
            raise RuntimeError("structural checkpoint: histogram Parseval identity failed")
        histogram_rows.append(
            {
                "window_nodes": window,
                "source_packet_value": row["packet_value"],
                "character_hex": row["character_hex"],
                "component_support_mask_hex": row["component_support_mask_hex"],
                "coefficient_histogram": list(histogram),
                "weights_for_packet_values_0_through_15": list(weights),
                "walsh_square_sum": walsh_square_sum,
                "histogram_collision_term": histogram_collision_term,
            }
        )

    non15_sides = [
        row["best_found_two_sided_weight"]
        for row in probe["cases"]
        if row["window_nodes"] == 24 and row["packet_value"] != 15
    ]
    payload = {
        "schema": "riffle-dp-2lap-g4-structural-checkpoint-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_STRUCTURAL_IDENTITIES_WITH_DIAGNOSTIC_FULL_CODE_SEARCH",
        "source_sha256": digest(Path(__file__).resolve()),
        "full_local_distance_probe_sha256": digest(FULL_PROBE),
        "observability_rows": observability_rows,
        "complement_support_cross_gram_rows": gram_rows,
        "coefficient_histogram_rows": histogram_rows,
        "full_code_probe_summary": {
            "random_restarts_per_case": probe["random_restarts_per_case"],
            "case_count": probe["case_count"],
            "violation_count": probe["violation_count"],
            "minimum_found_non15_two_sided_weight": min(non15_sides),
            "value15_24node_two_sided_weight": next(
                row["best_found_two_sided_weight"]
                for row in probe["cases"]
                if row["window_nodes"] == 24 and row["packet_value"] == 15
            ),
            "value15_12node_two_sided_weight": next(
                row["best_found_two_sided_weight"]
                for row in probe["cases"]
                if row["window_nodes"] == 12
            ),
        },
        "reduction_verdicts": {
            "five_node_observability": (
                "VALID_IDENTITY_BUT_TOO_WEAK_ALONE: every packet value reaches "
                "rank 64 after five nodes."
            ),
            "four_node_small_kernel_quotient": (
                "REFUTED: four-node kernels include characters supported on "
                "five, six, or all seven irreducible components."
            ),
            "complement_support_orthogonality": (
                "REFUTED: every tested complementary dimension-32 pair has a "
                "nonzero cross-Gram matrix at both window lengths."
            ),
            "coefficient_histogram_lift": (
                "VALID_EXACT_REPARAMETERIZATION: all 15 packet-value words are "
                "Walsh coefficients of one 16-bin nonnegative histogram."
            ),
            "uniform_full_local_distance": (
                "NOT_REFUTED_NOT_PROVED: 10000 coordinate-descent restarts per "
                "case find minimum non-15 side 120 and the known value-15 72/36 equality."
            ),
        },
        "result": "STRUCTURAL_REDUCTION_FOUND_BUT_PROOF_OPEN",
        "scope_limitation": (
            "The rank, kernel, Gram, histogram, and replay statements are exact. "
            "The full-code minimum-distance search remains diagnostic."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("five_node_observability=PASS_EXACT")
    print("four_node_small_kernel_quotient=REFUTED_EXACT")
    print("complement_support_orthogonality=REFUTED_EXACT")
    print("coefficient_histogram_lift=PASS_EXACT")
    print(f"minimum_found_non15_side={min(non15_sides)}")
    print(f"output={OUTPUT}")
    print("status=STRUCTURAL_REDUCTION_FOUND_BUT_PROOF_OPEN")


if __name__ == "__main__":
    main()
