#!/usr/bin/env python3
"""Exactly count the single-transposition basin of the support-33 cluster."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

from probe_riffle_dp_g4_support33_cluster import contribution_words, word_for_assignment


ROOT = Path(__file__).resolve().parents[1]
CLUSTER_RECEIPT = ROOT / "explorations" / "riffle_dp_g4_g2_support33_cluster_family.json"
DEFAULT_OUTPUT = ROOT / "explorations" / "riffle_dp_g4_g2_support33_cluster_basin.json"
SAMPLE_SPACE = 524_352
DISTANCE_THRESHOLD = 188_766
CLUSTER_NODES = 3
CELL_COUNT = 48


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assignment_from_words(node_words: list[int]) -> list[int]:
    return [
        (node_words[node] >> (4 * slot)) & 0xF
        for node in range(CLUSTER_NODES)
        for slot in range(16)
    ]


def terminal_translation_count(word: np.ndarray) -> int:
    prefix = np.cumsum(np.bitwise_count(word), dtype=np.uint32)
    last_good = int(np.searchsorted(prefix, DISTANCE_THRESHOLD, side="right"))
    return max(0, last_good - CLUSTER_NODES + 1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=CLUSTER_RECEIPT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    receipt = json.loads(args.receipt.read_text())
    maximum_length = receipt["window_length"]
    node_words = [int(value, 16) for value in receipt["node_words_hex"]]
    assignment = assignment_from_words(node_words)
    if Counter(value for value in assignment if value) != Counter(
        {int(value): count for value, count in receipt["packet_value_histogram"].items()}
    ):
        raise RuntimeError("cluster basin: receipt assignment mismatch")
    if any(all(assignment[16 * node + slot] == 0 for slot in range(16)) for node in range(3)):
        raise RuntimeError("cluster basin: an edge node is empty")

    words = contribution_words(maximum_length)
    base_word = word_for_assignment(words, assignment)
    base_weight = int(np.bitwise_count(base_word).sum())
    if base_weight != receipt["weight"]:
        raise RuntimeError("cluster basin: base replay mismatch")

    assignment_count = 1
    total_arrangements = terminal_translation_count(base_word)
    minimum_translations = total_arrangements
    maximum_translations = total_arrangements
    subthreshold_at_maximum_length = int(base_weight <= DISTANCE_THRESHOLD)
    swap_rows = []
    for left in range(CELL_COUNT):
        for right in range(left + 1, CELL_COUNT):
            left_value = assignment[left]
            right_value = assignment[right]
            if left_value == right_value:
                continue
            delta_value = left_value ^ right_value
            candidate = base_word ^ words[left, delta_value] ^ words[right, delta_value]
            translations = terminal_translation_count(candidate)
            candidate_weight = int(np.bitwise_count(candidate).sum())
            assignment_count += 1
            total_arrangements += translations
            minimum_translations = min(minimum_translations, translations)
            maximum_translations = max(maximum_translations, translations)
            subthreshold_at_maximum_length += int(candidate_weight <= DISTANCE_THRESHOLD)
            swap_rows.append((candidate_weight, left, right, translations))

    expected_assignment_count = 1 + sum(
        1
        for left in range(CELL_COUNT)
        for right in range(left + 1, CELL_COUNT)
        if assignment[left] != assignment[right]
    )
    if assignment_count != expected_assignment_count:
        raise RuntimeError("cluster basin: assignment count mismatch")
    if minimum_translations <= 0:
        raise RuntimeError("cluster basin: terminal translations unexpectedly empty")

    value_counts = Counter(value for value in assignment if value)
    label_multiplicity = math.prod(math.factorial(count) for count in value_counts.values())
    falling = math.prod(range(SAMPLE_SPACE - 32, SAMPLE_SPACE + 1))
    single_family = Fraction(total_arrangements * label_multiplicity, falling)
    family_count = len(receipt["compatible_outer_family_ids"])
    aggregate = family_count * single_family
    log2_single = math.log2(single_family.numerator) - math.log2(single_family.denominator)
    log2_aggregate = math.log2(aggregate.numerator) - math.log2(aggregate.denominator)
    best_swap = min(swap_rows)
    payload = {
        "schema": "riffle-dp-g4-g2-support33-cluster-basin-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "EXACT",
        "cluster_receipt_sha256": digest(args.receipt),
        "source_sha256": digest(Path(__file__).resolve()),
        "distance_threshold": DISTANCE_THRESHOLD,
        "maximum_window_length": maximum_length,
        "basin": "The authenticated assignment and every distinct swap of two cells with different values.",
        "assignment_count": assignment_count,
        "subthreshold_assignment_count_at_maximum_length": subthreshold_at_maximum_length,
        "minimum_terminal_translations_per_assignment": minimum_translations,
        "maximum_terminal_translations_per_assignment": maximum_translations,
        "total_disjoint_physical_arrangements": total_arrangements,
        "best_swap": {
            "weight_at_maximum_length": best_swap[0],
            "left_cell": best_swap[1],
            "right_cell": best_swap[2],
            "terminal_translation_count": best_swap[3],
        },
        "single_family_event_probability": {
            "numerator": str(single_family.numerator),
            "denominator": str(single_family.denominator),
            "log2": log2_single,
        },
        "compatible_outer_family_count": family_count,
        "aggregate_first_moment_contribution": {
            "numerator": str(aggregate.numerator),
            "denominator": str(aggregate.denominator),
            "log2": log2_aggregate,
        },
        "scope_limitation": (
            "This exact basin contains only assignments within one transposition "
            "of the authenticated cluster. Its contribution remains a lower bound."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"assignments={assignment_count}")
    print(f"subthreshold_at_length_{maximum_length}={subthreshold_at_maximum_length}")
    print(f"total_disjoint_arrangements={total_arrangements}")
    print(f"best_swap_weight={best_swap[0]}")
    print(f"aggregate_log2_first_moment_contribution={log2_aggregate:.12f}")
    print(f"output={args.output}")
    print("status=EXACT_SUPPORT33_CLUSTER_SINGLE_SWAP_BASIN")


if __name__ == "__main__":
    main()
