#!/usr/bin/env python3
"""Exact repeated two-packet-turnoff families for support-33 outer words."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
OUTER_RECEIPTS = EXPLORATIONS / "riffle_dp_g4_g1_support33_refutation.json"
ORBIT_RECEIPT = EXPLORATIONS / "riffle_dp_g4_g2_single_packet_orbits.json"
OUTPUT = EXPLORATIONS / "riffle_dp_g4_g2_multi_turnoff_lower_family.json"
PACKET_COUNT = 524_352
INNER_NODES = 32_772
DISTANCE_THRESHOLD = 188_766


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def falling(value: int, count: int) -> int:
    return math.prod(range(value - count + 1, value + 1))


def packet_values(codeword_hex: str) -> list[int]:
    word = int(codeword_hex, 16)
    return [
        (word >> (4 * packet)) & 0xF
        for packet in range(32)
        if ((word >> (4 * packet)) & 0xF) != 0
    ]


def sequence_counts(
    counts: Counter[int], transitions: list[dict]
) -> list[dict[int, int]]:
    values = tuple(sorted(counts))
    index = {value: position for position, value in enumerate(values)}
    initial = tuple(counts[value] for value in values)
    current: dict[tuple[tuple[int, ...], int], int] = {(initial, 0): 1}
    rows: list[dict[int, int]] = [defaultdict(int)]
    rows[0][0] = 1
    for episode_count in range(1, sum(counts.values()) // 2 + 1):
        following: dict[tuple[tuple[int, ...], int], int] = defaultdict(int)
        weight_counts: dict[int, int] = defaultdict(int)
        for (remaining_tuple, total_weight), sequence_count in current.items():
            for transition in transitions:
                first = transition["first_value"]
                second = transition["second_value"]
                if first not in index or second not in index:
                    continue
                remaining = list(remaining_tuple)
                first_count = remaining[index[first]]
                second_count = remaining[index[second]] - (first == second)
                choices = first_count * second_count
                if choices <= 0:
                    continue
                remaining[index[first]] -= 1
                remaining[index[second]] -= 1
                next_weight = (
                    total_weight + transition["emitted_weight_before_turnoff"]
                )
                ways = sequence_count * choices
                key = (tuple(remaining), next_weight)
                following[key] += ways
                weight_counts[next_weight] += ways
        if not following:
            break
        rows.append(weight_counts)
        current = following
    return rows


def family_probability(
    *,
    support: int,
    episode_count: int,
    episode_span: int,
    sequence_weight_counts: dict[int, int],
) -> Fraction:
    total = Fraction()
    for emitted_weight, sequence_count in sequence_weight_counts.items():
        terminal_nodes = (DISTANCE_THRESHOLD - emitted_weight) // 64
        if terminal_nodes < 0:
            continue
        prefix_nodes = INNER_NODES - terminal_nodes
        if prefix_nodes < episode_span * episode_count:
            continue
        start_sets = math.comb(
            prefix_nodes - (episode_span - 1) * episode_count,
            episode_count,
        )
        total += Fraction(
            start_sets
            * sequence_count
            * falling(16 * terminal_nodes, support - 2 * episode_count),
            falling(PACKET_COUNT, support),
        )
    return total


def main() -> None:
    outer = json.loads(OUTER_RECEIPTS.read_text())
    orbit = json.loads(ORBIT_RECEIPT.read_text())
    transitions = orbit["exact_one_packet_turnoffs_below_distance"]
    if len(outer["receipts"]) != 26 or len(transitions) != 6:
        raise RuntimeError("multi-turnoff family: unexpected input cardinality")
    zero_steps = {row["zero_steps_after_activation"] for row in transitions}
    if len(zero_steps) != 1:
        raise RuntimeError("multi-turnoff family: transition spans differ")
    episode_span = zero_steps.pop() + 2

    members = []
    total = Fraction()
    for receipt in outer["receipts"]:
        values = []
        for local_word in receipt["local_words"]:
            values.extend(packet_values(local_word["codeword_hex"]))
        if len(values) != 33:
            raise RuntimeError("multi-turnoff family: outer support changed")
        counts = Counter(values)
        episode_rows = sequence_counts(counts, transitions)
        candidates = []
        for episode_count in range(1, len(episode_rows)):
            probability = family_probability(
                support=33,
                episode_count=episode_count,
                episode_span=episode_span,
                sequence_weight_counts=episode_rows[episode_count],
            )
            if probability:
                candidates.append((probability, episode_count))
        if not candidates:
            continue
        best_probability, best_count = max(candidates)
        total += best_probability
        members.append(
            {
                "outer_family_id": receipt["family_id"],
                "active_packet_value_counts": {
                    str(value): count for value, count in sorted(counts.items())
                },
                "episode_probabilities": [
                    {
                        "episode_count": episode_count,
                        "log2": math.log2(probability.numerator)
                        - math.log2(probability.denominator),
                    }
                    for probability, episode_count in sorted(
                        candidates, key=lambda row: row[1]
                    )
                ],
                "selected_episode_count": best_count,
                "probability_numerator": str(best_probability.numerator),
                "probability_denominator": str(best_probability.denominator),
            }
        )

    payload = {
        "schema": "riffle-dp-g4-g2-multi-turnoff-lower-family-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "EXACT",
        "outer_receipt_sha256": digest(OUTER_RECEIPTS),
        "orbit_receipt_sha256": digest(ORBIT_RECEIPT),
        "outer_member_count": len(members),
        "episode_span_nodes": episode_span,
        "aggregate_probability": {
            "numerator": str(total.numerator),
            "denominator": str(total.denominator),
            "log2": math.log2(total.numerator) - math.log2(total.denominator),
            "exceeds_2^-40": total > Fraction(1, 1 << 40),
        },
        "members": members,
        "scope": (
            "For each included support-33 outer word, use the episode count with "
            "the largest exact probability. Each episode has the committed fixed "
            "span and resets the state. All other packets occupy a terminal suffix."
        ),
        "disjointness": (
            "Within one selected episode count, a placement uniquely determines "
            "its chronological episode starts, transition slots, and labeled packets."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"outer_words_with_repeated_turnoffs={len(members)}")
    print(f"aggregate_log2={payload['aggregate_probability']['log2']:.12f}")
    print(f"exceeds_2^-40={payload['aggregate_probability']['exceeds_2^-40']}")
    print(f"output={OUTPUT}")
    print("status=EXACT_G2_MULTI_TURNOFF_LOWER_FAMILY")


if __name__ == "__main__":
    main()
