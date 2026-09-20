#!/usr/bin/env python3
"""Exact support-33 lower families from same-node packet collisions."""

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
TURNOFF_RECEIPT = EXPLORATIONS / "riffle_dp_g4_g2_multi_packet_node_turnoffs.json"
OUTPUT = EXPLORATIONS / "riffle_dp_g4_g2_collision_turnoff_lower_family.json"
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


def labeled_count(counts: Counter[int], values: tuple[int, ...]) -> int:
    remaining = counts.copy()
    result = 1
    for value in values:
        result *= remaining[value]
        remaining[value] -= 1
    return result


def main() -> None:
    outer = json.loads(OUTER_RECEIPTS.read_text())
    turnoff = json.loads(TURNOFF_RECEIPT.read_text())
    witnesses = turnoff["below_threshold_witnesses"]
    if len(outer["receipts"]) != 26 or len(witnesses) != 539:
        raise RuntimeError("collision family: unexpected input cardinality")
    if turnoff["maximum_reset_node"] != INNER_NODES - 1:
        raise RuntimeError("collision family: search does not cover the chain")

    total = Fraction()
    members = []
    for receipt in outer["receipts"]:
        active_values = []
        for local_word in receipt["local_words"]:
            active_values.extend(packet_values(local_word["codeword_hex"]))
        if len(active_values) != 33:
            raise RuntimeError("collision family: outer support changed")
        counts = Counter(active_values)
        category_totals: dict[int, Fraction] = defaultdict(Fraction)
        category_matches: dict[int, int] = defaultdict(int)
        for witness in witnesses:
            values = tuple(
                nibble["value"]
                for group in (witness["initial_nibbles"], witness["reset_nibbles"])
                for nibble in group
            )
            outside_packets = len(values)
            choices = labeled_count(counts, values)
            if not choices:
                continue
            emitted_weight = witness["emitted_weight_before_turnoff"]
            terminal_nodes = (DISTANCE_THRESHOLD - emitted_weight) // 64
            prefix_nodes = INNER_NODES - terminal_nodes
            start_count = prefix_nodes - witness["reset_node"]
            if start_count <= 0:
                continue
            probability = Fraction(
                start_count
                * choices
                * falling(16 * terminal_nodes, 33 - outside_packets),
                falling(PACKET_COUNT, 33),
            )
            category_totals[outside_packets] += probability
            category_matches[outside_packets] += 1
        if not category_totals:
            continue
        selected_packets, selected_probability = max(
            category_totals.items(), key=lambda row: row[1]
        )
        total += selected_probability
        members.append(
            {
                "outer_family_id": receipt["family_id"],
                "active_packet_value_counts": {
                    str(value): count for value, count in sorted(counts.items())
                },
                "categories": [
                    {
                        "outside_packets": outside_packets,
                        "matched_witnesses": category_matches[outside_packets],
                        "probability_numerator": str(probability.numerator),
                        "probability_denominator": str(probability.denominator),
                        "log2": math.log2(probability.numerator)
                        - math.log2(probability.denominator),
                    }
                    for outside_packets, probability in sorted(category_totals.items())
                ],
                "selected_outside_packets": selected_packets,
                "selected_probability_numerator": str(selected_probability.numerator),
                "selected_probability_denominator": str(selected_probability.denominator),
            }
        )

    payload = {
        "schema": "riffle-dp-g4-g2-collision-turnoff-lower-family-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "EXACT",
        "outer_receipt_sha256": digest(OUTER_RECEIPTS),
        "turnoff_receipt_sha256": digest(TURNOFF_RECEIPT),
        "outer_member_count": len(members),
        "aggregate_probability": {
            "numerator": str(total.numerator),
            "denominator": str(total.denominator),
            "log2": math.log2(total.numerator) - math.log2(total.denominator),
            "exceeds_2^-40": total > Fraction(1, 1 << 40),
        },
        "members": members,
        "scope": (
            "For each outer word, select its strongest disjoint category with "
            "two, three, or four packets outside the terminal suffix. The outside "
            "packets form one exact full-chain turnoff episode."
        ),
        "disjointness": (
            "Within a fixed outside-packet category, a placement uniquely fixes "
            "the episode start, reset offset, nibble slots, and labeled packets."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"outer_words_with_collision_turnoff={len(members)}")
    print(f"aggregate_log2={payload['aggregate_probability']['log2']:.12f}")
    print(f"exceeds_2^-40={payload['aggregate_probability']['exceeds_2^-40']}")
    print(f"output={OUTPUT}")
    print("status=EXACT_G2_COLLISION_TURNOFF_LOWER_FAMILY")


if __name__ == "__main__":
    main()
