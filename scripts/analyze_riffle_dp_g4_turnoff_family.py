#!/usr/bin/env python3
"""Exact two-packet turnoff family for Riffle DP g=4 G2.

For each authenticated support-33 outer word, choose one exact local turnoff.
The selected packets form an isolated episode before a terminal suffix. Every
other active packet enters that suffix. Events for different starts and
labeled packet pairs are disjoint within one selected transition.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
EXPLORATIONS = REPOSITORY / "explorations"
OUTER_RECEIPTS = EXPLORATIONS / "riffle_dp_g4_g1_support33_refutation.json"
ORBIT_RECEIPT = EXPLORATIONS / "riffle_dp_g4_g2_single_packet_orbits.json"
OUTPUT = EXPLORATIONS / "riffle_dp_g4_g2_turnoff_lower_family.json"
PACKET_COUNT = 524_352
INNER_NODES = 32_772
DISTANCE_THRESHOLD = 188_766


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def falling(value: int, count: int) -> int:
    return math.prod(range(value - count + 1, value + 1))


def packet_values(codeword_hex: str) -> list[int]:
    codeword = int(codeword_hex, 16)
    return [
        (codeword >> (4 * packet)) & 0xF
        for packet in range(32)
        if ((codeword >> (4 * packet)) & 0xF) != 0
    ]


def transition_probability(
    *, support: int, pair_count: int, emitted_weight: int
) -> tuple[Fraction, int, int]:
    terminal_nodes = (DISTANCE_THRESHOLD - emitted_weight) // 64
    terminal_packets = 16 * terminal_nodes
    start_count = INNER_NODES - terminal_nodes - 2
    probability = Fraction(
        start_count * pair_count * falling(terminal_packets, support - 2),
        falling(PACKET_COUNT, support),
    )
    return probability, terminal_nodes, start_count


def main() -> None:
    outer = json.loads(OUTER_RECEIPTS.read_text())
    orbit = json.loads(ORBIT_RECEIPT.read_text())
    transitions = orbit["exact_one_packet_turnoffs_below_distance"]
    if len(outer["receipts"]) != 26 or len(transitions) != 6:
        raise RuntimeError("turnoff family inputs have unexpected cardinality")

    members = []
    total = Fraction()
    for receipt in outer["receipts"]:
        values = []
        for local_word in receipt["local_words"]:
            values.extend(packet_values(local_word["codeword_hex"]))
        if len(values) != 33:
            raise RuntimeError("support-33 receipt does not contain 33 active packets")
        counts = Counter(values)
        best = None
        for transition in transitions:
            pair_count = (
                counts[transition["first_value"]]
                * counts[transition["second_value"]]
            )
            if not pair_count:
                continue
            probability, terminal_nodes, start_count = transition_probability(
                support=33,
                pair_count=pair_count,
                emitted_weight=transition["emitted_weight_before_turnoff"],
            )
            candidate = (
                probability,
                transition,
                pair_count,
                terminal_nodes,
                start_count,
            )
            if best is None or candidate[0] > best[0]:
                best = candidate
        if best is None:
            continue
        probability, transition, pair_count, terminal_nodes, start_count = best
        if 64 * terminal_nodes + transition["emitted_weight_before_turnoff"] > DISTANCE_THRESHOLD:
            raise RuntimeError("turnoff family exceeds the distance threshold")
        total += probability
        members.append(
            {
                "outer_family_id": receipt["family_id"],
                "active_packet_value_counts": {
                    str(value): count for value, count in sorted(counts.items())
                },
                "selected_transition": transition,
                "ordered_labeled_packet_pairs": pair_count,
                "episode_start_nodes": start_count,
                "terminal_nodes": terminal_nodes,
                "probability_numerator": str(probability.numerator),
                "probability_denominator": str(probability.denominator),
            }
        )

    payload = {
        "schema": "riffle-dp-g4-g2-turnoff-lower-family-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "EXACT",
        "outer_receipt_sha256": digest(OUTER_RECEIPTS),
        "orbit_receipt_sha256": digest(ORBIT_RECEIPT),
        "member_count": len(members),
        "aggregate_probability": {
            "numerator": str(total.numerator),
            "denominator": str(total.denominator),
            "log2": math.log2(total.numerator) - math.log2(total.denominator),
            "exceeds_2^-40": total > Fraction(1, 1 << 40),
        },
        "members": members,
        "scope": (
            "One isolated exact two-packet turnoff per included support-33 word; "
            "all remaining active packets occupy a terminal suffix."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"outer_words_with_turnoff={len(members)}")
    print(f"aggregate_log2={payload['aggregate_probability']['log2']:.12f}")
    print(f"exceeds_2^-40={payload['aggregate_probability']['exceeds_2^-40']}")
    print(f"output={OUTPUT}")
    print("status=EXACT_G2_TURNOFF_LOWER_FAMILY")


if __name__ == "__main__":
    main()
