#!/usr/bin/env python3
"""Exact support-33 lower family from three-packet inner turnoffs."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
OUTER_RECEIPTS = EXPLORATIONS / "riffle_dp_g4_g1_support33_refutation.json"
TURNOFF_RECEIPT = EXPLORATIONS / "riffle_dp_g4_g2_three_packet_turnoffs_full.json"
OUTPUT = EXPLORATIONS / "riffle_dp_g4_g2_three_packet_lower_family.json"
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


def labeled_triple_count(counts: Counter[int], values: tuple[int, int, int]) -> int:
    remaining = counts.copy()
    result = 1
    for value in values:
        result *= remaining[value]
        remaining[value] -= 1
    return result


def witness_probability(
    *, support: int, labeled_triples: int, third_node: int, emitted_weight: int
) -> tuple[Fraction, int, int]:
    terminal_nodes = (DISTANCE_THRESHOLD - emitted_weight) // 64
    terminal_packets = 16 * terminal_nodes
    start_count = INNER_NODES - terminal_nodes - third_node
    if start_count <= 0:
        return Fraction(), terminal_nodes, 0
    probability = Fraction(
        start_count
        * labeled_triples
        * falling(terminal_packets, support - 3),
        falling(PACKET_COUNT, support),
    )
    return probability, terminal_nodes, start_count


def main() -> None:
    outer = json.loads(OUTER_RECEIPTS.read_text())
    turnoffs = json.loads(TURNOFF_RECEIPT.read_text())
    witnesses = turnoffs["below_threshold_witnesses"]
    if len(outer["receipts"]) != 26:
        raise RuntimeError("three-packet family: unexpected outer receipt count")
    if turnoffs["maximum_third_node"] != INNER_NODES - 1:
        raise RuntimeError("three-packet family: turnoff search does not cover the chain")
    if len(witnesses) != 27 or turnoffs["algebraic_turnoff_solutions"] != 27:
        raise RuntimeError("three-packet family: unexpected turnoff count")

    total = Fraction()
    members = []
    for receipt in outer["receipts"]:
        values = []
        for local_word in receipt["local_words"]:
            values.extend(packet_values(local_word["codeword_hex"]))
        if len(values) != 33:
            raise RuntimeError("three-packet family: outer word support is not 33")
        counts = Counter(values)
        rows = []
        outer_total = Fraction()
        for witness_index, witness in enumerate(witnesses):
            transition_values = tuple(
                witness[name]["value"] for name in ("first", "second", "third")
            )
            triple_count = labeled_triple_count(counts, transition_values)
            if not triple_count:
                continue
            probability, terminal_nodes, start_count = witness_probability(
                support=33,
                labeled_triples=triple_count,
                third_node=witness["third_node"],
                emitted_weight=witness["emitted_weight_before_turnoff"],
            )
            if not probability:
                continue
            if (
                64 * terminal_nodes
                + witness["emitted_weight_before_turnoff"]
                > DISTANCE_THRESHOLD
            ):
                raise RuntimeError("three-packet family exceeds distance threshold")
            outer_total += probability
            rows.append(
                {
                    "witness_index": witness_index,
                    "transition_values": list(transition_values),
                    "ordered_labeled_packet_triples": triple_count,
                    "episode_start_nodes": start_count,
                    "terminal_nodes": terminal_nodes,
                    "probability_numerator": str(probability.numerator),
                    "probability_denominator": str(probability.denominator),
                }
            )
        if not rows:
            continue
        total += outer_total
        members.append(
            {
                "outer_family_id": receipt["family_id"],
                "active_packet_value_counts": {
                    str(value): count for value, count in sorted(counts.items())
                },
                "matched_witness_count": len(rows),
                "probability_numerator": str(outer_total.numerator),
                "probability_denominator": str(outer_total.denominator),
                "witness_rows": rows,
            }
        )

    payload = {
        "schema": "riffle-dp-g4-g2-three-packet-lower-family-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "EXACT",
        "outer_receipt_sha256": digest(OUTER_RECEIPTS),
        "turnoff_receipt_sha256": digest(TURNOFF_RECEIPT),
        "outer_member_count": len(members),
        "matched_outer_witness_pairs": sum(
            member["matched_witness_count"] for member in members
        ),
        "aggregate_probability": {
            "numerator": str(total.numerator),
            "denominator": str(total.denominator),
            "log2": math.log2(total.numerator) - math.log2(total.denominator),
            "exceeds_2^-40": total > Fraction(1, 1 << 40),
        },
        "members": members,
        "scope": (
            "Exactly three labeled active packets occupy the authenticated "
            "distinct-node turnoff pattern. Every other active packet occupies "
            "a terminal suffix. The enumerated placement events are disjoint."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"outer_words_with_turnoff={len(members)}")
    print(f"matched_outer_witness_pairs={payload['matched_outer_witness_pairs']}")
    print(f"aggregate_log2={payload['aggregate_probability']['log2']:.12f}")
    print(f"exceeds_2^-40={payload['aggregate_probability']['exceeds_2^-40']}")
    print(f"output={OUTPUT}")
    print("status=EXACT_G2_THREE_PACKET_LOWER_FAMILY")


if __name__ == "__main__":
    main()
