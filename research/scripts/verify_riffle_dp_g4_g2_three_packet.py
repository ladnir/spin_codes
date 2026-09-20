#!/usr/bin/env python3
"""Replay three-packet turnoffs and their exact support-33 lower family."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
FULL = EXPLORATIONS / "riffle_dp_g4_g2_three_packet_turnoffs_full.json"
PREFIX = EXPLORATIONS / "riffle_dp_g4_g2_three_packet_turnoffs_span6000.json"
LOWER = EXPLORATIONS / "riffle_dp_g4_g2_three_packet_lower_family.json"
OUTER = EXPLORATIONS / "riffle_dp_g4_g1_support33_refutation.json"
PACKET_COUNT = 524_352
INNER_NODES = 32_772
DISTANCE_THRESHOLD = 188_766
MASK64 = (1 << 64) - 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accumulate(value: int) -> int:
    for shift in (1, 2, 4, 8, 16, 32):
        value ^= value << shift
    return value & MASK64


def apply(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        bit = value & -value
        result ^= columns[bit.bit_length() - 1]
        value ^= bit
    return result


def drive(row: dict) -> int:
    return row["value"] << (4 * row["slot"])


def falling(value: int, count: int) -> int:
    return math.prod(range(value - count + 1, value + 1))


def packet_values(codeword_hex: str) -> list[int]:
    word = int(codeword_hex, 16)
    return [
        (word >> (4 * packet)) & 0xF
        for packet in range(32)
        if ((word >> (4 * packet)) & 0xF) != 0
    ]


def labeled_count(counts: Counter[int], values: tuple[int, int, int]) -> int:
    remaining = counts.copy()
    result = 1
    for value in values:
        result *= remaining[value]
        remaining[value] -= 1
    return result


def main() -> None:
    full = json.loads(FULL.read_text())
    prefix = json.loads(PREFIX.read_text())
    witnesses = full["below_threshold_witnesses"]
    if full["maximum_third_node"] != INNER_NODES - 1:
        raise RuntimeError("three-packet verifier: full-chain coverage mismatch")
    expected_queries = (INNER_NODES - 2) * 240 * 240
    expected_left = (INNER_NODES - 2) * 240
    if full["lookup_queries"] != expected_queries or full["left_entries"] != expected_left:
        raise RuntimeError("three-packet verifier: enumeration cardinality mismatch")
    if witnesses != prefix["below_threshold_witnesses"]:
        raise RuntimeError("three-packet verifier: prefix and full witnesses differ")
    if len(witnesses) != 27 or full["algebraic_turnoff_solutions"] != 27:
        raise RuntimeError("three-packet verifier: turnoff count mismatch")

    columns = systematic_state_columns()
    identities = set()
    for witness in witnesses:
        first_node = 0
        second_node = witness["second_node"]
        third_node = witness["third_node"]
        if not 0 < second_node < third_node < INNER_NODES:
            raise RuntimeError("three-packet verifier: invalid node order")
        inputs = {
            first_node: drive(witness["first"]),
            second_node: drive(witness["second"]),
            third_node: drive(witness["third"]),
        }
        state = 0
        emitted_weight = 0
        for node in range(third_node + 1):
            emitted = accumulate(state ^ inputs.get(node, 0))
            if node < third_node:
                emitted_weight += emitted.bit_count()
            elif emitted != 0:
                raise RuntimeError("three-packet verifier: final input does not turn off")
            state = apply(columns, emitted)
        if state != 0 or emitted_weight != witness["emitted_weight_before_turnoff"]:
            raise RuntimeError("three-packet verifier: exact replay mismatch")
        identity = (
            witness["first"]["slot"],
            witness["first"]["value"],
            witness["second"]["slot"],
            witness["second"]["value"],
            witness["third"]["slot"],
            witness["third"]["value"],
            second_node,
            third_node,
        )
        if identity in identities:
            raise RuntimeError("three-packet verifier: duplicate witness")
        identities.add(identity)
    if max(witness["third_node"] for witness in witnesses) != 6:
        raise RuntimeError("three-packet verifier: last witness node changed")

    lower = json.loads(LOWER.read_text())
    outer = json.loads(OUTER.read_text())
    if lower["turnoff_receipt_sha256"] != digest(FULL):
        raise RuntimeError("three-packet verifier: lower-family input hash mismatch")
    outer_by_id = {row["family_id"]: row for row in outer["receipts"]}
    total = Fraction()
    matched = 0
    for member in lower["members"]:
        receipt = outer_by_id[member["outer_family_id"]]
        values = []
        for local_word in receipt["local_words"]:
            values.extend(packet_values(local_word["codeword_hex"]))
        counts = Counter(values)
        member_total = Fraction()
        for row in member["witness_rows"]:
            witness = witnesses[row["witness_index"]]
            transition_values = tuple(
                witness[name]["value"] for name in ("first", "second", "third")
            )
            triples = labeled_count(counts, transition_values)
            terminal_nodes = (
                DISTANCE_THRESHOLD - witness["emitted_weight_before_turnoff"]
            ) // 64
            start_count = INNER_NODES - terminal_nodes - witness["third_node"]
            probability = Fraction(
                start_count * triples * falling(16 * terminal_nodes, 30),
                falling(PACKET_COUNT, 33),
            )
            if triples != row["ordered_labeled_packet_triples"]:
                raise RuntimeError("three-packet verifier: labeled count mismatch")
            if probability != Fraction(
                int(row["probability_numerator"]),
                int(row["probability_denominator"]),
            ):
                raise RuntimeError("three-packet verifier: row probability mismatch")
            member_total += probability
            matched += 1
        if member_total != Fraction(
            int(member["probability_numerator"]),
            int(member["probability_denominator"]),
        ):
            raise RuntimeError("three-packet verifier: member probability mismatch")
        total += member_total
    committed_total = lower["aggregate_probability"]
    if total != Fraction(
        int(committed_total["numerator"]),
        int(committed_total["denominator"]),
    ):
        raise RuntimeError("three-packet verifier: aggregate mismatch")
    if matched != lower["matched_outer_witness_pairs"]:
        raise RuntimeError("three-packet verifier: matched-pair count mismatch")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"full_chain_algebraic_turnoffs={len(witnesses)}")
    print(f"last_turnoff_node={max(witness['third_node'] for witness in witnesses)}")
    print(f"matched_outer_witness_pairs={matched}")
    print(f"lower_family_log2={math.log2(total.numerator)-math.log2(total.denominator):.12f}")
    print("status=EXACT_G2_THREE_PACKET_ARTIFACTS_VERIFIED")


if __name__ == "__main__":
    main()
