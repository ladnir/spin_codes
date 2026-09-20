#!/usr/bin/env python3
"""Replay same-node collision turnoffs and their exact lower family."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
SEARCH = EXPLORATIONS / "riffle_dp_g4_g2_multi_packet_node_turnoffs.json"
LOWER = EXPLORATIONS / "riffle_dp_g4_g2_collision_turnoff_lower_family.json"
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


def drive(nibbles: list[dict]) -> int:
    result = 0
    for nibble in nibbles:
        result |= nibble["value"] << (4 * nibble["slot"])
    return result


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
    search = json.loads(SEARCH.read_text())
    witnesses = search["below_threshold_witnesses"]
    if len(witnesses) != 539 or search["maximum_reset_node"] != INNER_NODES - 1:
        raise RuntimeError("collision verifier: search coverage mismatch")
    columns = systematic_state_columns()
    counts = Counter()
    identities = set()
    for witness in witnesses:
        initial = drive(witness["initial_nibbles"])
        target = drive(witness["reset_nibbles"])
        state = initial
        weight = 0
        for _ in range(witness["reset_node"]):
            weight += accumulate(state).bit_count()
            state = apply(columns, accumulate(state))
        if state != target or weight != witness["emitted_weight_before_turnoff"]:
            raise RuntimeError("collision verifier: recurrence replay mismatch")
        if weight > DISTANCE_THRESHOLD:
            raise RuntimeError("collision verifier: witness exceeds threshold")
        key = (witness["initial_packets"], witness["reset_packets"])
        counts[key] += 1
        identity = (initial, target, witness["reset_node"])
        if identity in identities:
            raise RuntimeError("collision verifier: duplicate witness")
        identities.add(identity)
    expected = {(1, 1): 6, (1, 2): 24, (2, 1): 64, (2, 2): 445}
    if dict(counts) != expected or max(row["reset_node"] for row in witnesses) != 8:
        raise RuntimeError("collision verifier: category counts changed")

    lower = json.loads(LOWER.read_text())
    outer = json.loads(OUTER.read_text())
    if lower["turnoff_receipt_sha256"] != digest(SEARCH):
        raise RuntimeError("collision verifier: lower input hash mismatch")
    outer_by_id = {row["family_id"]: row for row in outer["receipts"]}
    aggregate = Fraction()
    for member in lower["members"]:
        receipt = outer_by_id[member["outer_family_id"]]
        values = []
        for local_word in receipt["local_words"]:
            values.extend(packet_values(local_word["codeword_hex"]))
        value_counts = Counter(values)
        category_totals: dict[int, Fraction] = defaultdict(Fraction)
        for witness in witnesses:
            transition_values = tuple(
                nibble["value"]
                for group in (witness["initial_nibbles"], witness["reset_nibbles"])
                for nibble in group
            )
            outside = len(transition_values)
            choices = labeled_count(value_counts, transition_values)
            if not choices:
                continue
            terminal_nodes = (
                DISTANCE_THRESHOLD - witness["emitted_weight_before_turnoff"]
            ) // 64
            start_count = INNER_NODES - terminal_nodes - witness["reset_node"]
            category_totals[outside] += Fraction(
                start_count * choices * falling(16 * terminal_nodes, 33 - outside),
                falling(PACKET_COUNT, 33),
            )
        selected_outside, selected = max(category_totals.items(), key=lambda row: row[1])
        if selected_outside != member["selected_outside_packets"]:
            raise RuntimeError("collision verifier: selected category mismatch")
        if selected != Fraction(
            int(member["selected_probability_numerator"]),
            int(member["selected_probability_denominator"]),
        ):
            raise RuntimeError("collision verifier: member probability mismatch")
        aggregate += selected
    committed = lower["aggregate_probability"]
    if aggregate != Fraction(int(committed["numerator"]), int(committed["denominator"])):
        raise RuntimeError("collision verifier: aggregate mismatch")
    print("candidate=Riffle DP g=4@g0-v1")
    print("turnoff_counts=1to1:6,1to2:24,2to1:64,2to2:445")
    print("last_turnoff_node=8")
    print(f"collision_lower_family_log2={math.log2(aggregate.numerator)-math.log2(aggregate.denominator):.12f}")
    print("status=EXACT_G2_COLLISION_TURNOFF_ARTIFACTS_VERIFIED")


if __name__ == "__main__":
    main()
