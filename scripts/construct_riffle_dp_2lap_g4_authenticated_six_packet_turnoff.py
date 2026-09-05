#!/usr/bin/env python3
"""Construct and count an authenticated DP-2Lap double-turnoff chain."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from collections import Counter
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns
from search_riffle_dp_2lap_g4_authenticated_double_turnoffs import (
    INNER_NODES,
    OUTPUT as OLD_SEARCH,
    build_apply,
    find_min_span_partition,
    load_episodes,
    nibble_drive,
    packet_values,
    run_lap,
    serialize_counts,
    value_counts,
)


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
OUTER = EXPLORATIONS / "riffle_dp_g4_g1_support33_refutation.json"
SIX_PACKET = EXPLORATIONS / "riffle_dp_2lap_g4_six_packet_two_node_turnoffs.json"
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_authenticated_six_packet_turnoff.json"
TARGET_FAMILY = "terminal-w3-s33-one-i0-tfbfc4e34c1f4398c"
PACKET_POSITIONS = 524_352
DISTANCE_THRESHOLD = 188_766


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log2_interval(value: Fraction, places: int = 12) -> list[str]:
    with localcontext() as context:
        context.prec = 100
        estimate = (
            Decimal(value.numerator).ln() - Decimal(value.denominator).ln()
        ) / Decimal(2).ln()
        unit = Decimal(1).scaleb(-places)
        return [
            format(estimate.quantize(unit, rounding=ROUND_FLOOR), "f"),
            format(estimate.quantize(unit, rounding=ROUND_CEILING), "f"),
        ]


def outer_target() -> tuple[dict, tuple[int, ...]]:
    payload = json.loads(OUTER.read_text())
    receipt = next(
        row for row in payload["receipts"] if row["family_id"] == TARGET_FAMILY
    )
    values: list[int] = []
    for local_word in receipt["local_words"]:
        values.extend(packet_values(local_word["codeword_hex"]))
    counts = value_counts(values)
    if sum(counts) != 33:
        raise RuntimeError("authenticated six-packet chain: target support changed")
    return receipt, counts


def six_packet_episode(systematic_right) -> tuple[list[int], dict]:
    payload = json.loads(SIX_PACKET.read_text())
    witness = next(
        row
        for row in payload["witnesses"]
        if row["packet_value_counts"] == {"5": 3, "13": 3}
    )
    inputs = [0] * (witness["reset_node"] + 1)
    inputs[0] = nibble_drive(witness["initial_nibbles"])
    inputs[witness["reset_node"]] = nibble_drive(witness["reset_nibbles"])
    first_output, first_terminal = run_lap(inputs, 0, systematic_right)
    second_output, second_terminal = run_lap(
        first_output, first_terminal, systematic_right
    )
    first_weight = sum(word.bit_count() for word in first_output)
    second_weight = sum(word.bit_count() for word in second_output)
    if (
        first_terminal != 0
        or second_terminal != 0
        or first_weight != witness["first_lap_weight"]
        or second_weight != witness["two_lap_weight"]
    ):
        raise RuntimeError("authenticated six-packet chain: witness replay failed")
    return inputs, witness


def main() -> None:
    systematic_right = build_apply(systematic_state_columns())
    old_episodes, _ = load_episodes(systematic_right)
    receipt, target = outer_target()
    new_inputs, new_witness = six_packet_episode(systematic_right)
    residual = tuple(3 if value in (5, 13) else 0 for value in range(1, 16))
    complement = tuple(left - right for left, right in zip(target, residual))
    if min(complement) < 0:
        raise RuntimeError("authenticated six-packet chain: invalid residual")
    partition = find_min_span_partition(complement, old_episodes)
    if partition is None:
        raise RuntimeError("authenticated six-packet chain: complement not partitioned")

    episode_inputs: list[tuple[list[int], dict]] = [
        (
            new_inputs,
            {
                "family": "six_packet_two_node",
                "source_index": 0,
                "packet_value_counts": new_witness["packet_value_counts"],
                "span_nodes": len(new_inputs),
                "two_lap_weight": new_witness["two_lap_weight"],
            },
        )
    ]
    for episode_index in partition:
        episode = old_episodes[episode_index]
        episode_inputs.append(
            (
                list(episode.inputs),
                {
                    "family": episode.family,
                    "source_index": episode.source_index,
                    "catalog_index": episode_index,
                    "packet_value_counts": serialize_counts(episode.counts),
                    "span_nodes": episode.span,
                    "two_lap_weight": episode.second_weight,
                },
            )
        )

    inputs: list[int] = []
    members = []
    for local_inputs, member in episode_inputs:
        member = dict(member)
        member["start_node"] = len(inputs)
        members.append(member)
        inputs.extend(local_inputs)
    occupied_span = len(inputs)
    padded = inputs + [0] * (INNER_NODES - occupied_span)
    first_output, first_terminal = run_lap(padded, 0, systematic_right)
    final_output, final_terminal = run_lap(
        first_output, first_terminal, systematic_right
    )
    final_weight = sum(word.bit_count() for word in final_output)
    if first_terminal != 0 or final_terminal != 0:
        raise RuntimeError("authenticated six-packet chain: full terminal mismatch")
    if final_weight != sum(member["two_lap_weight"] for member in members):
        raise RuntimeError("authenticated six-packet chain: additive weight mismatch")
    if final_weight >= DISTANCE_THRESHOLD:
        raise RuntimeError("authenticated six-packet chain: not a counterexample")

    observed_values: list[int] = []
    for word in inputs:
        observed_values.extend(
            (word >> (4 * slot)) & 15
            for slot in range(16)
            if (word >> (4 * slot)) & 15
        )
    if value_counts(observed_values) != target:
        raise RuntimeError("authenticated six-packet chain: packet multiset mismatch")

    episode_count = len(members)
    slack_nodes = INNER_NODES - occupied_span
    gap_placements = math.comb(slack_nodes + episode_count, episode_count)
    labeled_assignments_per_placement = math.prod(
        math.factorial(count) for count in target
    )
    favorable_assignments = gap_placements * labeled_assignments_per_placement
    probability = Fraction(
        favorable_assignments,
        math.prod(range(PACKET_POSITIONS - 33 + 1, PACKET_POSITIONS + 1)),
    )

    payload = {
        "schema": "riffle-dp-2lap-g4-authenticated-six-packet-turnoff-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_LOWER_FAMILY",
        "outer_receipt_sha256": digest(OUTER),
        "old_episode_search_sha256": digest(OLD_SEARCH),
        "six_packet_search_sha256": digest(SIX_PACKET),
        "outer_family_id": receipt["family_id"],
        "outer_parameterization": receipt["outer_family"]["parameterization"],
        "packet_value_counts": serialize_counts(target),
        "episode_count": episode_count,
        "occupied_span_nodes": occupied_span,
        "first_lap_terminal_hex": f"0x{first_terminal:016x}",
        "second_lap_terminal_hex": f"0x{final_terminal:016x}",
        "two_lap_weight": final_weight,
        "distance_threshold": DISTANCE_THRESHOLD,
        "below_distance": True,
        "episodes": members,
        "canonical_nonzero_input_nodes": [
            {"node": node, "drive_hex": f"0x{word:016x}"}
            for node, word in enumerate(inputs)
            if word
        ],
        "lower_family": {
            "description": (
                "Keep the declared episode order and internal packet cells. Insert "
                "arbitrary zero-node gaps before, between, and after the episodes."
            ),
            "slack_nodes": slack_nodes,
            "gap_placements": str(gap_placements),
            "labeled_assignments_per_placement": str(
                labeled_assignments_per_placement
            ),
            "favorable_labeled_assignments": str(favorable_assignments),
            "probability_numerator": str(probability.numerator),
            "probability_denominator": str(probability.denominator),
            "log2_interval": log2_interval(probability),
            "exceeds_2^-40": probability > Fraction(1, 1 << 40),
        },
        "scope_limitation": (
            "The lower family uses one authenticated outer word, one fixed episode "
            "order, and one representative for each episode multiset. It proves "
            "deterministic bad packet orders but need not be the strongest orbit."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"outer_family_id={receipt['family_id']}")
    print(f"episode_count={episode_count}")
    print(f"occupied_span_nodes={occupied_span}")
    print(f"two_lap_weight={final_weight}")
    print(f"lower_family_log2_interval={payload['lower_family']['log2_interval']}")
    print(f"exceeds_2^-40={payload['lower_family']['exceeds_2^-40']}")
    print(f"output={OUTPUT}")
    print("status=EXACT_AUTHENTICATED_SIX_PACKET_TURNOFF_FAMILY")


if __name__ == "__main__":
    main()
