#!/usr/bin/env python3
"""Independently verify the bounded authenticated double-turnoff search."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
OUTER = EXPLORATIONS / "riffle_dp_g4_g1_support33_refutation.json"
THREE = EXPLORATIONS / "riffle_dp_g4_g2_three_packet_turnoffs_full.json"
COLLISION = EXPLORATIONS / "riffle_dp_g4_g2_multi_packet_node_turnoffs.json"
SEARCH = EXPLORATIONS / "riffle_dp_2lap_g4_authenticated_double_turnoff_search.json"
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_authenticated_double_turnoff_verification.json"
MASK64 = (1 << 64) - 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accumulate(value: int) -> int:
    value ^= value << 1
    value ^= value << 2
    value ^= value << 4
    value ^= value << 8
    value ^= value << 16
    value ^= value << 32
    return value & MASK64


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        bit = (value & -value).bit_length() - 1
        result ^= columns[bit]
        value &= value - 1
    return result


def drive(nibbles: list[dict[str, int]]) -> int:
    result = 0
    for nibble in nibbles:
        if (result >> (4 * nibble["slot"])) & 15:
            raise RuntimeError("double-turnoff verification: duplicate slot")
        result |= nibble["value"] << (4 * nibble["slot"])
    return result


def terminal_after_two_laps(inputs: list[int], columns: tuple[int, ...]) -> tuple[int, int]:
    state = 0
    first_output = []
    for word in inputs:
        emitted = accumulate(state ^ word)
        first_output.append(emitted)
        state = apply_columns(columns, emitted)
    first_terminal = state
    for word in first_output:
        emitted = accumulate(state ^ word)
        state = apply_columns(columns, emitted)
    return first_terminal, state


def count_vector(nibbles: list[dict[str, int]]) -> tuple[int, ...]:
    counts = Counter(nibble["value"] for nibble in nibbles)
    return tuple(counts[value] for value in range(1, 16))


def episode_vectors(columns: tuple[int, ...]) -> tuple[set[tuple[int, ...]], dict[str, int]]:
    vectors: set[tuple[int, ...]] = set()
    totals = {"three_packet": 0, "collision": 0}
    three = json.loads(THREE.read_text())
    for witness in three["below_threshold_witnesses"]:
        inputs = [0] * (witness["third_node"] + 1)
        nibbles = [witness["first"], witness["second"], witness["third"]]
        inputs[0] = drive([nibbles[0]])
        inputs[witness["second_node"]] = drive([nibbles[1]])
        inputs[witness["third_node"]] = drive([nibbles[2]])
        if terminal_after_two_laps(inputs, columns) == (0, 0):
            totals["three_packet"] += 1
            vectors.add(count_vector(nibbles))

    collision = json.loads(COLLISION.read_text())
    for witness in collision["below_threshold_witnesses"]:
        inputs = [0] * (witness["reset_node"] + 1)
        inputs[0] = drive(witness["initial_nibbles"])
        inputs[witness["reset_node"]] = drive(witness["reset_nibbles"])
        nibbles = witness["initial_nibbles"] + witness["reset_nibbles"]
        if terminal_after_two_laps(inputs, columns) == (0, 0):
            totals["collision"] += 1
            vectors.add(count_vector(nibbles))
    return vectors, totals


def outer_vectors() -> dict[tuple[int, ...], list[str]]:
    result: dict[tuple[int, ...], list[str]] = {}
    payload = json.loads(OUTER.read_text())
    for receipt in payload["receipts"]:
        values = []
        for local_word in receipt["local_words"]:
            word = int(local_word["codeword_hex"], 16)
            values.extend((word >> shift) & 15 for shift in range(0, 128, 4) if (word >> shift) & 15)
        counts = Counter(values)
        vector = tuple(counts[value] for value in range(1, 16))
        result.setdefault(vector, []).append(receipt["family_id"])
    return result


def closure_data(
    target: tuple[int, ...], episodes: set[tuple[int, ...]]
) -> tuple[int, int, bool, int, list[dict[str, int]]]:
    allowed = [episode for episode in episodes if all(x <= y for x, y in zip(episode, target))]
    zero = (0,) * 15
    reached = {zero}
    frontier = {zero}
    while frontier:
        following = set()
        for state in frontier:
            for episode in allowed:
                candidate = tuple(x + y for x, y in zip(state, episode))
                if all(x <= y for x, y in zip(candidate, target)) and candidate not in reached:
                    following.add(candidate)
        reached.update(following)
        frontier = following
    residuals = [
        tuple(left - right for left, right in zip(target, state))
        for state in reached
        if state != target
    ]
    minimum = min(sum(residual) for residual in residuals)
    minimum_profiles = sorted(
        {
            residual
            for residual in residuals
            if sum(residual) == minimum
        }
    )
    serialized = [
        {
            str(value): count
            for value, count in enumerate(profile, 1)
            if count
        }
        for profile in minimum_profiles
    ]
    return len(allowed), len(reached), target in reached, minimum, serialized


def main() -> None:
    columns = systematic_state_columns()
    episodes, raw_counts = episode_vectors(columns)
    targets = outer_vectors()
    rows = []
    for target_index, (target, family_ids) in enumerate(sorted(targets.items())):
        allowed, reachable, hit, minimum, residual_profiles = closure_data(
            target, episodes
        )
        rows.append(
            {
                "target_index": target_index,
                "outer_family_ids": family_ids,
                "allowed_episode_multisets": allowed,
                "reachable_submultisets": reachable,
                "target_reachable": hit,
                "minimum_uncovered_packets": minimum,
                "minimum_residual_profiles": residual_profiles,
            }
        )

    search = json.loads(SEARCH.read_text())
    if raw_counts != search["raw_double_turnoff_counts"]:
        raise RuntimeError("double-turnoff verification: raw catalog mismatch")
    if len(episodes) != search["distinct_packet_multiset_episode_count"]:
        raise RuntimeError("double-turnoff verification: distinct catalog mismatch")
    if len(targets) != search["distinct_outer_packet_multiset_count"]:
        raise RuntimeError("double-turnoff verification: outer target mismatch")
    if any(row["target_reachable"] for row in rows):
        raise RuntimeError("double-turnoff verification: unexpected reachable target")
    if search["partitioned_outer_multiset_count"] != 0:
        raise RuntimeError("double-turnoff verification: search result mismatch")

    payload = {
        "schema": "riffle-dp-2lap-g4-authenticated-double-turnoff-verification-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_INDEPENDENT_VERIFICATION",
        "search_receipt_sha256": digest(SEARCH),
        "raw_double_turnoff_counts": raw_counts,
        "distinct_packet_multiset_episode_count": len(episodes),
        "distinct_outer_packet_multiset_count": len(targets),
        "reachable_outer_multiset_count": 0,
        "global_minimum_uncovered_packets": min(
            row["minimum_uncovered_packets"] for row in rows
        ),
        "method": (
            "Recompute both terminal states by direct state-column evaluation. Then "
            "enumerate the additive closure of allowed episode multisets below each "
            "authenticated target. This differs from the recursive subtraction search."
        ),
        "targets": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"raw_double_turnoffs={sum(raw_counts.values())}")
    print(f"distinct_episode_multisets={len(episodes)}")
    print(f"distinct_outer_multisets={len(targets)}")
    print("reachable_outer_multisets=0")
    print(f"output={OUTPUT}")
    print("status=EXACT_INDEPENDENT_VERIFICATION")


if __name__ == "__main__":
    main()
