#!/usr/bin/env python3
"""Search authenticated support-33 words for concatenated double turnoffs."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
OUTER_RECEIPTS = EXPLORATIONS / "riffle_dp_g4_g1_support33_refutation.json"
THREE_PACKET = EXPLORATIONS / "riffle_dp_g4_g2_three_packet_turnoffs_full.json"
COLLISION = EXPLORATIONS / "riffle_dp_g4_g2_multi_packet_node_turnoffs.json"
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_authenticated_double_turnoff_search.json"
INNER_NODES = 32_772
DISTANCE_THRESHOLD = 188_766
MASK64 = (1 << 64) - 1


@dataclass(frozen=True)
class Episode:
    family: str
    source_index: int
    counts: tuple[int, ...]
    inputs: tuple[int, ...]
    first_weight: int
    second_weight: int

    @property
    def packet_count(self) -> int:
        return sum(self.counts)

    @property
    def span(self) -> int:
        return len(self.inputs)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accumulate(value: int) -> int:
    for shift in (1, 2, 4, 8, 16, 32):
        value ^= value << shift
    return value & MASK64


def build_apply(columns: tuple[int, ...]):
    tables: list[list[int]] = []
    for byte_index in range(8):
        table: list[int] = []
        for byte in range(256):
            value = 0
            for bit in range(8):
                if (byte >> bit) & 1:
                    value ^= columns[8 * byte_index + bit]
            table.append(value)
        tables.append(table)

    def apply(value: int) -> int:
        return (
            tables[0][value & 0xFF]
            ^ tables[1][(value >> 8) & 0xFF]
            ^ tables[2][(value >> 16) & 0xFF]
            ^ tables[3][(value >> 24) & 0xFF]
            ^ tables[4][(value >> 32) & 0xFF]
            ^ tables[5][(value >> 40) & 0xFF]
            ^ tables[6][(value >> 48) & 0xFF]
            ^ tables[7][value >> 56]
        )

    return apply


def nibble_drive(nibbles: list[dict[str, int]]) -> int:
    result = 0
    for nibble in nibbles:
        value = nibble["value"]
        slot = nibble["slot"]
        if not 0 < value < 16 or not 0 <= slot < 16:
            raise RuntimeError("double-turnoff search: invalid nibble")
        if (result >> (4 * slot)) & 0xF:
            raise RuntimeError("double-turnoff search: duplicate slot")
        result |= value << (4 * slot)
    return result


def packet_values(codeword_hex: str) -> list[int]:
    word = int(codeword_hex, 16)
    return [
        (word >> (4 * packet)) & 0xF
        for packet in range(32)
        if ((word >> (4 * packet)) & 0xF) != 0
    ]


def value_counts(values: list[int]) -> tuple[int, ...]:
    counts = Counter(values)
    return tuple(counts[value] for value in range(1, 16))


def run_lap(
    inputs: tuple[int, ...] | list[int], initial_state: int, systematic_right
) -> tuple[list[int], int]:
    state = initial_state
    outputs: list[int] = []
    for value in inputs:
        emitted = accumulate(state ^ value)
        outputs.append(emitted)
        state = systematic_right(emitted)
    return outputs, state


def make_episode(
    family: str,
    source_index: int,
    inputs: list[int],
    values: list[int],
    committed_first_weight: int,
    systematic_right,
) -> Episode | None:
    first_output, first_terminal = run_lap(inputs, 0, systematic_right)
    first_weight = sum(word.bit_count() for word in first_output)
    if first_terminal != 0 or first_weight != committed_first_weight:
        raise RuntimeError("double-turnoff search: invalid first-lap receipt")
    second_output, second_terminal = run_lap(first_output, first_terminal, systematic_right)
    if second_terminal != 0:
        return None
    return Episode(
        family=family,
        source_index=source_index,
        counts=value_counts(values),
        inputs=tuple(inputs),
        first_weight=first_weight,
        second_weight=sum(word.bit_count() for word in second_output),
    )


def load_episodes(systematic_right) -> tuple[list[Episode], dict[str, int]]:
    three_payload = json.loads(THREE_PACKET.read_text())
    collision_payload = json.loads(COLLISION.read_text())
    episodes: list[Episode] = []
    raw_double_counts = {"three_packet": 0, "collision": 0}

    for source_index, witness in enumerate(three_payload["below_threshold_witnesses"]):
        inputs = [0] * (witness["third_node"] + 1)
        inputs[0] = nibble_drive([witness["first"]])
        inputs[witness["second_node"]] = nibble_drive([witness["second"]])
        inputs[witness["third_node"]] = nibble_drive([witness["third"]])
        values = [
            witness["first"]["value"],
            witness["second"]["value"],
            witness["third"]["value"],
        ]
        episode = make_episode(
            "three_packet",
            source_index,
            inputs,
            values,
            witness["emitted_weight_before_turnoff"],
            systematic_right,
        )
        if episode is not None:
            raw_double_counts["three_packet"] += 1
            episodes.append(episode)

    for source_index, witness in enumerate(collision_payload["below_threshold_witnesses"]):
        inputs = [0] * (witness["reset_node"] + 1)
        inputs[0] = nibble_drive(witness["initial_nibbles"])
        inputs[witness["reset_node"]] = nibble_drive(witness["reset_nibbles"])
        nibbles = witness["initial_nibbles"] + witness["reset_nibbles"]
        episode = make_episode(
            "collision",
            source_index,
            inputs,
            [nibble["value"] for nibble in nibbles],
            witness["emitted_weight_before_turnoff"],
            systematic_right,
        )
        if episode is not None:
            raw_double_counts["collision"] += 1
            episodes.append(episode)

    # A count vector can be reused at any translated node. Keep the representative
    # with the shortest span, then the smallest final weight.
    representatives: dict[tuple[int, ...], Episode] = {}
    for episode in episodes:
        old = representatives.get(episode.counts)
        if old is None or (episode.span, episode.second_weight) < (
            old.span,
            old.second_weight,
        ):
            representatives[episode.counts] = episode
    result = sorted(
        representatives.values(),
        key=lambda episode: (episode.packet_count, episode.counts, episode.span),
    )
    return result, raw_double_counts


def subtract_counts(
    remaining: tuple[int, ...], used: tuple[int, ...]
) -> tuple[int, ...] | None:
    difference = tuple(left - right for left, right in zip(remaining, used))
    return difference if min(difference) >= 0 else None


def find_min_span_partition(
    target: tuple[int, ...], episodes: list[Episode]
) -> tuple[int, ...] | None:
    by_value = [
        tuple(index for index, episode in enumerate(episodes) if episode.counts[value])
        for value in range(15)
    ]

    @lru_cache(maxsize=None)
    def solve(remaining: tuple[int, ...]) -> tuple[int, int, tuple[int, ...]] | None:
        if not any(remaining):
            return (0, 0, ())
        active_values = [value for value, count in enumerate(remaining) if count]
        pivot = min(active_values, key=lambda value: len(by_value[value]))
        best: tuple[int, int, tuple[int, ...]] | None = None
        for episode_index in by_value[pivot]:
            episode = episodes[episode_index]
            difference = subtract_counts(remaining, episode.counts)
            if difference is None:
                continue
            suffix = solve(difference)
            if suffix is None:
                continue
            candidate = (
                episode.span + suffix[0],
                episode.second_weight + suffix[1],
                (episode_index,) + suffix[2],
            )
            if candidate[:2] < (best[:2] if best is not None else (1 << 60, 1 << 60)):
                best = candidate
        return best

    answer = solve(target)
    if answer is None or answer[0] > INNER_NODES:
        return None
    return answer[2]


def serialize_counts(counts: tuple[int, ...]) -> dict[str, int]:
    return {str(value): count for value, count in enumerate(counts, 1) if count}


def replay_partition(indices: tuple[int, ...], episodes: list[Episode], systematic_right) -> dict:
    inputs: list[int] = []
    members = []
    for episode_index in indices:
        episode = episodes[episode_index]
        start_node = len(inputs)
        inputs.extend(episode.inputs)
        members.append(
            {
                "catalog_index": episode_index,
                "family": episode.family,
                "source_index": episode.source_index,
                "start_node": start_node,
                "span_nodes": episode.span,
                "packet_count": episode.packet_count,
                "packet_value_counts": serialize_counts(episode.counts),
                "first_lap_weight": episode.first_weight,
                "two_lap_weight": episode.second_weight,
            }
        )
    padded_inputs = tuple(inputs + [0] * (INNER_NODES - len(inputs)))
    first_output, first_terminal = run_lap(padded_inputs, 0, systematic_right)
    final_output, final_terminal = run_lap(first_output, first_terminal, systematic_right)
    first_weight = sum(word.bit_count() for word in first_output)
    final_weight = sum(word.bit_count() for word in final_output)
    if first_terminal != 0 or final_terminal != 0:
        raise RuntimeError("double-turnoff search: concatenated terminal mismatch")
    if first_weight != sum(member["first_lap_weight"] for member in members):
        raise RuntimeError("double-turnoff search: concatenated first weight mismatch")
    if final_weight != sum(member["two_lap_weight"] for member in members):
        raise RuntimeError("double-turnoff search: concatenated final weight mismatch")
    return {
        "episode_count": len(members),
        "occupied_span_nodes": len(inputs),
        "first_lap_terminal_hex": f"0x{first_terminal:016x}",
        "second_lap_terminal_hex": f"0x{final_terminal:016x}",
        "first_lap_weight": first_weight,
        "two_lap_weight": final_weight,
        "below_distance": final_weight < DISTANCE_THRESHOLD,
        "episodes": members,
        "nonzero_input_nodes": [
            {"node": node, "drive_hex": f"0x{word:016x}"}
            for node, word in enumerate(inputs)
            if word
        ],
    }


def main() -> None:
    systematic_right = build_apply(systematic_state_columns())
    episodes, raw_double_counts = load_episodes(systematic_right)
    outer = json.loads(OUTER_RECEIPTS.read_text())
    if len(outer["receipts"]) != 26:
        raise RuntimeError("double-turnoff search: outer receipt cardinality changed")

    distinct_targets: dict[tuple[int, ...], list[str]] = {}
    for receipt in outer["receipts"]:
        values: list[int] = []
        for local_word in receipt["local_words"]:
            values.extend(packet_values(local_word["codeword_hex"]))
        counts = value_counts(values)
        if sum(counts) != 33:
            raise RuntimeError("double-turnoff search: outer support changed")
        distinct_targets.setdefault(counts, []).append(receipt["family_id"])

    searches = []
    hit_count = 0
    for target_index, (counts, family_ids) in enumerate(sorted(distinct_targets.items())):
        partition = find_min_span_partition(counts, episodes)
        entry = {
            "target_index": target_index,
            "outer_family_ids": family_ids,
            "packet_value_counts": serialize_counts(counts),
            "partition_found": partition is not None,
        }
        if partition is not None:
            hit_count += 1
            entry["witness"] = replay_partition(partition, episodes, systematic_right)
        searches.append(entry)

    payload = {
        "schema": "riffle-dp-2lap-g4-authenticated-double-turnoff-search-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_BOUNDED_MECHANISM_SEARCH",
        "outer_receipt_sha256": digest(OUTER_RECEIPTS),
        "three_packet_receipt_sha256": digest(THREE_PACKET),
        "collision_receipt_sha256": digest(COLLISION),
        "inner_nodes": INNER_NODES,
        "distance_threshold": DISTANCE_THRESHOLD,
        "raw_double_turnoff_counts": raw_double_counts,
        "distinct_packet_multiset_episode_count": len(episodes),
        "outer_receipt_count": len(outer["receipts"]),
        "distinct_outer_packet_multiset_count": len(distinct_targets),
        "partitioned_outer_multiset_count": hit_count,
        "searches": searches,
        "interpretation": (
            "Each reported partition concatenates known episodes. Every episode starts "
            "and ends both laps in state zero. The exact replay pads the construction "
            "to the full chain and checks both terminal states and the final weight."
        ),
        "scope_limitation": (
            "A miss excludes only concatenations of the known two-to-four-packet "
            "turnoff catalog. It does not exclude larger primitive turnoffs or "
            "terminal-state cancellation between episodes."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"raw_double_turnoffs={sum(raw_double_counts.values())}")
    print(f"distinct_episode_multisets={len(episodes)}")
    print(f"distinct_outer_multisets={len(distinct_targets)}")
    print(f"partitioned_outer_multisets={hit_count}")
    print(f"output={OUTPUT}")
    print("status=EXACT_BOUNDED_MECHANISM_SEARCH")


if __name__ == "__main__":
    main()
