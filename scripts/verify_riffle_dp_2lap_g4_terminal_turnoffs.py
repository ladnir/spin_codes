#!/usr/bin/env python3
"""Replay known first-lap turnoffs under the DP-2Lap recurrence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
THREE_PACKET = EXPLORATIONS / "riffle_dp_g4_g2_three_packet_turnoffs_full.json"
COLLISION = EXPLORATIONS / "riffle_dp_g4_g2_multi_packet_node_turnoffs.json"
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_terminal_turnoff_replay.json"
DISTANCE_THRESHOLD = 188_766
MASK64 = (1 << 64) - 1


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
            raise RuntimeError("DP-2Lap turnoff replay: invalid nibble")
        if (result >> (4 * slot)) & 0xF:
            raise RuntimeError("DP-2Lap turnoff replay: duplicate slot")
        result |= value << (4 * slot)
    return result


def run_lap(
    inputs: list[int], initial_state: int, systematic_right
) -> tuple[list[int], int]:
    state = initial_state
    outputs: list[int] = []
    for value in inputs:
        emitted = accumulate(state ^ value)
        outputs.append(emitted)
        state = systematic_right(emitted)
    return outputs, state


def replay_family(
    name: str,
    episodes: list[tuple[list[int], int]],
    systematic_right,
) -> dict:
    first_weights: list[int] = []
    second_weights: list[int] = []
    second_turnoffs = 0
    for inputs, committed_first_weight in episodes:
        first_output, first_terminal = run_lap(inputs, 0, systematic_right)
        first_weight = sum(value.bit_count() for value in first_output)
        if first_terminal != 0 or first_weight != committed_first_weight:
            raise RuntimeError(f"DP-2Lap turnoff replay: {name} first-lap mismatch")
        second_output, second_terminal = run_lap(first_output, first_terminal, systematic_right)
        second_weight = sum(value.bit_count() for value in second_output)
        if second_weight >= DISTANCE_THRESHOLD:
            raise RuntimeError(f"DP-2Lap turnoff replay: {name} local weight overflow")
        first_weights.append(first_weight)
        second_weights.append(second_weight)
        second_turnoffs += second_terminal == 0

    return {
        "witness_count": len(episodes),
        "terminally_shifted_below_distance_count": sum(
            weight < DISTANCE_THRESHOLD for weight in second_weights
        ),
        "first_lap_weight_minimum": min(first_weights),
        "first_lap_weight_maximum": max(first_weights),
        "two_lap_weight_minimum": min(second_weights),
        "two_lap_weight_maximum": max(second_weights),
        "two_lap_weight_mean": sum(second_weights) / len(second_weights),
        "second_lap_terminal_zero_count": second_turnoffs,
    }


def main() -> None:
    three_payload = json.loads(THREE_PACKET.read_text())
    collision_payload = json.loads(COLLISION.read_text())
    systematic_right = build_apply(systematic_state_columns())

    three_episodes: list[tuple[list[int], int]] = []
    for witness in three_payload["below_threshold_witnesses"]:
        inputs = [0] * (witness["third_node"] + 1)
        inputs[0] = nibble_drive([witness["first"]])
        inputs[witness["second_node"]] = nibble_drive([witness["second"]])
        inputs[witness["third_node"]] = nibble_drive([witness["third"]])
        three_episodes.append((inputs, witness["emitted_weight_before_turnoff"]))

    collision_episodes: list[tuple[list[int], int]] = []
    for witness in collision_payload["below_threshold_witnesses"]:
        inputs = [0] * (witness["reset_node"] + 1)
        inputs[0] = nibble_drive(witness["initial_nibbles"])
        inputs[witness["reset_node"]] = nibble_drive(witness["reset_nibbles"])
        collision_episodes.append((inputs, witness["emitted_weight_before_turnoff"]))

    result = {
        "schema": "riffle-dp-2lap-g4-terminal-turnoff-replay-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT",
        "three_packet_receipt_sha256": digest(THREE_PACKET),
        "collision_receipt_sha256": digest(COLLISION),
        "distance_threshold": DISTANCE_THRESHOLD,
        "placement": (
            "Each isolated episode is shifted so its reset node is the final node. "
            "The preceding nodes and the retained first-lap terminal state are zero."
        ),
        "three_packet": replay_family(
            "three-packet", three_episodes, systematic_right
        ),
        "one_or_two_packet_collision": replay_family(
            "collision", collision_episodes, systematic_right
        ),
        "scope_limitation": (
            "The replay proves that isolated terminal turnoffs remain low-weight after "
            "two laps. It does not show that an authenticated support-33 outer word can "
            "consist of one isolated episode, and it does not reuse the old lower-family "
            "probabilities."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print("candidate=Riffle DP-2Lap g=4")
    print(f"three_packet_count={result['three_packet']['witness_count']}")
    print(
        "three_packet_two_lap_weight_range="
        f"{result['three_packet']['two_lap_weight_minimum']}.."
        f"{result['three_packet']['two_lap_weight_maximum']}"
    )
    print(
        "collision_count="
        f"{result['one_or_two_packet_collision']['witness_count']}"
    )
    print(
        "collision_two_lap_weight_range="
        f"{result['one_or_two_packet_collision']['two_lap_weight_minimum']}.."
        f"{result['one_or_two_packet_collision']['two_lap_weight_maximum']}"
    )
    print(f"output={OUTPUT}")
    print("status=EXACT_DP_2LAP_TERMINAL_TURNOFF_REPLAY")


if __name__ == "__main__":
    main()
