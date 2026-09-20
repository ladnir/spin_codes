#!/usr/bin/env python3
"""Deterministic local search for lifted-window low-weight witnesses."""

from __future__ import annotations

import json
import random
from pathlib import Path

from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    observation_rows,
)
from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "constructions"
    / "riffle_packetmul_wrapmul_2lap_g4"
    / "receipts"
    / "goal04_lifted_window_search.json"
)
WINDOWS = (8, 9, 10, 11, 12, 14, 16, 18, 24, 32, 48, 64)
RANDOM_RESTARTS = 512
SEED = 0x4C49465445444F52
SAT_RECEIPTS = tuple(
    ROOT
    / "constructions"
    / "riffle_packetmul_wrapmul_2lap_g4"
    / "receipts"
    / f"goal04_lifted_window_{nodes:02d}_sat.json"
    for nodes in (6, 7, 9)
)


def descend(word: int, state: int, rows: list[int]) -> tuple[int, int]:
    weight = word.bit_count()
    while True:
        best_weight = weight
        best_coordinate = -1
        for coordinate, row in enumerate(rows):
            candidate = (word ^ row).bit_count()
            if candidate < best_weight:
                best_weight = candidate
                best_coordinate = coordinate
        if best_coordinate < 0:
            return word, state
        word ^= rows[best_coordinate]
        state ^= 1 << best_coordinate
        weight = best_weight


def word_from_state(state: int, rows: list[int]) -> int:
    word = 0
    while state:
        bit = (state & -state).bit_length() - 1
        word ^= rows[bit]
        state &= state - 1
    return word


def main() -> None:
    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(max(WINDOWS), parity)
    rng = random.Random(SEED)
    known_states = []
    for path in SAT_RECEIPTS:
        if not path.exists():
            continue
        receipt = json.loads(path.read_text())
        if receipt.get("initial_a_hex") and receipt.get("initial_b_hex"):
            known_states.append(
                int(receipt["initial_a_hex"], 16)
                | (int(receipt["initial_b_hex"], 16) << 64)
            )
    results = []
    for nodes in WINDOWS:
        mask = (1 << (64 * nodes)) - 1
        rows = [row & mask for row in full_rows]
        starts = [1 << bit for bit in range(128)] + known_states
        starts.extend(rng.getrandbits(128) or 1 for _ in range(RANDOM_RESTARTS))
        best_word = mask
        best_state = 0
        for state in starts:
            word, final_state = descend(word_from_state(state, rows), state, rows)
            if 0 < word.bit_count() < best_word.bit_count():
                best_word = word
                best_state = final_state
        node_weights = [
            ((best_word >> (64 * node)) & ((1 << 64) - 1)).bit_count()
            for node in range(nodes)
        ]
        result = {
            "nodes": nodes,
            "best_weight": best_word.bit_count(),
            "weight_per_node": best_word.bit_count() / nodes,
            "initial_a_hex": f"0x{best_state & ((1 << 64) - 1):016x}",
            "initial_b_hex": f"0x{best_state >> 64:016x}",
            "node_weights": node_weights,
        }
        results.append(result)
        print(
            f"nodes={nodes} best_weight={result['best_weight']} "
            f"weight_per_node={result['weight_per_node']:.6f}"
        )

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-lifted-window-search-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "DIAGNOSTIC_LOCAL_SEARCH_WITH_EXACT_WITNESSES",
        "seed_hex": f"0x{SEED:016x}",
        "random_restarts_per_window": RANDOM_RESTARTS,
        "target_weight_per_node": 188_766 / 32_772,
        "windows": results,
        "scope_limitation": (
            "Every witness replays exactly. Coordinate descent does not prove "
            "minimum distance or exclude lower-weight words."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print("status=DIAGNOSTIC_LIFTED_WINDOW_SEARCH")


if __name__ == "__main__":
    main()
