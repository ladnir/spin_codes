#!/usr/bin/env python3
"""Deterministic search for repeated low windows on one lifted orbit."""

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
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
OUTPUT = CANDIDATE / "receipts" / "goal05_repeated_low_window_search.json"
WINDOW = 9
THRESHOLD = 52
OFFSETS = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 24, 32)
RESTARTS = 4096
SEED = 0x52455455524E4C4F
MASK64 = (1 << 64) - 1


def word_from_state(state: int, rows: list[int]) -> int:
    word = 0
    while state:
        low = state & -state
        word ^= rows[low.bit_length() - 1]
        state ^= low
    return word


def score(word: int, first_mask: int, second_mask: int, offset: int) -> tuple[int, int]:
    first = (word & first_mask).bit_count()
    second = ((word & second_mask) >> (64 * offset)).bit_count()
    return max(first, second), first + second


def main() -> None:
    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(max(OFFSETS) + WINDOW, parity)
    known_states = []
    for nodes in (6, 7, 9):
        path = (
            CANDIDATE
            / "receipts"
            / f"goal04_lifted_window_{nodes:02d}_sat.json"
        )
        if path.exists():
            receipt = json.loads(path.read_text())
            if receipt.get("initial_a_hex") and receipt.get("initial_b_hex"):
                known_states.append(
                    int(receipt["initial_a_hex"], 16)
                    | (int(receipt["initial_b_hex"], 16) << 64)
                )

    rng = random.Random(SEED)
    results = []
    for offset in OFFSETS:
        nodes = offset + WINDOW
        word_mask = (1 << (64 * nodes)) - 1
        rows = [row & word_mask for row in full_rows]
        first_mask = (1 << (64 * WINDOW)) - 1
        second_mask = first_mask << (64 * offset)
        starts = [1 << bit for bit in range(128)] + known_states
        starts.extend(rng.getrandbits(128) or 1 for _ in range(RESTARTS))
        best_score = (1 << 30, 1 << 30)
        best_state = 0
        best_word = 0

        for initial in starts:
            state = initial
            word = word_from_state(state, rows)
            current = score(word, first_mask, second_mask, offset)
            while True:
                next_score = current
                next_bit = -1
                for bit, row in enumerate(rows):
                    if state == (1 << bit):
                        continue
                    candidate = score(
                        word ^ row,
                        first_mask,
                        second_mask,
                        offset,
                    )
                    if candidate < next_score:
                        next_score = candidate
                        next_bit = bit
                if next_bit < 0:
                    break
                state ^= 1 << next_bit
                word ^= rows[next_bit]
                current = next_score
            if current < best_score:
                best_score = current
                best_state = state
                best_word = word

        node_weights = [
            ((best_word >> (64 * node)) & MASK64).bit_count()
            for node in range(nodes)
        ]
        first_weight = sum(node_weights[:WINDOW])
        second_weight = sum(node_weights[offset : offset + WINDOW])
        result = {
            "offset": offset,
            "first_window_weight": first_weight,
            "second_window_weight": second_weight,
            "maximum_window_weight": max(first_weight, second_weight),
            "both_at_most_threshold": max(first_weight, second_weight) <= THRESHOLD,
            "initial_a_hex": f"0x{best_state & MASK64:016x}",
            "initial_b_hex": f"0x{best_state >> 64:016x}",
            "node_weights": node_weights,
        }
        results.append(result)
        print(
            f"offset={offset} first={first_weight} second={second_weight} "
            f"maximum={result['maximum_window_weight']}"
        )

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-repeated-low-search-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "DIAGNOSTIC_LOCAL_SEARCH_WITH_EXACT_WITNESSES",
        "window_nodes": WINDOW,
        "threshold": THRESHOLD,
        "seed_hex": f"0x{SEED:016x}",
        "random_restarts_per_offset": RESTARTS,
        "offsets": results,
        "scope_limitation": (
            "Every displayed trajectory replays exactly. Local descent does not "
            "exclude better trajectories or prove an UNSAT return distance."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print("status=DIAGNOSTIC_REPEATED_LOW_WINDOW_SEARCH")


if __name__ == "__main__":
    main()
