#!/usr/bin/env python3
"""Deterministic local search for a counterexample to the C38 target."""

from __future__ import annotations

import hashlib
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
GOAL06 = CANDIDATE / "receipts" / "goal06_amortized_route_audit.json"
OUTPUT = CANDIDATE / "receipts" / "goal08_c38_refutation_search.json"
NODES = 38
TARGET = 228
RESTARTS = 10_000
SEED = 0x4333385245465554
MASK64 = (1 << 64) - 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def word_from_state(state: int, rows: list[int]) -> int:
    word = 0
    while state:
        bit = state & -state
        word ^= rows[bit.bit_length() - 1]
        state ^= bit
    return word


def single_descent(word: int, state: int, rows: list[int]) -> tuple[int, int]:
    weight = word.bit_count()
    while True:
        best_weight = weight
        best_bit = -1
        for bit, row in enumerate(rows):
            candidate = (word ^ row).bit_count()
            if candidate < best_weight:
                best_weight = candidate
                best_bit = bit
        if best_bit < 0:
            return word, state
        word ^= rows[best_bit]
        state ^= 1 << best_bit
        weight = best_weight


def pair_descent(word: int, state: int, rows: list[int]) -> tuple[int, int, int]:
    weight = word.bit_count()
    rounds = 0
    while True:
        best_weight = weight
        best_pair: tuple[int, int] | None = None
        for left in range(128):
            left_word = word ^ rows[left]
            for right in range(left + 1, 128):
                candidate = (left_word ^ rows[right]).bit_count()
                if candidate < best_weight:
                    best_weight = candidate
                    best_pair = (left, right)
        if best_pair is None:
            return word, state, rounds
        left, right = best_pair
        word ^= rows[left] ^ rows[right]
        state ^= (1 << left) | (1 << right)
        weight = best_weight
        rounds += 1


def main() -> None:
    parity = build_apply(systematic_state_columns())
    rows = observation_rows(NODES, parity)
    goal06 = json.loads(GOAL06.read_text())
    transient = goal06["goal05_counterexample_extension"]
    known_state = int(transient["initial_a_hex"], 16) | (
        int(transient["initial_b_hex"], 16) << 64
    )

    rng = random.Random(SEED)
    starts = [1 << bit for bit in range(128)]
    starts.append(known_state)
    starts.extend(rng.getrandbits(128) or 1 for _ in range(RESTARTS))

    best_word = (1 << (64 * NODES)) - 1
    best_state = 0
    for start in starts:
        word, state = single_descent(word_from_state(start, rows), start, rows)
        if 0 < word.bit_count() < best_word.bit_count():
            best_word = word
            best_state = state

    best_word, best_state, pair_rounds = pair_descent(best_word, best_state, rows)
    replay = word_from_state(best_state, rows)
    if replay != best_word:
        raise RuntimeError("C38 refutation search: witness replay failed")
    node_weights = [
        ((best_word >> (64 * node)) & MASK64).bit_count()
        for node in range(NODES)
    ]
    known_weights = [
        ((word_from_state(known_state, rows) >> (64 * node)) & MASK64).bit_count()
        for node in range(NODES)
    ]
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal08-c38-refutation-search-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "DIAGNOSTIC_LOCAL_SEARCH_WITH_EXACT_WITNESSES",
        "source_sha256": digest(Path(__file__).resolve()),
        "goal06_receipt_sha256": digest(GOAL06),
        "window_nodes": NODES,
        "distance_target": TARGET,
        "seed_hex": hex(SEED),
        "random_restarts": RESTARTS,
        "total_starts": len(starts),
        "pair_descent_rounds": pair_rounds,
        "known_transient_c38_weight": sum(known_weights),
        "best_weight": best_word.bit_count(),
        "counterexample_found": best_word.bit_count() < TARGET,
        "initial_a_hex": f"0x{best_state & MASK64:016x}",
        "initial_b_hex": f"0x{best_state >> 64:016x}",
        "node_weights": node_weights,
        "scope_limitation": (
            "Every displayed weight replays exactly. Random-restart coordinate descent "
            "and pair descent neither prove a lower bound nor exclude a counterexample."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"known_transient_c38_weight={sum(known_weights)}")
    print(f"best_weight={payload['best_weight']}")
    print(f"counterexample_found={payload['counterexample_found']}")
    print(f"output={OUTPUT}")
    print("status=DIAGNOSTIC_C38_REFUTATION_SEARCH")


if __name__ == "__main__":
    main()
