#!/usr/bin/env python3
"""Search authenticated support-33 profiles in three consecutive inner nodes."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
from collections import Counter
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import accumulate, build_apply


ROOT = Path(__file__).resolve().parents[1]
OUTER_RECEIPTS = ROOT / "explorations" / "riffle_dp_g4_g1_support33_refutation.json"
DEFAULT_OUTPUT = ROOT / "explorations" / "riffle_dp_g4_g2_support33_cluster_probe.json"
DISTANCE_THRESHOLD = 188_766
CLUSTER_NODES = 3
PACKET_SLOTS = 16
CELL_COUNT = CLUSTER_NODES * PACKET_SLOTS


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticated_profiles() -> list[dict]:
    outer = json.loads(OUTER_RECEIPTS.read_text())
    profiles: dict[tuple[int, ...], list[str]] = {}
    for receipt in outer["receipts"]:
        values = []
        for local_word in receipt["local_words"]:
            word = int(local_word["codeword_hex"], 16)
            values.extend(
                (word >> (4 * slot)) & 0xF
                for slot in range(32)
                if ((word >> (4 * slot)) & 0xF) != 0
            )
        key = tuple(sorted(values))
        if len(key) != 33:
            raise RuntimeError("cluster probe: authenticated profile is not support 33")
        profiles.setdefault(key, []).append(receipt["family_id"])
    return [
        {"values": key, "family_ids": family_ids}
        for key, family_ids in profiles.items()
    ]


def contribution_words(length: int) -> np.ndarray:
    apply_p = build_apply(systematic_state_columns())
    words = np.zeros((CELL_COUNT, 16, length), dtype=np.uint64)
    for cell in range(CELL_COUNT):
        node, slot = divmod(cell, PACKET_SLOTS)
        for bit in range(4):
            state = 0
            drive = 1 << (4 * slot + bit)
            outputs = words[cell, 1 << bit]
            for index in range(length):
                value = accumulate(state ^ drive) if index == node else accumulate(state)
                outputs[index] = value
                state = apply_p(value)
        for value in range(1, 16):
            bits = [bit for bit in range(4) if (value >> bit) & 1]
            words[cell, value] = np.bitwise_xor.reduce(
                words[cell, [1 << bit for bit in bits]], axis=0
            )
    return words


def word_for_assignment(words: np.ndarray, assignment: list[int]) -> np.ndarray:
    selected = [words[cell, value] for cell, value in enumerate(assignment) if value]
    return np.bitwise_xor.reduce(selected, axis=0)


def weight(word: np.ndarray) -> int:
    return int(np.bitwise_count(word).sum())


def anneal_profile(
    words: np.ndarray,
    values: tuple[int, ...],
    restarts: int,
    steps: int,
    rng: random.Random,
) -> tuple[int, list[int]]:
    base = list(values) + [0] * (CELL_COUNT - len(values))
    best_weight = 64 * words.shape[2] + 1
    best_assignment: list[int] = []
    for _ in range(restarts):
        assignment = base.copy()
        rng.shuffle(assignment)
        current_word = word_for_assignment(words, assignment)
        current_weight = weight(current_word)
        if current_weight < best_weight:
            best_weight = current_weight
            best_assignment = assignment.copy()
        initial_temperature = 96.0
        for step_index in range(steps):
            left = rng.randrange(CELL_COUNT)
            right = rng.randrange(CELL_COUNT - 1)
            if right >= left:
                right += 1
            left_value = assignment[left]
            right_value = assignment[right]
            if left_value == right_value:
                continue
            delta = words[left, left_value ^ right_value] ^ words[right, left_value ^ right_value]
            candidate_weight = weight(current_word ^ delta)
            temperature = initial_temperature * (1.0 - step_index / steps) + 0.25
            accept = candidate_weight <= current_weight
            if not accept:
                accept = rng.random() < math.exp((current_weight - candidate_weight) / temperature)
            if accept:
                current_word ^= delta
                assignment[left], assignment[right] = right_value, left_value
                current_weight = candidate_weight
                if current_weight < best_weight:
                    best_weight = current_weight
                    best_assignment = assignment.copy()
    return best_weight, best_assignment


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, default=5_956)
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--steps", type=int, default=5_000)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.length < CLUSTER_NODES or args.restarts <= 0 or args.steps <= 0:
        raise SystemExit("cluster probe: invalid parameters")

    profiles = authenticated_profiles()
    words = contribution_words(args.length)
    rng = random.Random(args.seed)
    rows = []
    for profile_index, profile in enumerate(profiles):
        best_weight, assignment = anneal_profile(
            words,
            profile["values"],
            args.restarts,
            args.steps,
            rng,
        )
        replay = word_for_assignment(words, assignment)
        if weight(replay) != best_weight:
            raise RuntimeError("cluster probe: best assignment replay failed")
        if Counter(assignment) != Counter(profile["values"] + (0,) * 15):
            raise RuntimeError("cluster probe: assignment changed the value multiset")
        rows.append(
            {
                "witness_evidence_label": "EXACT_REPLAY",
                "profile_index": profile_index,
                "family_ids": profile["family_ids"],
                "value_histogram": dict(sorted(Counter(profile["values"]).items())),
                "best_weight": best_weight,
                "margin_above_distance": best_weight - DISTANCE_THRESHOLD,
                "valid_terminal_translation_count_lower": (
                    args.length - CLUSTER_NODES + 1
                    if best_weight <= DISTANCE_THRESHOLD
                    else 0
                ),
                "node_words_hex": [
                    hex(sum(assignment[16 * node + slot] << (4 * slot) for slot in range(16)))
                    for node in range(CLUSTER_NODES)
                ],
            }
        )
        print(
            f"profile={profile_index} best_weight={best_weight} "
            f"margin={best_weight-DISTANCE_THRESHOLD}"
        )

    best = min(rows, key=lambda row: row["best_weight"])
    payload = {
        "schema": "riffle-dp-g4-g2-support33-cluster-probe-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "outer_receipt_sha256": digest(OUTER_RECEIPTS),
        "distance_threshold": DISTANCE_THRESHOLD,
        "window_length": args.length,
        "cluster_nodes": CLUSTER_NODES,
        "restarts_per_profile": args.restarts,
        "steps_per_restart": args.steps,
        "seed": args.seed,
        "unique_profile_count": len(profiles),
        "best_row": best,
        "profile_rows": rows,
        "scope_limitation": (
            "Each replayed placement and its reported weight are exact. Simulated "
            "annealing does not prove optimality. The search covers only three "
            "consecutive active nodes. A sub-threshold witness remains sub-threshold "
            "under every later terminal translation because translation truncates "
            "a suffix of the same output trajectory."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"best_profile={best['profile_index']} best_weight={best['best_weight']} "
        f"margin={best['margin_above_distance']}"
    )
    print(f"output={args.output}")
    print("status=DIAGNOSTIC_SUPPORT33_CLUSTER_SEARCH")


if __name__ == "__main__":
    main()
