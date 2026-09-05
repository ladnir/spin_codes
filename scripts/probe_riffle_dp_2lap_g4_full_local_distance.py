#!/usr/bin/env python3
"""Diagnostic full-character search for the endpoint-aware local distances."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import (
    accumulate,
    build_apply,
    krylov_polynomial,
    polynomial_mod,
)
from probe_riffle_dp_g4_component_mixing import transpose_columns


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "explorations" / "riffle_dp_2lap_g4_full_local_distance_probe.json"
)
FACTORS = (
    (1, 0x3),
    (2, 0x7),
    (4, 0x13),
    (9, 0x373),
    (10, 0x519),
    (18, 0x7C9C3),
    (20, 0x1E1FFF),
)
TIGHT_DEGREE_ONE_CHARACTER = 0xA685AAC60E5ACFB3


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def observation_states(step, value: int, window: int) -> np.ndarray:
    states = np.empty(16 * window, dtype=np.uint64)
    cursor = 0
    for slot in range(16):
        state = value << (4 * slot)
        for _ in range(window):
            state = step(state)
            states[cursor] = state
            cursor += 1
    return states


def bit_rows(states: np.ndarray) -> np.ndarray:
    rows = np.empty((64, len(states)), dtype=np.bool_)
    for bit in range(64):
        rows[bit] = ((states >> np.uint64(bit)) & np.uint64(1)).astype(np.bool_)
    return rows


def exact_weight(rows: np.ndarray, character: int) -> int:
    active = [bit for bit in range(64) if (character >> bit) & 1]
    if not active:
        return 0
    return int(np.count_nonzero(np.logical_xor.reduce(rows[active], axis=0)))


def hill_climb(rows: np.ndarray, start: int) -> tuple[int, int, int]:
    character = start or 1
    active = [bit for bit in range(64) if (character >> bit) & 1]
    word = np.logical_xor.reduce(rows[active], axis=0)
    length = rows.shape[1]
    weight = int(np.count_nonzero(word))
    side = min(weight, length - weight)
    steps = 0
    while True:
        neighbor_weights = np.count_nonzero(rows != word, axis=1)
        neighbor_sides = np.minimum(neighbor_weights, length - neighbor_weights)
        if character.bit_count() == 1:
            neighbor_sides[(character & -character).bit_length() - 1] = length + 1
        bit = int(np.argmin(neighbor_sides))
        next_side = int(neighbor_sides[bit])
        if next_side >= side:
            break
        word ^= rows[bit]
        character ^= 1 << bit
        weight = int(neighbor_weights[bit])
        side = next_side
        steps += 1
    return side, weight, character


def support_data(transpose_step, character: int) -> tuple[int, list[int], int]:
    polynomial = krylov_polynomial(transpose_step, character)
    degrees = [degree for degree, factor in FACTORS if polynomial_mod(polynomial, factor) == 0]
    mask = sum(1 << index for index, (degree, _) in enumerate(FACTORS) if degree in degrees)
    return mask, degrees, polynomial


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-restarts", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260820)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.random_restarts < 1:
        raise SystemExit("full local-distance probe: restarts must be positive")

    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    transpose_step = build_apply(
        transpose_columns(tuple(step(1 << bit) for bit in range(64)))
    )
    rng = random.Random(args.seed)
    cases = [(24, value) for value in range(1, 16)] + [(12, 15)]
    case_rows = []
    for case_index, (window, value) in enumerate(cases):
        rows = bit_rows(observation_states(step, value, window))
        required = 97 if value != 15 else (72 if window == 24 else 36)
        seeds = [TIGHT_DEGREE_ONE_CHARACTER]
        seeds.extend(1 << bit for bit in range(64))
        seeds.extend(rng.getrandbits(64) or 1 for _ in range(args.random_restarts))
        best = None
        for seed_index, seed in enumerate(seeds):
            side, weight, character = hill_climb(rows, seed)
            candidate = (side, weight, character, seed_index, seed)
            if best is None or candidate[:3] < best[:3]:
                best = candidate
        assert best is not None
        side, weight, character, seed_index, seed = best
        replay_weight = exact_weight(rows, character)
        if replay_weight != weight:
            raise RuntimeError("full local-distance probe: witness replay failed")
        support_mask, support_degrees, polynomial = support_data(
            transpose_step, character
        )
        row = {
            "window_nodes": window,
            "packet_value": value,
            "code_length": 16 * window,
            "required_two_sided_distance": required,
            "best_found_two_sided_weight": side,
            "weight": weight,
            "complement_weight": 16 * window - weight,
            "violates_required_distance": side < required,
            "character_hex": hex(character),
            "component_support_mask_hex": hex(support_mask),
            "component_degrees": support_degrees,
            "component_dimension": sum(support_degrees),
            "character_minimal_polynomial_hex": hex(polynomial),
            "seed_index": seed_index,
            "seed_character_hex": hex(seed),
        }
        case_rows.append(row)
        print(
            f"window={window} value={value} best_side={side} required={required} "
            f"support={hex(support_mask)} dimension={sum(support_degrees)}",
            flush=True,
        )

    violations = [row for row in case_rows if row["violates_required_distance"]]
    payload = {
        "schema": "riffle-dp-2lap-g4-full-local-distance-probe-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_WITNESSES_FROM_DIAGNOSTIC_SEARCH",
        "source_sha256": digest(Path(__file__).resolve()),
        "random_restarts_per_case": args.random_restarts,
        "seed": args.seed,
        "case_count": len(case_rows),
        "violation_count": len(violations),
        "cases": case_rows,
        "violations": violations,
        "result": "COUNTEREXAMPLE_FOUND" if violations else "NO_COUNTEREXAMPLE_FOUND",
        "scope_limitation": (
            "Every reported character and weight is exact. Coordinate-descent "
            "search does not prove a lower bound when it finds no violation."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"violations={len(violations)}")
    print(f"output={args.output}")
    print(f"status={payload['result']}")


if __name__ == "__main__":
    main()
