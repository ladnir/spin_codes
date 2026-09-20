#!/usr/bin/env python3
"""Diagnostic minimum-distance search for fixed autonomous windows."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections import Counter
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import (
    apply,
    inverse_columns,
    systematic_state_columns,
)
from analyze_riffle_dp_g4_autonomous_map import accumulate, build_apply


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "explorations" / "riffle_dp_g4_g2_observability_probe.json"
OUTER_RECEIPTS = ROOT / "explorations" / "riffle_dp_g4_g1_support33_refutation.json"
DISTANCE_THRESHOLD = 188_766
EXACT_WITNESS_STATE = 0x6331ABF1B618EFAB
EXACT_WITNESS_LENGTH = 5_956
EXACT_WITNESS_WEIGHT = 188_730


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trajectory_columns(maximum_length: int) -> np.ndarray:
    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    columns = np.empty((64, maximum_length), dtype=np.uint64)
    for bit in range(64):
        state = 1 << bit
        for index in range(maximum_length):
            columns[bit, index] = accumulate(state)
            state = step(state)
    return columns


def exact_weight(columns: np.ndarray, state: int, length: int) -> int:
    bits = [bit for bit in range(64) if (state >> bit) & 1]
    if not bits:
        return 0
    word = np.bitwise_xor.reduce(columns[bits, :length], axis=0)
    return int(np.bitwise_count(word).sum())


def hill_climb(
    columns: np.ndarray,
    length: int,
    restarts: int,
    seed: int,
) -> tuple[int, int]:
    rng = random.Random(seed)
    truncated = columns[:, :length]
    best_weight = 64 * length + 1
    best_state = 0
    for _ in range(restarts):
        state = rng.getrandbits(64) or 1
        bits = [bit for bit in range(64) if (state >> bit) & 1]
        word = np.bitwise_xor.reduce(truncated[bits], axis=0)
        weight = int(np.bitwise_count(word).sum())
        while True:
            neighbors = [
                int(np.bitwise_count(word ^ truncated[bit]).sum())
                for bit in range(64)
            ]
            best_bit = min(range(64), key=neighbors.__getitem__)
            if neighbors[best_bit] >= weight:
                break
            word ^= truncated[best_bit]
            state ^= 1 << best_bit
            weight = neighbors[best_bit]
        if weight < best_weight:
            best_weight = weight
            best_state = state
    return best_weight, best_state


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lengths", default="5900,5920,5940,5960,5980,6000")
    parser.add_argument("--restarts", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    lengths = tuple(int(value) for value in args.lengths.split(","))
    if not lengths or min(lengths) <= 0 or args.restarts <= 0:
        raise SystemExit("observability probe: invalid parameters")
    maximum_length = max(max(lengths), EXACT_WITNESS_LENGTH)
    columns = trajectory_columns(maximum_length)
    witness_weight = exact_weight(
        columns,
        EXACT_WITNESS_STATE,
        EXACT_WITNESS_LENGTH,
    )
    if witness_weight != EXACT_WITNESS_WEIGHT:
        raise RuntimeError("observability probe: exact witness replay failed")
    p_columns = systematic_state_columns()
    t_columns = tuple(apply(p_columns, accumulate(1 << bit)) for bit in range(64))
    initial_drive = apply(inverse_columns(t_columns), EXACT_WITNESS_STATE)
    initial_values = [
        (initial_drive >> (4 * slot)) & 0xF
        for slot in range(16)
        if ((initial_drive >> (4 * slot)) & 0xF) != 0
    ]
    required = Counter(initial_values)
    outer = json.loads(OUTER_RECEIPTS.read_text())
    compatible_families = []
    for receipt in outer["receipts"]:
        available = Counter()
        for local_word in receipt["local_words"]:
            word = int(local_word["codeword_hex"], 16)
            available.update(
                (word >> (4 * packet)) & 0xF
                for packet in range(32)
                if ((word >> (4 * packet)) & 0xF) != 0
            )
        if all(available[value] >= count for value, count in required.items()):
            compatible_families.append(receipt["family_id"])
    if compatible_families:
        raise RuntimeError("observability probe: witness unexpectedly matches support 33")
    rows = []
    for offset, length in enumerate(lengths):
        weight, state = hill_climb(
            columns,
            length,
            args.restarts,
            args.seed + offset,
        )
        if exact_weight(columns, state, length) != weight:
            raise RuntimeError("observability probe: hill witness replay failed")
        rows.append(
            {
                "length": length,
                "best_weight": weight,
                "margin_above_distance": weight - DISTANCE_THRESHOLD,
                "state_hex": hex(state),
            }
        )
    payload = {
        "schema": "riffle-dp-g4-g2-observability-probe-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "distance_threshold": DISTANCE_THRESHOLD,
        "restarts_per_length": args.restarts,
        "seed": args.seed,
        "exact_witness": {
            "evidence_label": "EXACT",
            "length": EXACT_WITNESS_LENGTH,
            "weight": witness_weight,
            "state_hex": hex(EXACT_WITNESS_STATE),
            "one_node_preimage_hex": hex(initial_drive),
            "one_node_preimage_nibble_support": len(initial_values),
            "one_node_preimage_values": initial_values,
            "compatible_authenticated_support33_families": compatible_families,
            "interpretation": (
                "A universal autonomous-window bound cannot force weight above d "
                "at every length at most 5939."
            ),
        },
        "hill_climb_rows": rows,
        "scope_limitation": (
            "Hill climbing supplies upper bounds on window minimum distance, not "
            "lower bounds. The exact witness is not known to be reachable from an "
            "authenticated sparse outer word."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(
        f"exact_witness_length={EXACT_WITNESS_LENGTH} "
        f"weight={witness_weight} margin={witness_weight-DISTANCE_THRESHOLD}"
    )
    for row in rows:
        print(
            f"length={row['length']} best_weight={row['best_weight']} "
            f"margin={row['margin_above_distance']} state={row['state_hex']}"
        )
    print(f"output={args.output}")
    print("status=DIAGNOSTIC_G2_OBSERVABILITY_DISTANCE")


if __name__ == "__main__":
    main()
