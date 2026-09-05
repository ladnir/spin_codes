#!/usr/bin/env python3
"""Probe a 12-node block-energy route to the full character caps."""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

import numpy as np

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import accumulate, build_apply
from probe_riffle_dp_g4_component_mixing import (
    INNER_NODES,
    PACKET_SLOTS,
    SAMPLE_SPACE,
    contribution_states,
)
from probe_riffle_dp_g4_full_character_bias import bit_rows
from probe_riffle_dp_g4_global_second_moment import (
    COMPONENT_MOMENT_RECEIPT,
    FULL_BIAS_RECEIPT,
    LOCAL_COUNTEREXAMPLE,
    parity_for_character,
    xor_selected,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "explorations" / "riffle_dp_g4_global_block_moment_probe.json"
FULL_LOCAL_DISTANCE_RECEIPT = (
    ROOT / "explorations" / "riffle_dp_2lap_g4_full_local_distance_probe.json"
)
BLOCK_NODES = 12
BLOCK_COUNT = INNER_NODES // BLOCK_NODES
BLOCK_CELLS = BLOCK_NODES * PACKET_SLOTS


def sums_and_objective(rows: np.ndarray, character: int) -> tuple[np.ndarray, int]:
    parity = parity_for_character(rows, character)
    signs = 1 - 2 * parity.astype(np.int16)
    block_sums = signs.reshape(PACKET_SLOTS, BLOCK_COUNT, BLOCK_NODES).sum(axis=(0, 2))
    wide = block_sums.astype(np.int64)
    return block_sums, int(np.dot(wide, wide))


def hill_climb(rows: np.ndarray, start: int) -> tuple[int, int]:
    character = start or 1
    parity = parity_for_character(rows, character).reshape(
        PACKET_SLOTS, INNER_NODES
    )
    signs = 1 - 2 * parity.astype(np.int16)
    block_sums = signs.reshape(PACKET_SLOTS, BLOCK_COUNT, BLOCK_NODES).sum(axis=(0, 2))
    wide = block_sums.astype(np.int64)
    objective = int(np.dot(wide, wide))
    row_blocks = rows.reshape(64, PACKET_SLOTS, BLOCK_COUNT, BLOCK_NODES)
    while True:
        best_objective = objective
        best_bit = None
        best_sums = None
        for bit in range(64):
            affected = (signs.reshape(PACKET_SLOTS, BLOCK_COUNT, BLOCK_NODES) * row_blocks[bit]).sum(
                axis=(0, 2)
            )
            candidate = block_sums - 2 * affected
            candidate_wide = candidate.astype(np.int64)
            candidate_objective = int(np.dot(candidate_wide, candidate_wide))
            if candidate_objective > best_objective:
                best_objective = candidate_objective
                best_bit = bit
                best_sums = candidate
        if best_bit is None:
            return objective, character
        mask = row_blocks[best_bit].reshape(PACKET_SLOTS, INNER_NODES)
        signs[mask] *= -1
        parity[mask] ^= True
        block_sums = best_sums
        objective = best_objective
        character ^= 1 << best_bit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-restarts", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260821)
    args = parser.parse_args()
    if args.random_restarts < 0:
        raise SystemExit("global block moment: negative restart count")

    bias_receipt = json.loads(FULL_BIAS_RECEIPT.read_text())
    bias_characters = {
        row["packet_value"]: int(row["character_hex"], 16)
        for row in bias_receipt["value_rows"]
    }
    component_receipt = json.loads(COMPONENT_MOMENT_RECEIPT.read_text())
    local_distance_receipt = json.loads(FULL_LOCAL_DISTANCE_RECEIPT.read_text())
    local_distance_characters = tuple(
        sorted({int(row["character_hex"], 16) for row in local_distance_receipt["cases"]})
    )
    component_characters = {
        value: [
            (
                f"component_{component['component_degree']}",
                int(
                    next(
                        row
                        for row in component["value_rows"]
                        if row["packet_value"] == value
                    )["maximizing_character_hex"],
                    16,
                ),
            )
            for component in component_receipt["component_rows"]
        ]
        for value in range(1, 16)
    }
    apply_p = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_p(accumulate(state))

    rng = random.Random(args.seed)
    value_rows = []
    any_refutation = False
    for packet_value in range(1, 16):
        rows = bit_rows(contribution_states(step, packet_value))
        pure = [character for _source, character in component_characters[packet_value]]
        structured = tuple(xor_selected(pure, subset) for subset in range(1, 1 << len(pure)))
        structured_best = max(
            structured,
            key=lambda character: sums_and_objective(rows, character)[1],
        )
        local_distance_best = max(
            local_distance_characters,
            key=lambda character: sums_and_objective(rows, character)[1],
        )
        seeds = [
            ("full_bias_probe", bias_characters[packet_value]),
            ("local_counterexample", LOCAL_COUNTEREXAMPLE),
            ("component_maximizer_subset", structured_best),
            ("full_local_distance_witness", local_distance_best),
            *component_characters[packet_value],
        ]
        seeds.extend(
            ("random", rng.getrandbits(64) or 1)
            for _ in range(args.random_restarts)
        )
        best = None
        for source, seed in seeds:
            objective, character = hill_climb(rows, seed)
            candidate = (objective, character, source)
            if best is None or candidate[0] > best[0]:
                best = candidate
        assert best is not None
        block_sums, objective = sums_and_objective(rows, best[1])
        mean_square = objective / BLOCK_COUNT
        cap = 120**2 if packet_value == 15 else 96**2
        refutes = mean_square > cap
        any_refutation |= refutes
        total = int(block_sums.sum())
        row = {
            "packet_value": packet_value,
            "proposed_block_mean_square_cap": cap,
            "refutes_proposed_cap": refutes,
            "maximum_found_block_mean_square": mean_square,
            "character_hex": hex(best[1]),
            "seed_source": best[2],
            "global_character_sum": total,
            "absolute_global_bias": abs(total) / SAMPLE_SPACE,
            "cauchy_bias_bound": math.sqrt(mean_square) / BLOCK_CELLS,
            "maximum_absolute_block_sum": int(np.abs(block_sums).max()),
        }
        value_rows.append(row)
        print(
            f"value={packet_value} max_block_mean_square={mean_square:.9f} "
            f"cap={cap} character={row['character_hex']} refutes={refutes}"
        )

    payload = {
        "schema": "riffle-dp-g4-global-block-moment-probe-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": (
            "DIAGNOSTIC_COUNTEREXAMPLE" if any_refutation else "DIAGNOSTIC"
        ),
        "block_nodes": BLOCK_NODES,
        "block_count": BLOCK_COUNT,
        "block_cells": BLOCK_CELLS,
        "random_restarts_per_value": args.random_restarts,
        "seed": args.seed,
        "value_rows": value_rows,
        "scope_limitation": (
            "Every metric is exact for its reported character. The structured "
            "seed scan and coordinate ascent do not prove a global maximum. A "
            "reported violation does refute the proposed block-moment cap."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print(f"status={payload['evidence_label']}")


if __name__ == "__main__":
    main()
