#!/usr/bin/env python3
"""Probe a Cauchy--Schwarz route to the full-orbit character caps."""

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


ROOT = Path(__file__).resolve().parents[1]
FULL_BIAS_RECEIPT = ROOT / "explorations" / "riffle_dp_g4_g2_full_character_bias_probe.json"
COMPONENT_MOMENT_RECEIPT = (
    ROOT / "explorations" / "riffle_dp_g4_component_second_moment_certificate.json"
)
DEFAULT_OUTPUT = ROOT / "explorations" / "riffle_dp_g4_global_second_moment_probe.json"
LOCAL_COUNTEREXAMPLE = 0x852AC8FCC6B27C67


def parity_for_character(rows: np.ndarray, character: int) -> np.ndarray:
    active = [bit for bit in range(64) if (character >> bit) & 1]
    if not active:
        return np.zeros(rows.shape[1], dtype=np.bool_)
    return np.logical_xor.reduce(rows[active], axis=0)


def observation_kernel(states: np.ndarray) -> tuple[int, ...]:
    """Return a basis for characters orthogonal to every observation state."""
    columns = tuple(
        sum(
            (int((state >> np.uint64(bit)) & np.uint64(1)) << coordinate)
            for coordinate, state in enumerate(states)
        )
        for bit in range(64)
    )
    pivots: dict[int, tuple[int, int]] = {}
    kernel = []
    for bit, column in enumerate(columns):
        value = column
        representation = 1 << bit
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (value, representation)
                break
            pivot_value, pivot_representation = pivots[pivot]
            value ^= pivot_value
            representation ^= pivot_representation
        if not value:
            kernel.append(representation)
    return tuple(kernel)


def xor_selected(values: list[int], selection: int) -> int:
    result = 0
    for index, value in enumerate(values):
        if (selection >> index) & 1:
            result ^= value
    return result


def metrics(rows: np.ndarray, character: int) -> dict[str, int | float | str]:
    parity = parity_for_character(rows, character)
    signs = 1 - 2 * parity.astype(np.int16)
    node_sums = signs.reshape(PACKET_SLOTS, INNER_NODES).sum(axis=0)
    total = int(node_sums.sum())
    square_sum = int(np.dot(node_sums.astype(np.int64), node_sums.astype(np.int64)))
    prefix = np.cumsum(node_sums, dtype=np.int64)
    block_sums = node_sums.reshape(-1, 12).sum(axis=1)
    return {
        "character_hex": hex(character),
        "global_character_sum": total,
        "absolute_global_bias": abs(total) / SAMPLE_SPACE,
        "node_second_moment": square_sum / INNER_NODES,
        "cauchy_bias_bound": math.sqrt(square_sum / INNER_NODES) / PACKET_SLOTS,
        "maximum_absolute_prefix_sum": int(np.abs(prefix).max()),
        "maximum_absolute_12_node_sum": int(np.abs(block_sums).max()),
        "blocks_over_non15_local_cap": int(np.count_nonzero(np.abs(block_sums) > 96)),
        "blocks_over_value15_local_cap": int(np.count_nonzero(np.abs(block_sums) > 120)),
    }


def second_moment_objective(rows: np.ndarray, character: int) -> int:
    node_count = rows.shape[1] // PACKET_SLOTS
    parity = parity_for_character(rows, character)
    signs = 1 - 2 * parity.astype(np.int16)
    node_sums = signs.reshape(PACKET_SLOTS, node_count).sum(axis=0).astype(np.int64)
    return int(np.dot(node_sums, node_sums))


def hill_climb_second_moment(rows: np.ndarray, start: int) -> tuple[int, int]:
    if rows.shape[1] % PACKET_SLOTS:
        raise ValueError("second moment: sample count is not slot aligned")
    node_count = rows.shape[1] // PACKET_SLOTS
    character = start or 1
    parity = parity_for_character(rows, character).reshape(PACKET_SLOTS, node_count)
    signs = 1 - 2 * parity.astype(np.int16)
    node_sums = signs.sum(axis=0)
    objective = int(np.dot(node_sums.astype(np.int64), node_sums.astype(np.int64)))
    row_blocks = rows.reshape(64, PACKET_SLOTS, node_count)
    while True:
        best_objective = objective
        best_bit = None
        best_node_sums = None
        for bit in range(64):
            affected_sum = (signs * row_blocks[bit]).sum(axis=0)
            candidate_node_sums = node_sums - 2 * affected_sum
            candidate_objective = int(
                np.dot(
                    candidate_node_sums.astype(np.int64),
                    candidate_node_sums.astype(np.int64),
                )
            )
            if candidate_objective > best_objective:
                best_objective = candidate_objective
                best_bit = bit
                best_node_sums = candidate_node_sums
        if best_bit is None:
            return objective, character
        mask = row_blocks[best_bit]
        signs[mask] *= -1
        parity[mask] ^= True
        node_sums = best_node_sums
        objective = best_objective
        character ^= 1 << best_bit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-restarts", type=int, default=20)
    parser.add_argument("--local-random-restarts", type=int, default=100)
    parser.add_argument("--local-windows", default="1,2,3,4,5,6,8,12,24")
    parser.add_argument("--seed", type=int, default=20260821)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.random_restarts < 0 or args.local_random_restarts < 0:
        raise SystemExit("global second moment: negative restart count")
    local_windows = tuple(int(item) for item in args.local_windows.split(",") if item)
    if not local_windows or min(local_windows) < 1 or max(local_windows) > INNER_NODES:
        raise SystemExit("global second moment: invalid local window")

    bias_receipt = json.loads(FULL_BIAS_RECEIPT.read_text())
    bias_characters = {
        row["packet_value"]: int(row["character_hex"], 16)
        for row in bias_receipt["value_rows"]
    }
    component_receipt = json.loads(COMPONENT_MOMENT_RECEIPT.read_text())
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
        states = contribution_states(step, packet_value)
        rows = bit_rows(states)
        seeds = [
            ("full_bias_probe", bias_characters[packet_value]),
            ("local_counterexample", LOCAL_COUNTEREXAMPLE),
        ]
        seeds.extend(component_characters[packet_value])
        pure_characters = [character for _source, character in component_characters[packet_value]]
        structured_characters = tuple(
            xor_selected(pure_characters, subset)
            for subset in range(1, 1 << len(pure_characters))
        )
        structured_best = max(
            structured_characters,
            key=lambda character: second_moment_objective(rows, character),
        )
        seeds.append(("component_maximizer_subset", structured_best))
        seeds.extend(
            ("random", rng.getrandbits(64) or 1)
            for _ in range(args.random_restarts)
        )
        best = None
        for source, seed in seeds:
            objective, character = hill_climb_second_moment(rows, seed)
            candidate = (objective, character, source)
            if best is None or candidate[0] > best[0]:
                best = candidate
        assert best is not None
        threshold = 100 if packet_value == 15 else 64
        best_metrics = metrics(rows, best[1])
        refutes = best_metrics["node_second_moment"] > threshold
        any_refutation |= refutes
        local_rows = []
        row_cube = rows.reshape(64, PACKET_SLOTS, INNER_NODES)
        for window in local_windows:
            local_bits = row_cube[:, :, :window].reshape(64, PACKET_SLOTS * window)
            local_seeds = [
                ("full_bias_probe", bias_characters[packet_value]),
                ("local_counterexample", LOCAL_COUNTEREXAMPLE),
            ]
            local_seeds.extend(component_characters[packet_value])
            local_states = states.reshape(PACKET_SLOTS, INNER_NODES)[:, :window].reshape(
                PACKET_SLOTS * window
            )
            local_seeds.extend(
                ("observation_kernel", character)
                for character in observation_kernel(local_states)
            )
            local_seeds.extend(
                ("random", rng.getrandbits(64) or 1)
                for _ in range(args.local_random_restarts)
            )
            local_best = None
            for source, seed in local_seeds:
                objective, character = hill_climb_second_moment(local_bits, seed)
                candidate = (objective, character, source)
                if local_best is None or candidate[0] > local_best[0]:
                    local_best = candidate
            assert local_best is not None
            local_mean_square = local_best[0] / window
            local_rows.append(
                {
                    "window_nodes": window,
                    "maximum_found_mean_square": local_mean_square,
                    "character_hex": hex(local_best[1]),
                    "seed_source": local_best[2],
                    "refutes_local_cap": local_mean_square > threshold,
                }
            )
        value_rows.append(
            {
                "packet_value": packet_value,
                "proposed_second_moment_cap": threshold,
                "refutes_proposed_cap": refutes,
                "seed_source": best[2],
                **best_metrics,
                "full_bias_extremizer_metrics": metrics(
                    rows, bias_characters[packet_value]
                ),
                "local_counterexample_metrics": metrics(rows, LOCAL_COUNTEREXAMPLE),
                "initial_window_rows": local_rows,
            }
        )
        print(
            f"value={packet_value} max_mean_square={best_metrics['node_second_moment']:.9f} "
            f"cap={threshold} character={best_metrics['character_hex']} "
            f"refutes={refutes}"
        )

    payload = {
        "schema": "riffle-dp-g4-global-second-moment-probe-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": (
            "DIAGNOSTIC_COUNTEREXAMPLE" if any_refutation else "DIAGNOSTIC"
        ),
        "node_count": INNER_NODES,
        "packet_slots": PACKET_SLOTS,
        "random_restarts_per_value": args.random_restarts,
        "local_random_restarts_per_value_and_window": args.local_random_restarts,
        "local_windows": list(local_windows),
        "seed": args.seed,
        "value_rows": value_rows,
        "scope_limitation": (
            "Every reported metric is exact for its character. Coordinate ascent "
            "does not prove a global maximum. A reported violation does refute the "
            "proposed second-moment cap."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={args.output}")
    print(f"status={payload['evidence_label']}")


if __name__ == "__main__":
    main()
