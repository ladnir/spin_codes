#!/usr/bin/env python3
"""Audit the bounded Goal 05 minimum-distance investigation."""

from __future__ import annotations

import hashlib
import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
PRIMARY_SOURCE = ROOT / "scripts" / "solve_riffle_packetmul_2lap_g4_goal05_native_xor.py"
INDEPENDENT_SOURCE = ROOT / "scripts" / "solve_riffle_packetmul_2lap_g4_goal05_pure_cnf.py"
PRIMARY = CANDIDATE / "receipts" / "goal05_distance_primary.json"
INDEPENDENT = CANDIDATE / "receipts" / "goal05_distance_independent.json"
GOAL04_PRIMARY = CANDIDATE / "receipts" / "goal04_four_node_primary.json"
GOAL04_AUDIT = CANDIDATE / "receipts" / "goal04_boundary_prefix_audit.json"
OUTPUT = CANDIDATE / "receipts" / "goal05_distance_audit.json"
GENERATOR = 0xF4845518B9582A1F
MASK64 = (1 << 64) - 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        low = value & -value
        result ^= columns[low.bit_length() - 1]
        value ^= low
    return result


def inverse_columns(columns: tuple[int, ...]) -> tuple[int, ...]:
    rows = []
    for output_bit in range(64):
        left = sum(
            ((columns[input_bit] >> output_bit) & 1) << input_bit
            for input_bit in range(64)
        )
        rows.append(left | (1 << (64 + output_bit)))
    for column in range(64):
        pivot = next(row for row in range(column, 64) if (rows[row] >> column) & 1)
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for row in range(64):
            if row != column and ((rows[row] >> column) & 1):
                rows[row] ^= rows[column]
    inverse_rows = tuple(row >> 64 for row in rows)
    return tuple(
        sum(
            ((inverse_rows[input_bit] >> output_bit) & 1) << input_bit
            for input_bit in range(64)
        )
        for output_bit in range(64)
    )


def state_columns() -> tuple[int, ...]:
    rows = tuple((GENERATOR << row) | (1 << 127) for row in range(64))
    left = tuple(row & MASK64 for row in rows)
    right = tuple((row >> 64) & MASK64 for row in rows)
    inverse = inverse_columns(left)
    if any(apply_columns(left, inverse[column]) != 1 << column for column in range(64)):
        raise RuntimeError("Goal 05 audit: systematic reconstruction failed")
    return tuple(apply_columns(right, inverse[column]) for column in range(64))


def accumulate(value: int) -> int:
    value ^= (value << 1) & MASK64
    value ^= (value << 2) & MASK64
    value ^= (value << 4) & MASK64
    value ^= (value << 8) & MASK64
    value ^= (value << 16) & MASK64
    value ^= (value << 32) & MASK64
    return value


def main() -> None:
    primary = load(PRIMARY)
    independent = load(INDEPENDENT)
    goal04_primary = load(GOAL04_PRIMARY)
    goal04_audit = load(GOAL04_AUDIT)

    if primary["source_sha256"] != digest(PRIMARY_SOURCE):
        raise RuntimeError("Goal 05 audit: primary source changed")
    if independent["source_sha256"] != digest(INDEPENDENT_SOURCE):
        raise RuntimeError("Goal 05 audit: independent source changed")
    expected_goal04_hash = digest(GOAL04_PRIMARY)
    if primary["goal04_primary_sha256"] != expected_goal04_hash:
        raise RuntimeError("Goal 05 audit: primary Goal 04 input changed")
    if independent["goal04_primary_sha256"] != expected_goal04_hash:
        raise RuntimeError("Goal 05 audit: independent Goal 04 input changed")
    if goal04_audit["status"] != "GOAL_04_PROVED_BOUNDARY_STRATA_CLOSED_FULL_CONSTRUCTION_OPEN":
        raise RuntimeError("Goal 05 audit: Goal 04 audit is not current")
    if primary["result"] != "UNKNOWN" or independent["result"] != "UNKNOWN":
        raise RuntimeError("Goal 05 audit: expected two bounded UNKNOWN results")
    if independent["solver_statistics"]["conflicts"] != independent["conflict_budget"]:
        raise RuntimeError("Goal 05 audit: independent conflict budget was not exhausted exactly")

    primary_witness = primary["authenticated_upper_bound_witness"]
    independent_witness = independent["authenticated_upper_bound_witness"]
    primary_states = [int(value, 16) for value in primary_witness["state_hex"]]
    independent_states = [int(value, 16) for value in independent_witness["block_hex"]]
    goal04_witness = goal04_primary["minimum_witness"]
    expected_states = [
        int(goal04_witness[key], 16)
        for key in ("left2_hex", "left1_hex", "anchor_hex", "right1_hex")
    ]
    if primary_states != expected_states or independent_states != expected_states:
        raise RuntimeError("Goal 05 audit: weight-42 witnesses disagree")

    columns = state_columns()

    def step(value: int) -> int:
        return accumulate(apply_columns(columns, value))

    replay = [expected_states[0]]
    for _ in range(3):
        replay.append(step(replay[-1]))
    weights = [value.bit_count() for value in replay]
    if replay != expected_states or sum(weights) != 42:
        raise RuntimeError("Goal 05 audit: recurrence replay failed")

    anchors_weight_9 = comb(64, 9)
    anchors_weight_10 = comb(64, 10)
    cutoff = goal04_audit["cutoff_improvement"]
    ledger = goal04_audit["support33_ledger_update"]
    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal05-distance-audit-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "evidence_label": "BOUNDED_EXACT_SOLVER_OBSTRUCTION",
        "source_sha256": digest(Path(__file__).resolve()),
        "authenticated_artifacts_sha256": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                PRIMARY_SOURCE,
                INDEPENDENT_SOURCE,
                PRIMARY,
                INDEPENDENT,
                GOAL04_PRIMARY,
                GOAL04_AUDIT,
            )
        },
        "code": {
            "definition": "C4={(o,U(o),U^2(o),U^3(o)):o in F_2^64}",
            "length": 256,
            "dimension": 64,
            "known_distance_interval": [36, 42],
            "exact_distance_resolved": False,
        },
        "authenticated_upper_bound": {
            "weight": 42,
            "block_hex": [hex(value) for value in replay],
            "block_weights": weights,
            "independent_recurrence_replay": True,
        },
        "unresolved_exact_decision": {
            "predicate": "exists nonzero o with wt(o,U(o),U^2(o),U^3(o)) <= 41",
            "primary_native_xor_result": primary["result"],
            "primary_time_budget_seconds": primary["timeout_seconds"],
            "independent_pure_cnf_result": independent["result"],
            "independent_conflict_budget": independent["conflict_budget"],
            "independent_conflicts_used": independent["solver_statistics"]["conflicts"],
        },
        "search_scale_diagnostic": {
            "candidate_block_weight_range_under_unresolved_predicate": [9, 14],
            "reason": (
                "Goal 04 excludes a block of weight at most 8 in a total-weight-at-most-41 "
                "word; four blocks of weight at least 9 force every block to have weight at most 14."
            ),
            "information_set_weight_9_anchors": anchors_weight_9,
            "information_set_weight_10_anchors": anchors_weight_10,
            "anchors_per_information_set": anchors_weight_9 + anchors_weight_10,
            "naive_four_information_set_anchor_visits": 4 * (anchors_weight_9 + anchors_weight_10),
        },
        "construction_ledger_effect": {
            "cutoff_improved_by_goal05": False,
            "retained_minimum_zero_prefix_nodes": cutoff["new_minimum_zero_prefix_nodes"],
            "retained_closed_partial_log2_interval": ledger[
                "retained_goal03_closed_partial_upper_bound_log2_interval"
            ],
            "retained_remaining_budget_log2_interval": ledger[
                "retained_remaining_numerical_budget_log2_interval"
            ],
            "conditional_independence_claim_added": False,
        },
        "status": "GOAL_05_BOUNDED_SOLVER_OBSTRUCTION_EXACT_DISTANCE_OPEN",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    print("known_distance_interval=[36,42]")
    print(f"weight_42_replay={weights}")
    print(f"anchors_per_information_set={anchors_weight_9 + anchors_weight_10}")
    print("status=GOAL_05_BOUNDED_SOLVER_OBSTRUCTION_EXACT_DISTANCE_OPEN")


if __name__ == "__main__":
    main()
