#!/usr/bin/env python3
"""Independently verify the Goal 02 tail-list and ledger certificate."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
PARENT = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
PRIMARY_SOURCE = (
    ROOT / "scripts" / "certify_riffle_packetmul_wrapmul_2lap_g4_goal02_tail_lists.py"
)
PRIMARY = CANDIDATE / "receipts" / "goal02_tail_list_primary.json"
GOAL04 = PARENT / "receipts" / "goal04_boundary_prefix_audit.json"
OUTPUT = CANDIDATE / "receipts" / "goal02_tail_list_independent.json"

GENERATOR = 0xF4845518B9582A1F
MASK64 = (1 << 64) - 1
TOTAL_NODES = 32_772
TOTAL_POSITIONS = TOTAL_NODES * 16
SUPPORT = 33
WORDS = 26
NONZERO_STATES = (1 << 64) - 1
BAD_MAXIMUM = 188_765
JOHNSON_CUTOFF = 7_060
UNIQUE_CUTOFF = 9_176
EMPTY_CUTOFF = 20_976


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    for bit in range(64):
        if (value >> bit) & 1:
            result ^= columns[bit]
    return result


def inverse_columns(columns: tuple[int, ...]) -> tuple[int, ...]:
    rows = []
    for output in range(64):
        left = sum(
            ((columns[input_bit] >> output) & 1) << input_bit
            for input_bit in range(64)
        )
        rows.append(left | (1 << (64 + output)))
    for column in range(64):
        pivot = next(row for row in range(column, 64) if (rows[row] >> column) & 1)
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for row in range(64):
            if row != column and ((rows[row] >> column) & 1):
                rows[row] ^= rows[column]
    inverse_rows = tuple(row >> 64 for row in rows)
    result = tuple(
        sum(
            ((inverse_rows[input_bit] >> output_bit) & 1) << input_bit
            for input_bit in range(64)
        )
        for output_bit in range(64)
    )
    if any(apply_columns(columns, result[bit]) != 1 << bit for bit in range(64)):
        raise RuntimeError("Goal 02 independent: inverse reconstruction failed")
    return result


def systematic_right_columns() -> tuple[int, ...]:
    rows = tuple((GENERATOR << row) | (1 << 127) for row in range(64))
    left = tuple(row & MASK64 for row in rows)
    right = tuple((row >> 64) & MASK64 for row in rows)
    inverse = inverse_columns(left)
    columns = tuple(apply_columns(right, inverse[bit]) for bit in range(64))
    if any(apply_columns(left, inverse[bit]) != 1 << bit for bit in range(64)):
        raise RuntimeError("Goal 02 independent: systematic map failed")
    return columns


def accumulate_by_bits(value: int) -> int:
    running = 0
    result = 0
    for bit in range(64):
        running ^= (value >> bit) & 1
        result |= running << bit
    return result


def rank_by_rows(columns: tuple[int, ...]) -> int:
    rows = [
        sum(((columns[column] >> row) & 1) << column for column in range(64))
        for row in range(64)
    ]
    rank = 0
    for column in range(64):
        pivot = next(
            (row for row in range(rank, 64) if (rows[row] >> column) & 1),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for row in range(64):
            if row != rank and ((rows[row] >> column) & 1):
                rows[row] ^= rows[rank]
        rank += 1
    return rank


def integer_list_bound(first_node: int) -> tuple[int | None, bool, tuple[int, int, int]]:
    blocks = first_node // 4
    remaining = TOTAL_NODES // 4 - blocks
    length = 256 * remaining
    distance = 36 * remaining
    radius = BAD_MAXIMUM - 36 * blocks
    denominator = distance * length - 2 * radius * (length - radius)
    bound = None if denominator <= 0 else (distance * length) // denominator
    unique = distance > 2 * radius
    return (1 if unique else bound), unique, (length, distance, radius)


def as_fraction(record: dict) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def main() -> None:
    primary = load(PRIMARY)
    goal04 = load(GOAL04)
    if primary["source_sha256"] != digest(PRIMARY_SOURCE):
        raise RuntimeError("Goal 02 independent: primary source changed")

    parity = systematic_right_columns()
    state_update = tuple(
        apply_columns(parity, accumulate_by_bits(1 << bit)) for bit in range(64)
    )
    output_update = tuple(
        accumulate_by_bits(parity[bit]) for bit in range(64)
    )
    if rank_by_rows(parity) != 64:
        raise RuntimeError("Goal 02 independent: parity map is singular")
    if rank_by_rows(state_update) != 64 or rank_by_rows(output_update) != 64:
        raise RuntimeError("Goal 02 independent: recurrence map is singular")

    # Replay common-drive cancellation on every pair of basis directions.
    checked_pairs = 0
    for state_bit in range(64):
        state = 1 << state_bit
        for drive_bit in range(64):
            drive = 1 << drive_bit
            left_output = accumulate_by_bits(state ^ drive)
            right_output = accumulate_by_bits(drive)
            if left_output ^ right_output != accumulate_by_bits(state):
                raise RuntimeError("Goal 02 independent: emission cancellation failed")
            left_state = apply_columns(parity, left_output)
            right_state = apply_columns(parity, right_output)
            if left_state ^ right_state != state_update[state_bit]:
                raise RuntimeError("Goal 02 independent: state cancellation failed")
            checked_pairs += 1

    rows = [integer_list_bound(r) for r in range(EMPTY_CUTOFF)]
    first_johnson = next(r for r, row in enumerate(rows) if row[0] is not None)
    first_unique = next(r for r, row in enumerate(rows) if row[1])
    maximum_late = max(row[0] for row in rows[JOHNSON_CUTOFF:] if row[0] is not None)
    if (first_johnson, first_unique, maximum_late) != (
        JOHNSON_CUTOFF,
        UNIQUE_CUTOFF,
        13_749,
    ):
        raise RuntimeError("Goal 02 independent: cutoff replay failed")

    denominator = comb(TOTAL_POSITIONS, SUPPORT)
    late_placement = Fraction(
        comb(TOTAL_POSITIONS - 16 * JOHNSON_CUTOFF, SUPPORT)
        - comb(TOTAL_POSITIONS - 16 * EMPTY_CUTOFF, SUPPORT),
        denominator,
    )
    early_placement = Fraction(
        denominator - comb(TOTAL_POSITIONS - 16 * JOHNSON_CUTOFF, SUPPORT),
        denominator,
    )
    late_charge = WORDS * late_placement * Fraction(maximum_late, NONZERO_STATES)
    old_budget = as_fraction(
        goal04["support33_ledger_update"]["retained_remaining_numerical_budget"]
    )
    new_budget = old_budget - late_charge
    early_cap = new_budget * NONZERO_STATES // (WORDS * early_placement)

    primary_late = primary["late_band_ledger"]
    primary_early = primary["early_region_target"]
    if as_fraction(primary_late["placement_probability"]) != late_placement:
        raise RuntimeError("Goal 02 independent: late placement differs")
    if as_fraction(primary_late["aggregate_charge"]) != late_charge:
        raise RuntimeError("Goal 02 independent: late charge differs")
    if as_fraction(primary_late["new_remaining_budget"]) != new_budget:
        raise RuntimeError("Goal 02 independent: remaining budget differs")
    if primary_early["maximum_uniform_bad_states_per_fixed_drive"] != early_cap:
        raise RuntimeError("Goal 02 independent: early cap differs")

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-tail-independent-v1",
        "candidate": primary["candidate"],
        "candidate_id": primary["candidate_id"],
        "evidence_label": "INDEPENDENT_EXACT_RECONSTRUCTION",
        "source_sha256": digest(Path(__file__).resolve()),
        "primary_source_sha256": digest(PRIMARY_SOURCE),
        "primary_receipt_sha256": digest(PRIMARY),
        "goal04_audit_sha256": digest(GOAL04),
        "recurrence_reconstruction": {
            "method": "binary Gaussian elimination and bit-by-bit accumulation",
            "parity_rank": rank_by_rows(parity),
            "stored_state_update_rank": rank_by_rows(state_update),
            "autonomous_output_rank": rank_by_rows(output_update),
            "basis_drive_pairs_checked": checked_pairs,
            "common_drive_cancellation": True,
        },
        "cutoffs": {
            "first_johnson_node": first_johnson,
            "maximum_late_list": maximum_late,
            "first_unique_node": first_unique,
            "cutoff_parameters": {
                "tail_length_bits": rows[JOHNSON_CUTOFF][2][0],
                "tail_distance_lower_bound": rows[JOHNSON_CUTOFF][2][1],
                "bad_tail_radius": rows[JOHNSON_CUTOFF][2][2],
            },
        },
        "ledger_replay": {
            "late_placement_probability": {
                "numerator": str(late_placement.numerator),
                "denominator": str(late_placement.denominator),
            },
            "late_aggregate_charge": {
                "numerator": str(late_charge.numerator),
                "denominator": str(late_charge.denominator),
            },
            "new_remaining_budget": {
                "numerator": str(new_budget.numerator),
                "denominator": str(new_budget.denominator),
            },
            "early_uniform_list_cap": early_cap,
            "sufficient_difference_spectrum_cap": early_cap - 1,
        },
        "primary_and_independent_agree": True,
        "status": "GOAL_02_PART_I_INDEPENDENTLY_VERIFIED",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    print(f"basis_drive_pairs_checked={checked_pairs}")
    print(f"first_johnson_node={first_johnson}")
    print(f"maximum_late_list={maximum_late}")
    print(f"first_unique_node={first_unique}")
    print(f"early_list_cap={early_cap}")
    print("status=GOAL_02_PART_I_INDEPENDENTLY_VERIFIED")


if __name__ == "__main__":
    main()
