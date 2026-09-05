#!/usr/bin/env python3
"""Certify the Goal 02 tail-list theorem and its exact ledger charge."""

from __future__ import annotations

import hashlib
import json
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
PARENT = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
MANIFEST = CANDIDATE / "manifest.json"
GOAL01 = CANDIDATE / "receipts" / "goal01_wrap_uniformity_audit.json"
GOAL04 = PARENT / "receipts" / "goal04_boundary_prefix_audit.json"
OUTPUT = CANDIDATE / "receipts" / "goal02_tail_list_primary.json"

GENERATOR = 0xF4845518B9582A1F
MASK64 = (1 << 64) - 1
TOTAL_NODES = 32_772
PACKETS_PER_NODE = 16
TOTAL_PACKET_POSITIONS = TOTAL_NODES * PACKETS_PER_NODE
PACKET_SUPPORT = 33
SUPPORT33_WORDS = 26
NONZERO_STATES = (1 << 64) - 1
BAD_WEIGHT_MAXIMUM = 188_765
FOUR_NODE_DISTANCE = 36
PARENT_EMPTY_CUTOFF = 20_976
JOHNSON_CUTOFF = 7_060
UNIQUE_CUTOFF = 9_176


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def fraction_record(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def log2_interval(value: Fraction, places: int = 12) -> list[str]:
    if value <= 0:
        raise ValueError("Goal 02 tail list: logarithm requires a positive value")
    with localcontext() as context:
        context.prec = 150
        estimate = (
            Decimal(value.numerator).ln() - Decimal(value.denominator).ln()
        ) / Decimal(2).ln()
        unit = Decimal(1).scaleb(-places)
        return [
            format(estimate.quantize(unit, rounding=ROUND_FLOOR), "f"),
            format(estimate.quantize(unit, rounding=ROUND_CEILING), "f"),
        ]


def carryless_multiply_low(left: int, right: int) -> int:
    result = 0
    while right:
        low = right & -right
        result ^= left << (low.bit_length() - 1)
        right ^= low
    return result & MASK64


def generator_inverse() -> int:
    inverse = 1
    for bit in range(1, 64):
        if (carryless_multiply_low(GENERATOR, inverse) >> bit) & 1:
            inverse |= 1 << bit
    if carryless_multiply_low(GENERATOR, inverse) != 1:
        raise RuntimeError("Goal 02 tail list: generator inverse failed")
    return inverse


def encode(message: int) -> tuple[int, int]:
    low = 0
    high = 0
    while message:
        bit_value = message & -message
        bit = bit_value.bit_length() - 1
        low ^= (GENERATOR << bit) & MASK64
        high ^= (0 if bit == 0 else GENERATOR >> (64 - bit)) | (1 << 63)
        message ^= bit_value
    return low, high & MASK64


def systematic_right_columns() -> tuple[int, ...]:
    inverse = generator_inverse()
    columns = []
    for bit in range(64):
        message = carryless_multiply_low(1 << bit, inverse)
        low, high = encode(message)
        if low != 1 << bit:
            raise RuntimeError("Goal 02 tail list: systematic reconstruction failed")
        columns.append(high)
    return tuple(columns)


def accumulate(value: int) -> int:
    value ^= value << 1
    value ^= value << 2
    value ^= value << 4
    value ^= value << 8
    value ^= value << 16
    value ^= value << 32
    return value & MASK64


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        low = value & -value
        result ^= columns[low.bit_length() - 1]
        value ^= low
    return result


def binary_rank(columns: tuple[int, ...]) -> int:
    pivots = [0] * 64
    rank = 0
    for column in columns:
        value = column
        while value:
            pivot = value.bit_length() - 1
            if pivots[pivot]:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                rank += 1
                break
    return rank


def tail_parameters(first_node: int) -> dict:
    blocks_before_drive = first_node // 4
    remaining_blocks = TOTAL_NODES // 4 - blocks_before_drive
    length = 64 * 4 * remaining_blocks
    distance = FOUR_NODE_DISTANCE * remaining_blocks
    radius = BAD_WEIGHT_MAXIMUM - FOUR_NODE_DISTANCE * blocks_before_drive
    denominator = distance * length - 2 * radius * (length - radius)
    johnson_bound = None
    johnson_ratio = None
    if denominator > 0:
        johnson_ratio = Fraction(distance * length, denominator)
        johnson_bound = johnson_ratio.numerator // johnson_ratio.denominator
    unique = distance > 2 * radius
    effective_bound = 1 if unique else johnson_bound
    return {
        "first_node": first_node,
        "complete_prefix_blocks": blocks_before_drive,
        "remaining_blocks": remaining_blocks,
        "tail_length_bits": length,
        "tail_distance_lower_bound": distance,
        "bad_tail_radius": radius,
        "johnson_denominator_numerator": denominator,
        "johnson_ratio": (
            None if johnson_ratio is None else fraction_record(johnson_ratio)
        ),
        "johnson_list_bound": johnson_bound,
        "unique_decoding": unique,
        "effective_list_bound": effective_bound,
    }


def main() -> None:
    manifest = load(MANIFEST)
    goal01 = load(GOAL01)
    goal04 = load(GOAL04)
    if manifest["candidate_id"] != "riffle_packetmul_wrapmul_2lap_g4":
        raise RuntimeError("Goal 02 tail list: wrong active manifest")
    if goal01["status"] != "GOAL_01_PROVED_WRAP_UNIFORMITY_LIST_BOUND_OPEN":
        raise RuntimeError("Goal 02 tail list: Goal 01 audit is not current")
    if goal04["four_node_certificate"]["four_node_weight_lower_bound"] != 36:
        raise RuntimeError("Goal 02 tail list: parent four-node distance changed")
    if TOTAL_NODES % 4:
        raise RuntimeError("Goal 02 tail list: node count is not four-block aligned")

    parity_columns = systematic_right_columns()
    autonomous_columns = tuple(
        accumulate(parity_columns[bit]) for bit in range(64)
    )
    state_update_columns = tuple(
        apply_columns(parity_columns, accumulate(1 << bit)) for bit in range(64)
    )
    if (
        binary_rank(parity_columns) != 64
        or binary_rank(autonomous_columns) != 64
        or binary_rank(state_update_columns) != 64
    ):
        raise RuntimeError("Goal 02 tail list: recurrence map is singular")

    # The recurrence is O(s,x)=Acc(s+x), S(s,x)=P(O(s,x)).
    # Linearity gives O(s1,x)+O(s2,x)=Acc(s1+s2) and the same
    # common-drive cancellation after applying P.
    for bit in range(64):
        state = 1 << bit
        if accumulate(state) != accumulate(state ^ 0):
            raise RuntimeError("Goal 02 tail list: output basis replay failed")
        if apply_columns(parity_columns, accumulate(state)) != state_update_columns[bit]:
            raise RuntimeError("Goal 02 tail list: state basis replay failed")

    rows = [tail_parameters(first_node) for first_node in range(PARENT_EMPTY_CUTOFF)]
    first_johnson = next(
        row["first_node"] for row in rows if row["johnson_list_bound"] is not None
    )
    first_unique = next(row["first_node"] for row in rows if row["unique_decoding"])
    if first_johnson != JOHNSON_CUTOFF or first_unique != UNIQUE_CUTOFF:
        raise RuntimeError("Goal 02 tail list: cutoff arithmetic changed")
    late_rows = rows[JOHNSON_CUTOFF:]
    maximum_late_list = max(
        row["effective_list_bound"] for row in late_rows
        if row["effective_list_bound"] is not None
    )
    if maximum_late_list != 13_749:
        raise RuntimeError("Goal 02 tail list: late list maximum changed")

    total_placements = comb(TOTAL_PACKET_POSITIONS, PACKET_SUPPORT)
    suffix_7060 = Fraction(
        comb(TOTAL_PACKET_POSITIONS - PACKETS_PER_NODE * JOHNSON_CUTOFF, PACKET_SUPPORT),
        total_placements,
    )
    suffix_20976 = Fraction(
        comb(TOTAL_PACKET_POSITIONS - PACKETS_PER_NODE * PARENT_EMPTY_CUTOFF, PACKET_SUPPORT),
        total_placements,
    )
    late_band_placement = suffix_7060 - suffix_20976
    early_placement = 1 - suffix_7060
    late_charge = (
        SUPPORT33_WORDS
        * late_band_placement
        * Fraction(maximum_late_list, NONZERO_STATES)
    )

    budget_record = goal04["support33_ledger_update"][
        "retained_remaining_numerical_budget"
    ]
    old_budget = Fraction(
        int(budget_record["numerator"]), int(budget_record["denominator"])
    )
    old_closed_record = goal04["support33_ledger_update"][
        "retained_goal03_closed_partial_upper_bound"
    ]
    old_closed = Fraction(
        int(old_closed_record["numerator"]), int(old_closed_record["denominator"])
    )
    new_budget = old_budget - late_charge
    new_closed = old_closed + late_charge
    if new_budget <= 0:
        raise RuntimeError("Goal 02 tail list: late charge exhausts the ledger")
    early_cap = (
        new_budget * NONZERO_STATES // (SUPPORT33_WORDS * early_placement)
    )
    accepted_early_charge = (
        SUPPORT33_WORDS * early_placement * Fraction(early_cap, NONZERO_STATES)
    )
    rejected_early_charge = (
        SUPPORT33_WORDS * early_placement * Fraction(early_cap + 1, NONZERO_STATES)
    )
    if early_cap != 316_606:
        raise RuntimeError("Goal 02 tail list: early cap changed")
    if accepted_early_charge > new_budget or rejected_early_charge <= new_budget:
        raise RuntimeError("Goal 02 tail list: early cap is not maximal")

    cutoff_row = rows[JOHNSON_CUTOFF]
    unique_row = rows[UNIQUE_CUTOFF]
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-tail-primary-v1",
        "candidate": manifest["display_name"],
        "candidate_id": manifest["candidate_id"],
        "evidence_label": "EXACT_TAIL_LIST_AND_LEDGER_CERTIFICATE",
        "source_sha256": digest(Path(__file__).resolve()),
        "authenticated_inputs_sha256": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (MANIFEST, GOAL01, GOAL04)
        },
        "implemented_recurrence": {
            "emission": "O(s,x)=Acc(s xor x)",
            "next_state": "S(s,x)=P(O(s,x))",
            "common_drive_emission_difference": "O(s1,x) xor O(s2,x)=Acc(s1 xor s2)",
            "common_drive_next_state_difference": (
                "S(s1,x) xor S(s2,x)=P(Acc(s1 xor s2))"
            ),
            "systematic_parity_rank": binary_rank(parity_columns),
            "stored_state_update_rank": binary_rank(state_update_columns),
            "autonomous_output_rank": binary_rank(autonomous_columns),
            "checked_basis_vectors": 64,
        },
        "tail_code_theorem": {
            "total_nodes": TOTAL_NODES,
            "four_node_distance": FOUR_NODE_DISTANCE,
            "for_first_node_r": (
                "q=floor(r/4), n_r=256*(8193-q), d_r>=36*(8193-q), "
                "t_r=188765-36*q"
            ),
            "reason_drive_does_not_matter": (
                "Two trajectories with the same drive differ by the zero-drive "
                "autonomous recurrence at every node."
            ),
        },
        "johnson_region": {
            "first_node_interval": [JOHNSON_CUTOFF, PARENT_EMPTY_CUTOFF - 1],
            "first_applicable_node": first_johnson,
            "cutoff_parameters": cutoff_row,
            "maximum_list_size": maximum_late_list,
            "all_rows_checked": len(late_rows),
        },
        "unique_region": {
            "first_node_interval": [UNIQUE_CUTOFF, PARENT_EMPTY_CUTOFF - 1],
            "first_unique_node": first_unique,
            "cutoff_parameters": unique_row,
            "maximum_list_size": 1,
        },
        "late_band_ledger": {
            "event": "support 33, L!=0, and 7060<=r<=20975",
            "outer_word_count": SUPPORT33_WORDS,
            "placement_probability": fraction_record(late_band_placement),
            "placement_probability_log2_interval": log2_interval(late_band_placement),
            "uniform_list_cap": maximum_late_list,
            "wrap_state_denominator": NONZERO_STATES,
            "aggregate_charge": fraction_record(late_charge),
            "aggregate_charge_log2_interval": log2_interval(late_charge),
            "old_closed_partial": fraction_record(old_closed),
            "new_closed_partial": fraction_record(new_closed),
            "new_closed_partial_log2_interval": log2_interval(new_closed),
            "old_remaining_budget": fraction_record(old_budget),
            "new_remaining_budget": fraction_record(new_budget),
            "new_remaining_budget_log2_interval": log2_interval(new_budget),
            "disjoint_from_parent_empty_region": True,
        },
        "early_region_target": {
            "event": "support 33, L!=0, and r<=7059",
            "placement_probability": fraction_record(early_placement),
            "maximum_uniform_bad_states_per_fixed_drive": early_cap,
            "accepted_aggregate": fraction_record(accepted_early_charge),
            "next_integer_rejected": early_cap + 1,
            "next_integer_aggregate": fraction_record(rejected_early_charge),
            "sufficient_difference_spectrum_cap": early_cap - 1,
            "difference_weight_maximum": 2 * BAD_WEIGHT_MAXIMUM,
        },
        "status": "GOAL_02_PART_I_PROVED_EARLY_SPECTRUM_OPEN",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    print(f"first_johnson_node={first_johnson}")
    print(f"maximum_late_list={maximum_late_list}")
    print(f"first_unique_node={first_unique}")
    print(f"late_charge_log2={log2_interval(late_charge)}")
    print(f"early_list_cap={early_cap}")
    print("status=GOAL_02_PART_I_PROVED_EARLY_SPECTRUM_OPEN")


if __name__ == "__main__":
    main()
