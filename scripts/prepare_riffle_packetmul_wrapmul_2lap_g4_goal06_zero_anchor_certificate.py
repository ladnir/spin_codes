#!/usr/bin/env python3
"""Prepare systematic rows for the exact zero-anchor C24 certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
from pathlib import Path

from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    observation_rows,
)
from analyze_systematic_group_kernel import systematic_state_columns
from probe_riffle_packetmul_wrapmul_2lap_g4_goal05_zero_anchor import (
    MASK64,
    delete_node,
    kernel_basis,
    xor_selected,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
NODES = 24
DIMENSION = 64
WORDS = 23
BOUND = 143


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank(values: list[int]) -> int:
    pivots: dict[int, int] = {}
    result = 0
    for original in values:
        value = original
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                result += 1
                break
    return result


def inverse_columns(columns: tuple[int, ...]) -> tuple[int, ...]:
    if len(columns) != DIMENSION or rank(list(columns)) != DIMENSION:
        raise RuntimeError("zero-anchor certificate input: singular information set")
    pivot_values: dict[int, int] = {}
    pivot_representations: dict[int, int] = {}
    for input_bit, original in enumerate(columns):
        value = original
        representation = 1 << input_bit
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivot_values:
                value ^= pivot_values[pivot]
                representation ^= pivot_representations[pivot]
            else:
                pivot_values[pivot] = value
                pivot_representations[pivot] = representation
                break
    result = []
    for output_bit in range(DIMENSION):
        value = 1 << output_bit
        representation = 0
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivot_values:
                raise RuntimeError("zero-anchor certificate input: inverse failed")
            value ^= pivot_values[pivot]
            representation ^= pivot_representations[pivot]
        result.append(representation)
    return tuple(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor-node", type=int, default=0)
    args = parser.parse_args()
    if not 0 <= args.anchor_node < NODES:
        raise SystemExit("zero-anchor certificate input: anchor must be in [0,23]")
    sets_receipt_path = (
        CANDIDATE
        / "receipts"
        / f"goal06_zero_anchor_node_{args.anchor_node:02d}_information_sets.json"
    )
    binary_output = (
        CANDIDATE
        / "receipts"
        / f"goal06_zero_anchor_node_{args.anchor_node:02d}_c24.bin"
    )
    json_output = binary_output.with_suffix(".json")
    sets_receipt = json.loads(sets_receipt_path.read_text())
    information_sets = sets_receipt["full_information_sets"]
    set_count = len(information_sets)
    if set_count < 1 or any(
        len(row) != DIMENSION for row in information_sets
    ):
        raise RuntimeError("zero-anchor certificate input: set dimensions changed")
    base_weight, remainder = divmod(BOUND, set_count)
    base_weight_set_count = remainder + 1
    if base_weight > 7:
        raise RuntimeError("zero-anchor certificate input: enumeration exceeds weight 7")
    expected_candidates = (
        set_count
        * sum(math.comb(DIMENSION, weight) for weight in range(base_weight))
        + base_weight_set_count * math.comb(DIMENSION, base_weight)
    )

    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(NODES, parity)
    state_basis = kernel_basis(
        tuple((row >> (64 * args.anchor_node)) & MASK64 for row in full_rows),
        64,
    )
    if len(state_basis) != DIMENSION:
        raise RuntimeError("zero-anchor certificate input: shortening changed")
    generators = tuple(
        delete_node(xor_selected(full_rows, state), args.anchor_node)
        for state in state_basis
    )

    systematic_rows = []
    systematic_states = []
    for information_set in information_sets:
        selected_map_columns = tuple(
            sum(
                ((generators[input_bit] >> coordinate) & 1) << position
                for position, coordinate in enumerate(information_set)
            )
            for input_bit in range(DIMENSION)
        )
        inverse = inverse_columns(selected_map_columns)
        rows = tuple(xor_selected(generators, message) for message in inverse)
        states = tuple(xor_selected(state_basis, message) for message in inverse)
        for row_index, row in enumerate(rows):
            selected = sum(
                ((row >> coordinate) & 1) << position
                for position, coordinate in enumerate(information_set)
            )
            if selected != 1 << row_index:
                raise RuntimeError("zero-anchor certificate input: identity failed")
        systematic_rows.append(rows)
        systematic_states.append(states)

    with binary_output.open("wb") as handle:
        handle.write(b"RZ24Z02\0")
        handle.write(
            struct.pack(
                "<IIIIIIQ",
                set_count,
                DIMENSION,
                WORDS,
                BOUND,
                base_weight,
                base_weight_set_count,
                expected_candidates,
            )
        )
        for rows, states in zip(systematic_rows, systematic_states, strict=True):
            for row, state in zip(rows, states, strict=True):
                handle.write(
                    struct.pack(
                        "<23Q",
                        *(row >> (64 * word) & MASK64 for word in range(WORDS)),
                    )
                )
                handle.write(struct.pack("<2Q", state & MASK64, state >> 64))

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal06-zero-anchor-input-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_SYSTEMATIC_CERTIFICATE_INPUT",
        "source_sha256": digest(Path(__file__).resolve()),
        "sets_receipt_sha256": digest(sets_receipt_path),
        "anchor_node": args.anchor_node,
        "shortened_code": {"length": WORDS * 64, "dimension": DIMENSION},
        "information_set_count": set_count,
        "base_restriction_weight": base_weight,
        "fixed_base_weight_set_count": base_weight_set_count,
        "bound": BOUND,
        "expected_candidate_count": expected_candidates,
        "binary_file": binary_output.name,
        "binary_sha256": digest(binary_output),
        "binary_bytes": binary_output.stat().st_size,
    }
    json_output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"binary={binary_output}")
    print(f"output={json_output}")
    print("status=EXACT_SYSTEMATIC_CERTIFICATE_INPUT")


if __name__ == "__main__":
    main()
