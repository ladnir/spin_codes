#!/usr/bin/env python3
"""Prepare the zero-anchor odd-coset distance certificate input."""

from __future__ import annotations

import hashlib
import json
import random
import struct
from pathlib import Path

from analyze_riffle_packetmul_wrapmul_2lap_g4_goal05_zero_anchor_even_trace import (
    EVEN_DIMENSION,
    EVEN_NODES,
    ODD_NODES,
    image_basis,
    select_nodes,
)
from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    observation_rows,
)
from analyze_systematic_group_kernel import systematic_state_columns
from probe_riffle_packetmul_wrapmul_2lap_g4_goal05_zero_anchor import (
    MASK64,
    kernel_basis,
    xor_selected,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
BINARY_OUTPUT = CANDIDATE / "receipts" / "goal05_zero_anchor_odd_cosets.bin"
JSON_OUTPUT = BINARY_OUTPUT.with_suffix(".json")
SETS = 14
DIMENSION = 40
WORDS = 9
LENGTH = 576
THRESHOLD = 105
SEED = 0x4F44445F434F5331


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
    if rank(list(columns)) != DIMENSION:
        raise RuntimeError("odd cosets: information matrix is singular")
    pivot_values = [0] * DIMENSION
    pivot_representations = [0] * DIMENSION
    for input_bit, original in enumerate(columns):
        value = original
        representation = 1 << input_bit
        while value:
            pivot = value.bit_length() - 1
            if pivot_values[pivot]:
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
            if not pivot_values[pivot]:
                raise RuntimeError("odd cosets: inverse reduction failed")
            value ^= pivot_values[pivot]
            representation ^= pivot_representations[pivot]
        result.append(representation)
    return tuple(result)


def disjoint_information_sets(columns: list[int]) -> tuple[list[list[int]], int]:
    for attempt in range(10_000):
        rng = random.Random(SEED + attempt)
        pool = list(range(LENGTH))
        result = []
        success = True
        for _ in range(SETS):
            rng.shuffle(pool)
            selected = []
            selected_values = []
            selected_set = set()
            current_rank = 0
            for coordinate in pool:
                candidate_values = selected_values + [columns[coordinate]]
                candidate_rank = rank(candidate_values)
                if candidate_rank == current_rank + 1:
                    selected.append(coordinate)
                    selected_values.append(columns[coordinate])
                    selected_set.add(coordinate)
                    current_rank = candidate_rank
                    if current_rank == DIMENSION:
                        break
            if current_rank != DIMENSION:
                success = False
                break
            result.append(selected)
            pool = [coordinate for coordinate in pool if coordinate not in selected_set]
        if success:
            return result, attempt
    raise RuntimeError("odd cosets: no 14-set packing in 10000 attempts")


def main() -> None:
    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(18, parity)
    zero_anchor_basis = kernel_basis(tuple(row & MASK64 for row in full_rows), 64)
    full_words = tuple(xor_selected(full_rows, state) for state in zero_anchor_basis)
    even_columns = tuple(select_nodes(word, EVEN_NODES) for word in full_words)
    even_basis, quotient_state_basis = image_basis(even_columns, zero_anchor_basis)
    if len(even_basis) != EVEN_DIMENSION:
        raise RuntimeError("odd cosets: even rank changed")
    even_kernel_coordinates = kernel_basis(even_columns, 64 * len(EVEN_NODES))
    even_kernel_states = tuple(
        xor_selected(zero_anchor_basis, coordinates)
        for coordinates in even_kernel_coordinates
    )
    if len(even_kernel_states) != DIMENSION:
        raise RuntimeError("odd cosets: kernel dimension changed")
    odd_generators = tuple(
        select_nodes(xor_selected(full_rows, state), ODD_NODES)
        for state in even_kernel_states
    )
    if rank(list(odd_generators)) != DIMENSION:
        raise RuntimeError("odd cosets: odd kernel observation is not injective")

    coordinate_columns = [
        sum(
            ((generator >> coordinate) & 1) << bit
            for bit, generator in enumerate(odd_generators)
        )
        for coordinate in range(LENGTH)
    ]
    information_sets, successful_attempt = disjoint_information_sets(
        coordinate_columns
    )
    if len({coordinate for row in information_sets for coordinate in row}) != SETS * DIMENSION:
        raise RuntimeError("odd cosets: information sets overlap")

    systematic_generators = []
    systematic_states = []
    for information_set in information_sets:
        selected_map_columns = tuple(
            sum(
                ((odd_generators[input_bit] >> coordinate) & 1) << position
                for position, coordinate in enumerate(information_set)
            )
            for input_bit in range(DIMENSION)
        )
        inverse = inverse_columns(selected_map_columns)
        rows = tuple(xor_selected(odd_generators, message) for message in inverse)
        states = tuple(
            xor_selected(even_kernel_states, message) for message in inverse
        )
        for row_index, row in enumerate(rows):
            selected = sum(
                ((row >> coordinate) & 1) << position
                for position, coordinate in enumerate(information_set)
            )
            if selected != 1 << row_index:
                raise RuntimeError("odd cosets: systematic identity failed")
        systematic_generators.append(rows)
        systematic_states.append(states)

    low_cosets = [(0, 0, 0)]
    previous_gray = 0
    even_word = 0
    state = 0
    for index in range(1, 1 << EVEN_DIMENSION):
        gray = index ^ (index >> 1)
        changed = gray ^ previous_gray
        bit = changed.bit_length() - 1
        even_word ^= even_basis[bit]
        state ^= quotient_state_basis[bit]
        previous_gray = gray
        even_weight = even_word.bit_count()
        if even_weight <= THRESHOLD:
            odd_word = select_nodes(xor_selected(full_rows, state), ODD_NODES)
            low_cosets.append((even_weight, odd_word, state))
    low_cosets.sort(key=lambda row: (row[0], row[1], row[2]))
    if len(low_cosets) != 1025:
        raise RuntimeError(f"odd cosets: low coset count {len(low_cosets)} changed")

    with BINARY_OUTPUT.open("wb") as handle:
        handle.write(b"RZOC01\0\0")
        handle.write(struct.pack("<IIIII", SETS, DIMENSION, WORDS, len(low_cosets), THRESHOLD))
        for information_set in information_sets:
            handle.write(struct.pack("<40H", *information_set))
        for rows, states in zip(
            systematic_generators, systematic_states, strict=True
        ):
            for row, state in zip(rows, states, strict=True):
                handle.write(
                    struct.pack(
                        "<9Q", *(row >> (64 * word) & MASK64 for word in range(WORDS))
                    )
                )
                handle.write(struct.pack("<2Q", state & MASK64, state >> 64))
        for even_weight, odd_word, state in low_cosets:
            handle.write(struct.pack("<H6x", even_weight))
            handle.write(
                struct.pack(
                    "<9Q",
                    *(odd_word >> (64 * word) & MASK64 for word in range(WORDS)),
                )
            )
            handle.write(struct.pack("<2Q", state & MASK64, state >> 64))

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal05-odd-cosets-input-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_AFFINE_COSET_CERTIFICATE_INPUT",
        "source_sha256": digest(Path(__file__).resolve()),
        "zero_anchor_node": 0,
        "odd_trace_code": {"length": LENGTH, "dimension": DIMENSION},
        "low_even_coset_count": len(low_cosets),
        "information_set_count": SETS,
        "information_set_size": DIMENSION,
        "information_sets_are_disjoint": True,
        "unused_coordinate_count": LENGTH - SETS * DIMENSION,
        "packing_seed_hex": hex(SEED),
        "packing_successful_attempt": successful_attempt,
        "binary_file": BINARY_OUTPUT.name,
        "binary_sha256": digest(BINARY_OUTPUT),
        "binary_bytes": BINARY_OUTPUT.stat().st_size,
        "certificate_search_count": 482_617_685,
    }
    JSON_OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"packing_successful_attempt={successful_attempt}")
    print(f"low_even_coset_count={len(low_cosets)}")
    print(f"binary={BINARY_OUTPUT}")
    print(f"output={JSON_OUTPUT}")
    print("status=EXACT_AFFINE_COSET_CERTIFICATE_INPUT")


if __name__ == "__main__":
    main()
