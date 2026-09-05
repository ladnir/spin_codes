#!/usr/bin/env python3
"""Pack disjoint information sets in the zero-anchor C24 shortening."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import deque
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
LENGTH = 23 * 64
SET_COUNT = 22


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def eliminator(elements: list[int], columns: list[int]):
    pivot_values: dict[int, int] = {}
    pivot_representations: dict[int, int] = {}
    for position, element in enumerate(elements):
        value = columns[element]
        representation = 1 << position
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivot_values:
                value ^= pivot_values[pivot]
                representation ^= pivot_representations[pivot]
            else:
                pivot_values[pivot] = value
                pivot_representations[pivot] = representation
                break
        if value == 0:
            raise RuntimeError("zero-anchor sets: maintained set is dependent")

    def reduce(value: int) -> tuple[bool, int]:
        representation = 0
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivot_values:
                return True, representation
            value ^= pivot_values[pivot]
            representation ^= pivot_representations[pivot]
        return False, representation

    return reduce


def augment(start: int, sets: list[list[int]], columns: list[int]) -> bool:
    reducers = [eliminator(row, columns) for row in sets]
    queue = deque([start])
    parent: dict[int, tuple[int, int] | None] = {start: None}
    terminal: tuple[int, int] | None = None
    while queue and terminal is None:
        element = queue.popleft()
        for set_index, (row, reduce) in enumerate(zip(sets, reducers, strict=True)):
            independent, circuit = reduce(columns[element])
            if independent:
                terminal = (element, set_index)
                break
            while circuit:
                bit = circuit & -circuit
                displaced = row[bit.bit_length() - 1]
                circuit ^= bit
                if displaced not in parent:
                    parent[displaced] = (element, set_index)
                    queue.append(displaced)
    if terminal is None:
        return False

    path_elements = [terminal[0]]
    path_sets = []
    current = terminal[0]
    while parent[current] is not None:
        previous, set_index = parent[current]
        path_elements.append(previous)
        path_sets.append(set_index)
        current = previous
    path_elements.reverse()
    path_sets.reverse()

    sets[terminal[1]].append(path_elements[-1])
    for edge in range(len(path_sets) - 1, -1, -1):
        set_index = path_sets[edge]
        incoming = path_elements[edge]
        displaced = path_elements[edge + 1]
        position = sets[set_index].index(displaced)
        sets[set_index][position] = incoming

    for row in sets:
        eliminator(row, columns)
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor-node", type=int, default=0)
    args = parser.parse_args()
    if not 0 <= args.anchor_node < NODES:
        raise SystemExit("zero-anchor sets: anchor node must be in [0,23]")
    output = (
        CANDIDATE
        / "receipts"
        / f"goal06_zero_anchor_node_{args.anchor_node:02d}_information_sets.json"
    )
    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(NODES, parity)
    zero_anchor_basis = kernel_basis(
        tuple((row >> (64 * args.anchor_node)) & MASK64 for row in full_rows),
        64,
    )
    if len(zero_anchor_basis) != DIMENSION:
        raise RuntimeError("zero-anchor sets: shortening dimension changed")
    generators = tuple(
        delete_node(xor_selected(full_rows, state), args.anchor_node)
        for state in zero_anchor_basis
    )
    coordinate_columns = [
        sum(
            ((generator >> coordinate) & 1) << bit
            for bit, generator in enumerate(generators)
        )
        for coordinate in range(LENGTH)
    ]

    sets = [[] for _ in range(SET_COUNT)]
    selected: set[int] = set()
    rejected = []
    for element in range(LENGTH):
        if augment(element, sets, coordinate_columns):
            selected.add(element)
        else:
            rejected.append(element)
    selected = {element for row in sets for element in row}
    if len(selected) != sum(len(row) for row in sets):
        raise RuntimeError("zero-anchor sets: packed sets overlap")
    maximum_union_size = len(selected)
    full_set_count = sum(len(row) == DIMENSION for row in sets)
    if any(len(row) > DIMENSION for row in sets):
        raise RuntimeError("zero-anchor sets: independent set exceeds rank")

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal06-zero-anchor-sets-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_LINEAR_MATROID_UNION",
        "source_sha256": digest(Path(__file__).resolve()),
        "anchor_node": args.anchor_node,
        "shortened_code": {"length": LENGTH, "dimension": DIMENSION},
        "requested_set_count": SET_COUNT,
        "set_sizes": [len(row) for row in sets],
        "full_information_set_count": full_set_count,
        "full_information_sets": [row for row in sets if len(row) == DIMENSION],
        "maximum_union_size_for_22_MATROID_COPIES": maximum_union_size,
        "rejected_coordinate_count": LENGTH - maximum_union_size,
        "information_sets": sets if full_set_count == SET_COUNT else None,
        "result": (
            "TWENTY_TWO_DISJOINT_INFORMATION_SETS"
            if full_set_count == SET_COUNT
            else "MAXIMUM_UNION_SMALLER_THAN_TWENTY_TWO_BASES"
        ),
        "scope_limitation": (
            "The augmenting-path computation certifies maximum cardinality in "
            "the union of 22 copies of the represented linear matroid."
        ),
    }
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"set_sizes={payload['set_sizes']}")
    print(f"maximum_union_size={maximum_union_size}")
    print(f"result={payload['result']}")
    print(f"output={output}")
    print("status=EXACT_LINEAR_MATROID_UNION")


if __name__ == "__main__":
    main()
