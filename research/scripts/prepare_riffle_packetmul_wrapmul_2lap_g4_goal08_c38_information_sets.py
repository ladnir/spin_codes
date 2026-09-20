#!/usr/bin/env python3
"""Find the exact information-set packing number of a zero-anchor C38 shortening."""

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
NODES = 38
DIMENSION = 64
LENGTH = (NODES - 1) * 64
TARGET = 228


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
            raise RuntimeError("C38 packing: maintained set is dependent")

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
    return True


def pack(set_count: int, columns: list[int]) -> dict:
    sets = [[] for _ in range(set_count)]
    rejected = 0
    target_size = set_count * DIMENSION
    for element in range(LENGTH):
        if augment(element, sets, columns):
            if sum(map(len, sets)) == target_size:
                break
        else:
            rejected += 1
    for row in sets:
        eliminator(row, columns)
    union_size = sum(map(len, sets))
    return {
        "requested_set_count": set_count,
        "maximum_union_size": union_size,
        "full_information_set_count": sum(len(row) == DIMENSION for row in sets),
        "set_sizes": [len(row) for row in sets],
        "rejected_before_stop": rejected,
        "success": union_size == target_size,
        "sets": sets,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor-node", type=int, required=True)
    args = parser.parse_args()
    if not 0 <= args.anchor_node < NODES:
        raise SystemExit("C38 packing: anchor node must be in [0,37]")

    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(NODES, parity)
    zero_anchor_basis = kernel_basis(
        tuple((row >> (64 * args.anchor_node)) & MASK64 for row in full_rows),
        64,
    )
    if len(zero_anchor_basis) != DIMENSION:
        raise RuntimeError("C38 packing: zero-anchor dimension changed")
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
    if any(column == 0 for column in coordinate_columns):
        zero_coordinate_count = sum(column == 0 for column in coordinate_columns)
    else:
        zero_coordinate_count = 0

    theoretical_maximum = LENGTH // DIMENSION
    trials: list[dict] = []
    low = 0
    high = theoretical_maximum + 1
    successful: dict | None = None
    first_infeasible: dict | None = None
    while low + 1 < high:
        middle = (low + high) // 2
        trial = pack(middle, coordinate_columns)
        trials.append({key: value for key, value in trial.items() if key != "sets"})
        if trial["success"]:
            low = middle
            successful = trial
        else:
            high = middle
            first_infeasible = trial

    if low == 0:
        raise RuntimeError("C38 packing: no information set found")
    if successful is None or successful["requested_set_count"] != low:
        successful = pack(low, coordinate_columns)
        trials.append({key: value for key, value in successful.items() if key != "sets"})
    if high <= theoretical_maximum and (
        first_infeasible is None or first_infeasible["requested_set_count"] != high
    ):
        first_infeasible = pack(high, coordinate_columns)
        trials.append({key: value for key, value in first_infeasible.items() if key != "sets"})
    if not successful["success"]:
        raise RuntimeError("C38 packing: final feasible packing failed")
    if first_infeasible is not None and first_infeasible["success"]:
        raise RuntimeError("C38 packing: claimed infeasible packing succeeded")

    sets = successful["sets"]
    selected = [coordinate for row in sets for coordinate in row]
    if len(selected) != len(set(selected)):
        raise RuntimeError("C38 packing: information sets overlap")

    output = (
        CANDIDATE
        / "receipts"
        / f"goal08_c38_anchor_{args.anchor_node:02d}_information_sets.json"
    )
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal08-c38-information-sets-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_LINEAR_MATROID_BASE_PACKING",
        "source_sha256": digest(Path(__file__).resolve()),
        "anchor_node": args.anchor_node,
        "window_nodes": NODES,
        "distance_target": TARGET,
        "counterexample_weight_maximum": TARGET - 1,
        "shortened_code": {"length": LENGTH, "dimension": DIMENSION},
        "zero_coordinate_count": zero_coordinate_count,
        "theoretical_maximum_information_sets": theoretical_maximum,
        "maximum_disjoint_information_sets": low,
        "information_sets": sets,
        "search_trials": trials,
        "infeasible_successor": (
            None
            if first_infeasible is None
            else {key: value for key, value in first_infeasible.items() if key != "sets"}
        ),
        "result": "EXACT_MAXIMUM_DISJOINT_INFORMATION_SET_COUNT",
        "scope_limitation": (
            "The matroid-union computation determines the exact number of disjoint "
            "information sets in this zero-anchor shortening. It does not establish "
            "the C38 distance target."
        ),
    }
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"anchor_node={args.anchor_node}")
    print(f"maximum_disjoint_information_sets={low}")
    if first_infeasible is not None:
        print(
            "successor_union_size="
            f"{first_infeasible['maximum_union_size']}/"
            f"{first_infeasible['requested_set_count'] * DIMENSION}"
        )
    print(f"output={output}")
    print("status=EXACT_LINEAR_MATROID_BASE_PACKING")


if __name__ == "__main__":
    main()
