#!/usr/bin/env python3
"""Native-XOR SAT probe for one systematic C24/C38 chart."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import deque
from pathlib import Path

from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

from audit_riffle_packetmul_wrapmul_2lap_bch_family import build_instance, observe
from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    observation_rows,
)
from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
RECEIPTS = CANDIDATE / "receipts"
WIDTH16_WITNESS = 0xCD7E70D6


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def xor_selected(rows: list[int] | tuple[int, ...], selector: int) -> int:
    result = 0
    while selector:
        bit = selector & -selector
        result ^= rows[bit.bit_length() - 1]
        selector ^= bit
    return result


def delete_packet(word: int, packet: int, width: int) -> int:
    low_width = width * packet
    low_mask = (1 << low_width) - 1
    return (word & low_mask) | ((word >> (low_width + width)) << low_width)


def kernel_basis(columns: tuple[int, ...], output_width: int) -> tuple[int, ...]:
    pivot_values = [0] * output_width
    pivot_representations = [0] * output_width
    result = []
    for column, original in enumerate(columns):
        value = original
        representation = 1 << column
        while value:
            pivot = value.bit_length() - 1
            if pivot_values[pivot]:
                value ^= pivot_values[pivot]
                representation ^= pivot_representations[pivot]
            else:
                pivot_values[pivot] = value
                pivot_representations[pivot] = representation
                break
        if value == 0:
            result.append(representation)
    return tuple(result)


def rank(columns: list[int]) -> int:
    pivots: dict[int, int] = {}
    for original in columns:
        value = original
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


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
            raise RuntimeError("Goal 09: maintained information set is dependent")

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


def pack_information_sets(
    count: int, dimension: int, columns: list[int]
) -> list[list[int]] | None:
    sets = [[] for _ in range(count)]
    target = count * dimension
    for coordinate in range(len(columns)):
        augment(coordinate, sets, columns)
        if sum(map(len, sets)) == target:
            break
    if sum(map(len, sets)) != target:
        return None
    for row in sets:
        if len(row) != dimension:
            raise RuntimeError("Goal 09: information-set size mismatch")
        eliminator(row, columns)
    return sets


def width16_instance() -> dict:
    width = 16
    nodes = 24
    instance = build_instance(5, 7)
    full_rows = [
        sum(value << (width * node) for node, value in enumerate(
            observe(1 << bit, width, instance["columns"], nodes)
        ))
        for bit in range(2 * width)
    ]
    witness_word = xor_selected(full_rows, WIDTH16_WITNESS)
    if witness_word.bit_count() != 121:
        raise RuntimeError("Goal 09: width-16 minimum witness replay changed")

    selected = None
    for anchor in range(nodes):
        anchor_value = (witness_word >> (width * anchor)) & ((1 << width) - 1)
        if anchor_value.bit_count() > 5:
            continue
        anchor_columns = tuple(
            (row >> (width * anchor)) & ((1 << width) - 1) for row in full_rows
        )
        basis = kernel_basis(anchor_columns, width)
        if len(basis) != width:
            raise RuntimeError("Goal 09: width-16 anchor kernel dimension changed")
        shortened_generators = tuple(
            delete_packet(xor_selected(full_rows, state), anchor, width)
            for state in basis
        )
        shortened_length = (nodes - 1) * width
        coordinate_columns = [
            sum(
                ((generator >> coordinate) & 1) << bit
                for bit, generator in enumerate(shortened_generators)
            )
            for coordinate in range(shortened_length)
        ]
        sets = None
        information_set_count = 0
        for count in range(nodes - 1, 0, -1):
            sets = pack_information_sets(count, width, coordinate_columns)
            if sets is not None:
                information_set_count = count
                break
        if sets is None:
            raise RuntimeError("Goal 09: width-16 shortening has no information set")
        shortened_witness = delete_packet(witness_word, anchor, width)
        for set_index, information_set in enumerate(sets):
            restriction_weight = sum(
                (shortened_witness >> coordinate) & 1
                for coordinate in information_set
            )
            if restriction_weight <= 5:
                selected = (
                    anchor,
                    set_index,
                    information_set,
                    restriction_weight,
                    information_set_count,
                )
                break
        if selected is not None:
            break
    if selected is None:
        raise RuntimeError("Goal 09: no width-16 validation chart covers the witness")
    anchor, set_index, information_set, restriction_weight, information_set_count = selected
    return {
        "label": "width16_c24",
        "width": width,
        "nodes": nodes,
        "dimension": 2 * width,
        "full_rows": full_rows,
        "anchor": anchor,
        "set_index": set_index,
        "information_set": information_set,
        "anchor_cap": 5,
        "restriction_cap": 5,
        "known_witness": WIDTH16_WITNESS,
        "known_witness_restriction_weight": restriction_weight,
        "chart_source": (
            f"deterministic {information_set_count}-base packing computed by this script"
        ),
    }


def width64_instance(anchor: int, set_index: int) -> dict:
    width = 64
    nodes = 38
    receipt_path = RECEIPTS / f"goal08_c38_anchor_{anchor:02d}_uniform35.json"
    receipt = json.loads(receipt_path.read_text())
    sets = receipt["information_sets"]
    if not 0 <= set_index < len(sets):
        raise SystemExit(f"Goal 09: set index must be in [0,{len(sets) - 1}]")
    parity = build_apply(systematic_state_columns())
    return {
        "label": "width64_c38",
        "width": width,
        "nodes": nodes,
        "dimension": 2 * width,
        "full_rows": observation_rows(nodes, parity),
        "anchor": anchor,
        "set_index": set_index,
        "information_set": sets[set_index],
        "anchor_cap": 5,
        "restriction_cap": 6,
        "known_witness": None,
        "known_witness_restriction_weight": None,
        "chart_source": str(receipt_path.relative_to(ROOT)),
        "chart_source_sha256": digest(receipt_path),
    }


def full_coordinates(instance: dict) -> list[int]:
    width = instance["width"]
    anchor = instance["anchor"]
    return [
        coordinate if coordinate < anchor * width else coordinate + width
        for coordinate in instance["information_set"]
    ]


def coordinate_columns(rows: list[int], coordinates: list[int]) -> list[int]:
    return [
        sum(((row >> coordinate) & 1) << bit for bit, row in enumerate(rows))
        for coordinate in coordinates
    ]


def chart_basis_states(rows: list[int], coordinates: list[int]) -> list[int]:
    dimension = len(rows)
    forms = coordinate_columns(rows, coordinates)
    map_columns = [
        sum(((form >> bit) & 1) << coordinate for coordinate, form in enumerate(forms))
        for bit in range(dimension)
    ]
    pivot_values: dict[int, int] = {}
    pivot_states: dict[int, int] = {}
    for bit, original in enumerate(map_columns):
        value = original
        state = 1 << bit
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivot_values:
                value ^= pivot_values[pivot]
                state ^= pivot_states[pivot]
            else:
                pivot_values[pivot] = value
                pivot_states[pivot] = state
                break
    if len(pivot_values) != dimension:
        raise RuntimeError("Goal 09: chart map is singular")

    result = []
    for coordinate in range(dimension):
        value = 1 << coordinate
        state = 0
        while value:
            pivot = value.bit_length() - 1
            value ^= pivot_values[pivot]
            state ^= pivot_states[pivot]
        result.append(state)
    for coordinate, state in enumerate(result):
        image = sum(
            ((xor_selected(rows, state) >> output_coordinate) & 1) << bit
            for bit, output_coordinate in enumerate(coordinates)
        )
        if image != 1 << coordinate:
            raise RuntimeError("Goal 09: chart inverse replay failed")
    return result


def add_at_most(
    clauses: list[list[int]], literals: list[int], bound: int, top_variable: int
) -> int:
    cardinality = CardEnc.atmost(
        lits=literals,
        bound=bound,
        top_id=top_variable,
        encoding=EncType.seqcounter,
    )
    clauses.extend(cardinality.clauses)
    return cardinality.nv


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instance", choices=("width16", "width64"), required=True)
    parser.add_argument("--weight-bound", type=int, required=True)
    parser.add_argument("--anchor-node", type=int, default=0)
    parser.add_argument("--set-index", type=int, default=0)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    args = parser.parse_args()
    if args.timeout_seconds <= 0 or args.weight_bound < 0:
        raise SystemExit("Goal 09: timeout must be positive and bound nonnegative")
    if args.instance == "width16":
        instance = width16_instance()
    else:
        if not 0 <= args.anchor_node < 38:
            raise SystemExit("Goal 09: width-64 anchor must be in [0,37]")
        instance = width64_instance(args.anchor_node, args.set_index)

    width = instance["width"]
    nodes = instance["nodes"]
    dimension = instance["dimension"]
    length = width * nodes
    rows = instance["full_rows"]
    anchor = instance["anchor"]
    restriction_coordinates = full_coordinates(instance)
    anchor_coordinates = list(range(width * anchor, width * (anchor + 1)))
    chart_coordinates = anchor_coordinates + restriction_coordinates
    chart_rank = rank(coordinate_columns(rows, chart_coordinates))
    if chart_rank != dimension:
        raise RuntimeError(
            f"Goal 09: selected chart rank {chart_rank} is not {dimension}"
        )

    basis_states = chart_basis_states(rows, chart_coordinates)
    chart_generators = [xor_selected(rows, state) for state in basis_states]
    information_variables = list(range(1, dimension + 1))
    chart_variable_by_coordinate = {
        coordinate: information_variables[bit]
        for bit, coordinate in enumerate(chart_coordinates)
    }
    next_variable = dimension + 1
    output_literals = []
    nonsystematic_outputs: list[tuple[int, int]] = []
    for coordinate in range(length):
        if coordinate in chart_variable_by_coordinate:
            output_literals.append(chart_variable_by_coordinate[coordinate])
        else:
            output_literals.append(next_variable)
            nonsystematic_outputs.append((coordinate, next_variable))
            next_variable += 1
    clauses: list[list[int]] = [information_variables]
    top_variable = next_variable - 1
    top_variable = add_at_most(
        clauses, output_literals, args.weight_bound, top_variable
    )
    top_variable = add_at_most(
        clauses,
        information_variables[:width],
        instance["anchor_cap"],
        top_variable,
    )
    top_variable = add_at_most(
        clauses,
        information_variables[width:],
        instance["restriction_cap"],
        top_variable,
    )

    encoding_seconds = 0.0
    solve_seconds = 0.0
    started = time.perf_counter()
    with Solver(name="cms", bootstrap_with=clauses) as solver:
        xor_widths = []
        for coordinate, output_variable in nonsystematic_outputs:
            inputs = [
                bit + 1
                for bit, generator in enumerate(chart_generators)
                if (generator >> coordinate) & 1
            ]
            xor_widths.append(len(inputs))
            solver.add_xor_clause(inputs + [output_variable], value=False)
        encoding_seconds = time.perf_counter() - started
        solver.solver.time_budget(args.timeout_seconds)
        solver.conf_budget((1 << 31) - 1)
        solve_started = time.perf_counter()
        satisfiable = solver.solve_limited()
        solve_seconds = time.perf_counter() - solve_started
        model = solver.get_model() if satisfiable is True else None
        try:
            statistics = solver.accum_stats()
        except NotImplementedError:
            statistics = {}

    result = "SAT" if satisfiable is True else "UNSAT" if satisfiable is False else "UNKNOWN"
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal09-chart-sat-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_SAT_RESULT_OR_TIMEOUT_DIAGNOSTIC",
        "source_sha256": digest(Path(__file__).resolve()),
        "instance": instance["label"],
        "code": {"length": length, "dimension": dimension, "packet_width": width},
        "chart": {
            "anchor_node": anchor,
            "information_set_index": instance["set_index"],
            "information_set_shortened_coordinates": instance["information_set"],
            "rank": chart_rank,
            "source": instance["chart_source"],
            "source_sha256": instance.get("chart_source_sha256"),
        },
        "decision_problem": {
            "state_nonzero": True,
            "total_output_weight_upper_bound": args.weight_bound,
            "anchor_packet_weight_upper_bound": instance["anchor_cap"],
            "information_set_restriction_weight_upper_bound": instance["restriction_cap"],
        },
        "encoding": {
            "solver": "CryptoMiniSat through python-sat",
            "linear_constraints": "native XOR",
            "cardinality_encoding": "sequential counter",
            "information_variables": dimension,
            "systematic_output_aliases": dimension,
            "nonsystematic_output_variables": len(nonsystematic_outputs),
            "cnf_variables": top_variable,
            "cnf_clauses_before_xors": len(clauses),
            "native_xor_constraints": len(nonsystematic_outputs),
            "xor_width_minimum": min(xor_widths),
            "xor_width_maximum": max(xor_widths),
        },
        "known_width16_witness": (
            None
            if instance["known_witness"] is None
            else {
                "state_hex": hex(instance["known_witness"]),
                "total_weight": xor_selected(rows, instance["known_witness"]).bit_count(),
                "anchor_weight": (
                    (xor_selected(rows, instance["known_witness"]) >> (width * anchor))
                    & ((1 << width) - 1)
                ).bit_count(),
                "restriction_weight": instance["known_witness_restriction_weight"],
            }
        ),
        "timeout_seconds": args.timeout_seconds,
        "encoding_seconds": encoding_seconds,
        "solve_seconds": solve_seconds,
        "elapsed_seconds": encoding_seconds + solve_seconds,
        "solver_statistics": statistics,
        "result": result,
        "scope_limitation": (
            "SAT is replayed exactly. UNSAT is an exact solver conclusion but this "
            "receipt does not contain an independently checkable proof. UNKNOWN "
            "records only a timeout. A single width-64 chart is only one of 1330 charts."
        ),
    }
    if model is not None:
        positive = {literal for literal in model if literal > 0}
        chart_selector = sum(
            1 << bit
            for bit, variable in enumerate(information_variables)
            if variable in positive
        )
        state = xor_selected(basis_states, chart_selector)
        word = xor_selected(rows, state)
        node_weights = [
            ((word >> (width * node)) & ((1 << width) - 1)).bit_count()
            for node in range(nodes)
        ]
        restriction_weight = sum(
            (word >> coordinate) & 1 for coordinate in restriction_coordinates
        )
        if (
            state == 0
            or sum(node_weights) > args.weight_bound
            or node_weights[anchor] > instance["anchor_cap"]
            or restriction_weight > instance["restriction_cap"]
        ):
            raise RuntimeError("Goal 09: SAT witness replay failed")
        payload["witness"] = {
            "state_hex": hex(state),
            "initial_a_hex": hex(state & ((1 << width) - 1)),
            "initial_b_hex": hex(state >> width),
            "node_weights": node_weights,
            "total_weight": sum(node_weights),
            "anchor_weight": node_weights[anchor],
            "restriction_weight": restriction_weight,
        }

    stem = (
        f"goal09_chart_sat_{instance['label']}_anchor_{anchor:02d}_"
        f"set_{instance['set_index']:02d}_bound_{args.weight_bound:03d}.json"
    )
    output = RECEIPTS / stem
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"instance={instance['label']}", flush=True)
    print(f"chart=anchor_{anchor:02d}_set_{instance['set_index']:02d}", flush=True)
    print(f"result={result}", flush=True)
    print(f"encoding_seconds={encoding_seconds:.6f}", flush=True)
    print(f"solve_seconds={solve_seconds:.6f}", flush=True)
    if model is not None:
        print(f"witness={payload['witness']}", flush=True)
    print(f"output={output}", flush=True)
    print("status=SYSTEMATIC_CHART_NATIVE_XOR_DECISION", flush=True)


if __name__ == "__main__":
    main()
