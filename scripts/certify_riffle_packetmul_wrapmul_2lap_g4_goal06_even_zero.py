#!/usr/bin/env python3
"""Exactly exclude C24 weight at most 143 in the even-zero subcode."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import random
import time
from pathlib import Path

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
from probe_riffle_packetmul_wrapmul_2lap_g4_goal06_even_zero import select_nodes


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
OUTPUT = CANDIDATE / "receipts" / "goal06_even_zero_c24_exact.json"
NODES = 24
BOUND = 143
DIMENSION = 31
LENGTH = 768
INFORMATION_SETS = 24
MAXIMUM_RESTRICTION_WEIGHT = BOUND // INFORMATION_SETS
SEED = 0x4556454E5A45524F


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
        raise RuntimeError("even-zero exact: singular information set")
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
                raise RuntimeError("even-zero exact: inverse reduction failed")
            value ^= pivot_values[pivot]
            representation ^= pivot_representations[pivot]
        result.append(representation)
    return tuple(result)


def pack_information_sets(columns: list[int]) -> tuple[list[list[int]], int]:
    for attempt in range(10_000):
        rng = random.Random(SEED + attempt)
        pool = list(range(LENGTH))
        result = []
        for _ in range(INFORMATION_SETS):
            rng.shuffle(pool)
            selected = []
            selected_values = []
            selected_set = set()
            current_rank = 0
            for coordinate in pool:
                candidate = selected_values + [columns[coordinate]]
                candidate_rank = rank(candidate)
                if candidate_rank == current_rank + 1:
                    selected.append(coordinate)
                    selected_values.append(columns[coordinate])
                    selected_set.add(coordinate)
                    current_rank = candidate_rank
                    if current_rank == DIMENSION:
                        break
            if current_rank != DIMENSION:
                break
            result.append(selected)
            pool = [coordinate for coordinate in pool if coordinate not in selected_set]
        if len(result) == INFORMATION_SETS:
            return result, attempt
    raise RuntimeError("even-zero exact: no information-set packing found")


def main() -> None:
    started = time.perf_counter()
    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(NODES, parity)
    even_nodes = tuple(range(0, NODES, 2))
    odd_nodes = tuple(range(1, NODES, 2))
    even_columns = tuple(select_nodes(row, even_nodes) for row in full_rows)
    state_basis = kernel_basis(even_columns, 64 * len(even_nodes))
    if len(state_basis) != DIMENSION:
        raise RuntimeError("even-zero exact: subcode dimension changed")
    generators = tuple(
        select_nodes(xor_selected(full_rows, state), odd_nodes)
        for state in state_basis
    )
    if rank(list(generators)) != DIMENSION:
        raise RuntimeError("even-zero exact: odd observation is not injective")
    coordinate_columns = [
        sum(
            ((generator >> coordinate) & 1) << bit
            for bit, generator in enumerate(generators)
        )
        for coordinate in range(LENGTH)
    ]
    information_sets, packing_attempt = pack_information_sets(coordinate_columns)
    used = [coordinate for row in information_sets for coordinate in row]
    if len(used) != len(set(used)):
        raise RuntimeError("even-zero exact: information sets overlap")

    expected_per_set = sum(
        math.comb(DIMENSION, weight)
        for weight in range(MAXIMUM_RESTRICTION_WEIGHT + 1)
    )
    expected_candidates = INFORMATION_SETS * expected_per_set
    candidates = 0
    minimum_checked_weight = 1 << 30
    minimum_checked_state = 0
    counterexample = None
    for information_set in information_sets:
        selected_map_columns = tuple(
            sum(
                ((generators[input_bit] >> coordinate) & 1) << position
                for position, coordinate in enumerate(information_set)
            )
            for input_bit in range(DIMENSION)
        )
        inverse = inverse_columns(selected_map_columns)
        systematic_words = tuple(
            xor_selected(generators, message) for message in inverse
        )
        systematic_states = tuple(
            xor_selected(state_basis, message) for message in inverse
        )
        for row, word in enumerate(systematic_words):
            selected = sum(
                ((word >> coordinate) & 1) << position
                for position, coordinate in enumerate(information_set)
            )
            if selected != 1 << row:
                raise RuntimeError("even-zero exact: systematic identity failed")
        for weight in range(MAXIMUM_RESTRICTION_WEIGHT + 1):
            for support in itertools.combinations(range(DIMENSION), weight):
                candidates += 1
                word = 0
                state = 0
                for bit in support:
                    word ^= systematic_words[bit]
                    state ^= systematic_states[bit]
                if state == 0:
                    continue
                exact_weight = word.bit_count()
                if exact_weight < minimum_checked_weight:
                    minimum_checked_weight = exact_weight
                    minimum_checked_state = state
                if exact_weight <= BOUND:
                    counterexample = {
                        "initial_a_hex": f"0x{state & MASK64:016x}",
                        "initial_b_hex": f"0x{state >> 64:016x}",
                        "weight": exact_weight,
                    }
                    break
            if counterexample is not None:
                break
        if counterexample is not None:
            break
    result = "COUNTEREXAMPLE" if counterexample is not None else "EXHAUSTED"
    if counterexample is None and candidates != expected_candidates:
        raise RuntimeError("even-zero exact: candidate count changed")

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal06-even-zero-exact-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_DISJOINT_INFORMATION_SET_ENUMERATION",
        "source_sha256": digest(Path(__file__).resolve()),
        "nodes": NODES,
        "bound": BOUND,
        "forced_zero_nodes": list(even_nodes),
        "retained_nodes": list(odd_nodes),
        "subcode": {"length": LENGTH, "dimension": DIMENSION},
        "information_set_count": INFORMATION_SETS,
        "information_set_size": DIMENSION,
        "information_sets_are_disjoint": True,
        "used_coordinate_count": len(used),
        "unused_coordinate_count": LENGTH - len(used),
        "packing_seed_hex": hex(SEED),
        "packing_attempt": packing_attempt,
        "maximum_enumerated_restriction_weight": MAXIMUM_RESTRICTION_WEIGHT,
        "coverage_inequality": (
            f"{INFORMATION_SETS} * ({MAXIMUM_RESTRICTION_WEIGHT} + 1) "
            f"= {INFORMATION_SETS * (MAXIMUM_RESTRICTION_WEIGHT + 1)} > {BOUND}"
        ),
        "expected_candidates_if_exhausted": expected_candidates,
        "checked_candidates": candidates,
        "minimum_checked_nonzero_weight": minimum_checked_weight,
        "minimum_checked_witness_state": {
            "a_hex": f"0x{minimum_checked_state & MASK64:016x}",
            "b_hex": f"0x{minimum_checked_state >> 64:016x}",
        },
        "result": result,
        "counterexample": counterexample,
        "elapsed_seconds": time.perf_counter() - started,
        "scope_limitation": (
            "The certificate excludes low words only in the subcode whose even "
            "node outputs are all zero. It does not prove D(24) >= 144."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"packing_attempt={packing_attempt}")
    print(f"checked_candidates={candidates}")
    print(f"minimum_checked_nonzero_weight={minimum_checked_weight}")
    print(f"result={result}")
    if counterexample is not None:
        print(f"counterexample={counterexample}")
    print(f"output={OUTPUT}")
    print(f"elapsed_seconds={payload['elapsed_seconds']:.3f}")
    print("status=EXACT_EVEN_ZERO_SUBCODE_ENUMERATION")


if __name__ == "__main__":
    main()
