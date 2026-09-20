#!/usr/bin/env python3
"""Exact Goal 05 audit on low-dimensional lifted primary supports.

The lifted zero-input recurrence has two copies of the 64-bit autonomous map on
its diagonal.  This script constructs the primary spaces ker(f_i(R)^2), verifies
that they form a direct sum, and exhausts every exact primary support whose
dimension is at most 22.  For each state it measures the two adjacent nine-node
windows, hence the 18-node return cost.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    lifted_step,
    observation_rows,
)
from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
OUTPUT = CANDIDATE / "receipts" / "goal05_primary_components.json"

DEGREES = (1, 2, 4, 9, 10, 18, 20)
FACTORS = (0x3, 0x7, 0x13, 0x373, 0x519, 0x7C9C3, 0x1E1FFF)
DIMENSION_CAP = 22
WINDOW_NODES = 9
PAIR_NODES = 18
LOW_WINDOW_MAXIMUM = 52
RETURN_TARGET = 106
MASK64 = (1 << 64) - 1
FIRST_WINDOW_MASK = (1 << (64 * WINDOW_NODES)) - 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lifted_state_step(state: int, parity) -> int:
    a = state & MASK64
    b = state >> 64
    next_a, next_b, _ = lifted_step(a, b, parity)
    return next_a | (next_b << 64)


def apply_polynomial(step, polynomial: int, value: int) -> int:
    result = 0
    current = value
    while polynomial:
        if polynomial & 1:
            result ^= current
        polynomial >>= 1
        current = step(current)
    return result


def polynomial_square(polynomial: int) -> int:
    result = 0
    exponent = 0
    while polynomial:
        if polynomial & 1:
            result |= 1 << (2 * exponent)
        polynomial >>= 1
        exponent += 1
    return result


def kernel_basis(step, polynomial: int, width: int) -> tuple[int, ...]:
    """Return dependency vectors among the columns of polynomial(step)."""
    pivot_values = [0] * width
    pivot_representations = [0] * width
    result: list[int] = []
    for column in range(width):
        value = apply_polynomial(step, polynomial, 1 << column)
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


def coordinate_solver(basis: tuple[int, ...], width: int):
    pivot_values = [0] * width
    pivot_coordinates = [0] * width
    for coordinate, original in enumerate(basis):
        value = original
        representation = 1 << coordinate
        while value:
            pivot = value.bit_length() - 1
            if pivot_values[pivot]:
                value ^= pivot_values[pivot]
                representation ^= pivot_coordinates[pivot]
            else:
                pivot_values[pivot] = value
                pivot_coordinates[pivot] = representation
                break
        if value == 0:
            raise RuntimeError("goal05 components: dependent primary basis")

    def coordinates(value: int) -> int:
        result = 0
        while value:
            pivot = value.bit_length() - 1
            if not pivot_values[pivot]:
                raise RuntimeError("goal05 components: state outside primary basis")
            value ^= pivot_values[pivot]
            result ^= pivot_coordinates[pivot]
        return result

    return coordinates


def xor_selected(rows: list[int] | tuple[int, ...], selector: int) -> int:
    result = 0
    while selector:
        bit = selector & -selector
        result ^= rows[bit.bit_length() - 1]
        selector ^= bit
    return result


def enumerate_support(
    support_mask: int,
    component_bases: tuple[tuple[int, ...], ...],
    physical_observation_rows: list[int],
) -> dict:
    basis: list[int] = []
    segment_masks: list[int] = []
    component_degrees: list[int] = []
    exact_count = 1
    for component, component_basis in enumerate(component_bases):
        if not ((support_mask >> component) & 1):
            continue
        offset = len(basis)
        basis.extend(component_basis)
        dimension = len(component_basis)
        segment_masks.append(((1 << dimension) - 1) << offset)
        component_degrees.append(DEGREES[component])
        exact_count *= (1 << dimension) - 1

    observation_basis = [
        xor_selected(physical_observation_rows, state) for state in basis
    ]
    dimension = len(basis)
    previous_gray = 0
    state = 0
    word = 0
    visited = 0
    low_first_count = 0
    violating_count = 0
    minimum_total = 1 << 30
    minimum_total_state = 0
    minimum_first = 1 << 30
    minimum_first_state = 0
    minimum_second_given_low_first = 1 << 30
    minimum_second_given_low_first_state = 0
    minimum_return_given_low_first = 1 << 30
    minimum_return_given_low_first_state = 0

    for index in range(1, 1 << dimension):
        gray = index ^ (index >> 1)
        changed = gray ^ previous_gray
        bit = changed.bit_length() - 1
        state ^= basis[bit]
        word ^= observation_basis[bit]
        previous_gray = gray
        if any((gray & segment_mask) == 0 for segment_mask in segment_masks):
            continue
        visited += 1
        first = (word & FIRST_WINDOW_MASK).bit_count()
        second = (word >> (64 * WINDOW_NODES)).bit_count()
        total = first + second
        if total < minimum_total:
            minimum_total = total
            minimum_total_state = state
        if first < minimum_first:
            minimum_first = first
            minimum_first_state = state
        if first <= LOW_WINDOW_MAXIMUM:
            low_first_count += 1
            if second < minimum_second_given_low_first:
                minimum_second_given_low_first = second
                minimum_second_given_low_first_state = state
            if total < minimum_return_given_low_first:
                minimum_return_given_low_first = total
                minimum_return_given_low_first_state = state
            if total < RETURN_TARGET:
                violating_count += 1

    if visited != exact_count:
        raise RuntimeError("goal05 components: exact-support state count mismatch")

    def witness(value: int) -> dict[str, str]:
        return {
            "a_hex": f"0x{value & MASK64:016x}",
            "b_hex": f"0x{value >> 64:016x}",
        }

    row = {
        "component_support_mask_hex": hex(support_mask),
        "component_degrees": component_degrees,
        "lifted_dimension": dimension,
        "exact_nonzero_state_count": exact_count,
        "minimum_18_node_weight": minimum_total,
        "minimum_18_node_witness": witness(minimum_total_state),
        "minimum_first_9_node_weight": minimum_first,
        "minimum_first_9_node_witness": witness(minimum_first_state),
        "first_9_weight_at_most_52_count": low_first_count,
        "paired_return_violating_state_count": violating_count,
    }
    if low_first_count:
        row.update(
            {
                "minimum_second_9_weight_given_first_at_most_52": (
                    minimum_second_given_low_first
                ),
                "minimum_second_given_low_first_witness": witness(
                    minimum_second_given_low_first_state
                ),
                "minimum_18_weight_given_first_at_most_52": (
                    minimum_return_given_low_first
                ),
                "minimum_return_given_low_first_witness": witness(
                    minimum_return_given_low_first_state
                ),
            }
        )
    return row


def component_support(
    state: int,
    all_coordinates,
    component_bases: tuple[tuple[int, ...], ...],
) -> dict:
    coordinates = all_coordinates(state)
    offset = 0
    support_mask = 0
    coordinate_chunks = []
    for component, basis in enumerate(component_bases):
        chunk = (coordinates >> offset) & ((1 << len(basis)) - 1)
        coordinate_chunks.append(hex(chunk))
        if chunk:
            support_mask |= 1 << component
        offset += len(basis)
    return {
        "component_support_mask_hex": hex(support_mask),
        "coordinate_chunks_hex": coordinate_chunks,
        "support_lifted_dimension": sum(
            len(component_bases[component])
            for component in range(len(component_bases))
            if (support_mask >> component) & 1
        ),
    }


def main() -> None:
    started = time.perf_counter()
    source = Path(__file__).resolve()
    parity = build_apply(systematic_state_columns())
    step = lambda state: lifted_state_step(state, parity)

    component_bases = tuple(
        kernel_basis(step, polynomial_square(factor), 128) for factor in FACTORS
    )
    dimensions = tuple(len(basis) for basis in component_bases)
    expected_dimensions = tuple(2 * degree for degree in DEGREES)
    if dimensions != expected_dimensions:
        raise RuntimeError(
            f"goal05 components: dimensions {dimensions} != {expected_dimensions}"
        )

    all_basis = tuple(value for basis in component_bases for value in basis)
    if len(all_basis) != 128:
        raise RuntimeError("goal05 components: primary dimensions do not sum to 128")
    all_coordinates = coordinate_solver(all_basis, 128)
    for bit in range(128):
        all_coordinates(step(1 << bit))
    for component, (factor, basis) in enumerate(
        zip(FACTORS, component_bases, strict=True)
    ):
        annihilator = polynomial_square(factor)
        if any(apply_polynomial(step, annihilator, value) for value in basis):
            raise RuntimeError(
                f"goal05 components: component {component} annihilation failed"
            )
        if any(
            component_support(step(value), all_coordinates, component_bases)[
                "component_support_mask_hex"
            ]
            not in ("0x0", hex(1 << component))
            for value in basis
        ):
            raise RuntimeError(
                f"goal05 components: component {component} is not invariant"
            )

    physical_observation_rows = observation_rows(PAIR_NODES, parity)
    known_state = int("961e658de73963d2", 16) << 64 | int(
        "ba22ae96294ba476", 16
    )
    known_word = xor_selected(physical_observation_rows, known_state)
    known_first = (known_word & FIRST_WINDOW_MASK).bit_count()
    known_second = (known_word >> (64 * WINDOW_NODES)).bit_count()
    if (known_first, known_second) != (44, 275):
        raise RuntimeError("goal05 components: known witness replay changed")

    rows = []
    for support_mask in range(1, 1 << len(DEGREES)):
        dimension = sum(
            dimensions[component]
            for component in range(len(DEGREES))
            if (support_mask >> component) & 1
        )
        if dimension > DIMENSION_CAP:
            continue
        row = enumerate_support(
            support_mask, component_bases, physical_observation_rows
        )
        rows.append(row)
        print(
            f"support={row['component_support_mask_hex']} "
            f"dim={row['lifted_dimension']} min18={row['minimum_18_node_weight']} "
            f"low9={row['first_9_weight_at_most_52_count']} "
            f"violations={row['paired_return_violating_state_count']}"
        )

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal05-components-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_LOW_DIMENSIONAL_PRIMARY_SUPPORT_ENUMERATION",
        "source_sha256": digest(source),
        "autonomous_irreducible_degrees": list(DEGREES),
        "autonomous_irreducible_factors_hex": [hex(value) for value in FACTORS],
        "lifted_primary_dimensions": list(dimensions),
        "dimension_cap": DIMENSION_CAP,
        "window_nodes": WINDOW_NODES,
        "low_window_maximum": LOW_WINDOW_MAXIMUM,
        "paired_return_target": RETURN_TARGET,
        "known_weight_44_witness": {
            "a_hex": f"0x{known_state & MASK64:016x}",
            "b_hex": f"0x{known_state >> 64:016x}",
            "first_9_weight": known_first,
            "second_9_weight": known_second,
            **component_support(known_state, all_coordinates, component_bases),
        },
        "rows": rows,
        "enumerated_exact_nonzero_states": sum(
            row["exact_nonzero_state_count"] for row in rows
        ),
        "enumerated_low_first_window_states": sum(
            row["first_9_weight_at_most_52_count"] for row in rows
        ),
        "enumerated_paired_return_violations": sum(
            row["paired_return_violating_state_count"] for row in rows
        ),
        "global_minimum_18_node_weight": min(
            row["minimum_18_node_weight"] for row in rows
        ),
        "elapsed_seconds": time.perf_counter() - started,
        "scope_limitation": (
            "This is exhaustive only for exact primary supports of total lifted "
            "dimension at most 22; larger and mixed high-degree supports remain."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print(f"elapsed_seconds={payload['elapsed_seconds']:.3f}")
    print("status=EXACT_LOW_DIMENSIONAL_PRIMARY_SUPPORT_ENUMERATION")


if __name__ == "__main__":
    main()
