#!/usr/bin/env python3
"""Independently replay the Goal 02 dimension-at-most-20 response spectrum."""

from __future__ import annotations

import hashlib
import json
import time
from array import array
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
PRIMARY_SOURCE = (
    ROOT
    / "scripts"
    / "certify_riffle_packetmul_wrapmul_2lap_g4_goal02_low_components.cpp"
)
PRIMARY_BINARY = PRIMARY_SOURCE.with_suffix(".exe")
PRIMARY = CANDIDATE / "receipts" / "goal02_low_components_primary.json"
OUTPUT = CANDIDATE / "receipts" / "goal02_low_components_independent.json"

GENERATOR = 0xF4845518B9582A1F
MASK64 = (1 << 64) - 1
RESPONSE_NODES = 32_772
THRESHOLD = 377_530
DIMENSION_CAP = 22
DEGREES = (1, 2, 4, 9, 10, 18, 20)
FACTORS = (0x3, 0x7, 0x13, 0x373, 0x519, 0x7C9C3, 0x1E1FFF)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def apply_columns(columns: tuple[int, ...] | list[int], value: int) -> int:
    result = 0
    while value:
        bit = value & -value
        result ^= columns[bit.bit_length() - 1]
        value ^= bit
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
        raise RuntimeError("low components independent: inverse reconstruction failed")
    return result


def systematic_right_columns() -> tuple[int, ...]:
    rows = tuple((GENERATOR << row) | (1 << 127) for row in range(64))
    left = tuple(row & MASK64 for row in rows)
    right = tuple((row >> 64) & MASK64 for row in rows)
    inverse = inverse_columns(left)
    result = tuple(apply_columns(right, inverse[bit]) for bit in range(64))
    if any(apply_columns(left, inverse[bit]) != 1 << bit for bit in range(64)):
        raise RuntimeError("low components independent: systematic map failed")
    return result


def accumulate_by_bits(value: int) -> int:
    running = 0
    result = 0
    for bit in range(64):
        running ^= (value >> bit) & 1
        result |= running << bit
    return result


def accumulate_fast(value: int) -> int:
    value ^= value << 1
    value ^= value << 2
    value ^= value << 4
    value ^= value << 8
    value ^= value << 16
    value ^= value << 32
    return value & MASK64


def build_apply(columns: tuple[int, ...]):
    tables = []
    for byte_index in range(8):
        table = []
        for byte in range(256):
            table.append(
                apply_columns(columns[8 * byte_index : 8 * byte_index + 8], byte)
            )
        tables.append(tuple(table))

    def apply(value: int) -> int:
        return (
            tables[0][value & 0xFF]
            ^ tables[1][(value >> 8) & 0xFF]
            ^ tables[2][(value >> 16) & 0xFF]
            ^ tables[3][(value >> 24) & 0xFF]
            ^ tables[4][(value >> 32) & 0xFF]
            ^ tables[5][(value >> 40) & 0xFF]
            ^ tables[6][(value >> 48) & 0xFF]
            ^ tables[7][value >> 56]
        )

    return apply


def apply_polynomial(step, polynomial: int, value: int) -> int:
    result = 0
    current = value
    while polynomial:
        if polynomial & 1:
            result ^= current
        polynomial >>= 1
        current = step(current)
    return result


def kernel_basis(step, polynomial: int, expected_dimension: int) -> tuple[int, ...]:
    pivot_values = [0] * 64
    pivot_representations = [0] * 64
    result = []
    for column in range(64):
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
    if len(result) != expected_dimension:
        raise RuntimeError("low components independent: kernel dimension mismatch")
    return tuple(result)


def coordinate_solver(basis: tuple[int, ...]):
    pivot_values = [0] * 64
    pivot_coordinates = [0] * 64
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
            raise RuntimeError("low components independent: dependent basis")

    def coordinates(value: int) -> int:
        result = 0
        while value:
            pivot = value.bit_length() - 1
            if not pivot_values[pivot]:
                raise RuntimeError("low components independent: state outside basis")
            value ^= pivot_values[pivot]
            result ^= pivot_coordinates[pivot]
        return result

    return coordinates


def enumerate_support(
    support_mask: int,
    step,
    component_bases: tuple[tuple[int, ...], ...],
) -> dict:
    basis: list[int] = []
    segment_masks = []
    exact_count = 1
    support_degrees = []
    for component, degree in enumerate(DEGREES):
        if not ((support_mask >> component) & 1):
            continue
        offset = len(basis)
        basis.extend(component_bases[component])
        segment_masks.append(((1 << degree) - 1) << offset)
        exact_count *= (1 << degree) - 1
        support_degrees.append(degree)
    dimension = len(basis)
    size = 1 << dimension
    coordinates = coordinate_solver(tuple(basis))
    next_columns = [coordinates(step(value)) for value in basis]

    next_states = array("I", [0]) * size
    physical_states = array("Q", [0]) * size
    emitted_weights = bytearray(size)
    for index in range(1, size):
        low_bit = index & -index
        bit = low_bit.bit_length() - 1
        rest = index ^ low_bit
        next_states[index] = next_states[rest] ^ next_columns[bit]
        physical_states[index] = physical_states[rest] ^ basis[bit]
        emitted_weights[index] = accumulate_fast(physical_states[index]).bit_count()

    visited = bytearray(size)
    cycle = array("I")
    cycle_count = 0
    visited_count = 0
    low_count = 0
    minimum = 1 << 62
    minimum_index = 0
    minimum_state = 0
    minimum_period = 0
    maximum_absolute_bias = 0
    maximum_bias_index = 0
    maximum_bias_state = 0
    maximum_bias_response_weight = 0
    for start in range(1, size):
        if visited[start] or any((start & mask) == 0 for mask in segment_masks):
            continue
        del cycle[:]
        current = start
        cycle_total = 0
        while True:
            if visited[current] or any(
                (current & mask) == 0 for mask in segment_masks
            ):
                raise RuntimeError(
                    "low components independent: exact-support orbit escaped or merged"
                )
            visited[current] = 1
            cycle.append(current)
            cycle_total += emitted_weights[current]
            current = next_states[current]
            if current == start:
                break
        cycle_count += 1
        visited_count += len(cycle)
        period = len(cycle)
        quotient, remainder = divmod(RESPONSE_NODES, period)
        rolling = sum(emitted_weights[index] for index in cycle[:remainder])
        for cycle_index, state_index in enumerate(cycle):
            response_weight = quotient * cycle_total + rolling
            if response_weight <= THRESHOLD:
                low_count += 1
            if response_weight < minimum:
                minimum = response_weight
                minimum_index = state_index
                minimum_state = physical_states[state_index]
                minimum_period = period
            absolute_bias = abs(RESPONSE_NODES * 64 - 2 * response_weight)
            if absolute_bias > maximum_absolute_bias:
                maximum_absolute_bias = absolute_bias
                maximum_bias_index = state_index
                maximum_bias_state = physical_states[state_index]
                maximum_bias_response_weight = response_weight
            if remainder:
                rolling -= emitted_weights[state_index]
                entering = cycle_index + remainder
                if entering >= period:
                    entering -= period
                rolling += emitted_weights[cycle[entering]]
    if visited_count != exact_count:
        raise RuntimeError("low components independent: exact state count mismatch")
    return {
        "component_support_mask_hex": hex(support_mask),
        "component_degrees": support_degrees,
        "dimension": dimension,
        "exact_nonzero_state_count": exact_count,
        "visited_state_count": visited_count,
        "cycle_count": cycle_count,
        "low_response_state_count": low_count,
        "minimum_response_weight": minimum,
        "minimum_witness_coordinates_hex": hex(minimum_index),
        "minimum_witness_state_hex": hex(minimum_state),
        "minimum_witness_cycle_period": minimum_period,
        "maximum_absolute_response_bias": maximum_absolute_bias,
        "maximum_bias_response_weight": maximum_bias_response_weight,
        "maximum_bias_witness_coordinates_hex": hex(maximum_bias_index),
        "maximum_bias_witness_state_hex": hex(maximum_bias_state),
    }


def main() -> None:
    started = time.perf_counter()
    primary = json.loads(PRIMARY.read_text())
    if primary["source_sha256"] != digest(PRIMARY_SOURCE):
        raise RuntimeError("low components independent: primary source changed")
    if primary["executable_sha256"] != digest(PRIMARY_BINARY):
        raise RuntimeError("low components independent: primary executable changed")

    parity_columns = systematic_right_columns()
    apply_parity = build_apply(parity_columns)
    for bit in range(64):
        if accumulate_by_bits(1 << bit) != accumulate_fast(1 << bit):
            raise RuntimeError("low components independent: accumulation mismatch")
    step_columns = tuple(
        apply_parity(accumulate_by_bits(1 << bit)) for bit in range(64)
    )
    step = build_apply(step_columns)
    component_bases = tuple(
        kernel_basis(step, factor, degree)
        for factor, degree in zip(FACTORS, DEGREES, strict=True)
    )
    all_coordinates = coordinate_solver(tuple(sum((list(row) for row in component_bases), [])))
    for bit in range(64):
        all_coordinates(step_columns[bit])

    rows = []
    for support_mask in range(1, 1 << len(DEGREES)):
        if sum(
            degree
            for component, degree in enumerate(DEGREES)
            if (support_mask >> component) & 1
        ) <= DIMENSION_CAP:
            row = enumerate_support(support_mask, step, component_bases)
            rows.append(row)
            print(
                f"support={row['component_support_mask_hex']} "
                f"dimension={row['dimension']} minimum={row['minimum_response_weight']} "
                f"low={row['low_response_state_count']}"
            )

    primary_rows = primary["rows"]
    if len(rows) != len(primary_rows):
        raise RuntimeError("low components independent: support count differs")
    comparison_fields = (
        "component_support_mask_hex",
        "component_degrees",
        "dimension",
        "exact_nonzero_state_count",
        "visited_state_count",
        "cycle_count",
        "low_response_state_count",
        "minimum_response_weight",
        "minimum_witness_state_hex",
        "minimum_witness_cycle_period",
        "maximum_absolute_response_bias",
        "maximum_bias_response_weight",
        "maximum_bias_witness_state_hex",
    )
    for actual, expected in zip(rows, primary_rows, strict=True):
        for field in comparison_fields:
            if actual[field] != expected[field]:
                raise RuntimeError(
                    f"low components independent: row field differs: {field}"
                )
    total_states = sum(row["exact_nonzero_state_count"] for row in rows)
    total_low = sum(row["low_response_state_count"] for row in rows)
    global_minimum = min(row["minimum_response_weight"] for row in rows)
    global_maximum_absolute_bias = max(
        row["maximum_absolute_response_bias"] for row in rows
    )
    if (
        total_states != primary["total_exact_nonzero_states"]
        or total_low != primary["total_low_response_states"]
        or global_minimum != primary["global_minimum_response_weight"]
        or global_maximum_absolute_bias
        != primary["global_maximum_absolute_response_bias"]
    ):
        raise RuntimeError("low components independent: aggregate differs")

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-low-components-independent-v1",
        "candidate": primary["candidate"],
        "evidence_label": "INDEPENDENT_EXACT_EXHAUSTIVE_REPLAY",
        "source_sha256": digest(Path(__file__).resolve()),
        "primary_source_sha256": digest(PRIMARY_SOURCE),
        "primary_executable_sha256": digest(PRIMARY_BINARY),
        "primary_receipt_sha256": digest(PRIMARY),
        "recurrence_reconstruction": {
            "systematic_map_method": "binary Gaussian elimination",
            "accumulation_method": "bit-by-bit prefix parity",
            "component_kernel_dimensions": list(DEGREES),
            "direct_sum_dimension": 64,
        },
        "response_nodes": RESPONSE_NODES,
        "low_weight_threshold": THRESHOLD,
        "component_dimension_cap": DIMENSION_CAP,
        "support_count": len(rows),
        "total_exact_nonzero_states": total_states,
        "total_low_response_states": total_low,
        "global_minimum_response_weight": global_minimum,
        "global_maximum_absolute_response_bias": global_maximum_absolute_bias,
        "rows_replayed": len(rows),
        "elapsed_seconds": time.perf_counter() - started,
        "result": "PASS" if total_low == 0 else "LOW_RESPONSE_STATES_FOUND",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print(f"total_states={total_states}")
    print(f"global_minimum={global_minimum}")
    print(f"global_maximum_absolute_bias={global_maximum_absolute_bias}")
    print(f"status={payload['result']}")


if __name__ == "__main__":
    main()
