#!/usr/bin/env python3
"""Enumerate the 24-dimensional even trace after a zero anchor."""

from __future__ import annotations

import hashlib
import json
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


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
OUTPUT = CANDIDATE / "receipts" / "goal05_zero_anchor_even_trace.json"
ANCHOR_NODE = 0
EVEN_NODES = tuple(range(2, 18, 2))
ODD_NODES = tuple(range(1, 18, 2))
EVEN_DIMENSION = 24
KERNEL_DIMENSION = 40
LOW_TOTAL_MAXIMUM = 105


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def select_nodes(word: int, nodes: tuple[int, ...]) -> int:
    result = 0
    for output_node, input_node in enumerate(nodes):
        result |= ((word >> (64 * input_node)) & MASK64) << (64 * output_node)
    return result


def image_basis(
    image_columns: tuple[int, ...], state_basis: tuple[int, ...]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    pivots: dict[int, int] = {}
    image_result = []
    state_result = []
    for image, state in zip(image_columns, state_basis, strict=True):
        reduced = image
        reduced_state = state
        while reduced:
            pivot = reduced.bit_length() - 1
            if pivot in pivots:
                owner = pivots[pivot]
                reduced ^= image_result[owner]
                reduced_state ^= state_result[owner]
            else:
                pivots[pivot] = len(image_result)
                image_result.append(reduced)
                state_result.append(reduced_state)
                break
    return tuple(image_result), tuple(state_result)


def main() -> None:
    started = time.perf_counter()
    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(18, parity)
    anchor_columns = tuple(row & MASK64 for row in full_rows)
    zero_anchor_basis = kernel_basis(anchor_columns, 64)
    if len(zero_anchor_basis) != 64:
        raise RuntimeError("even trace: zero-anchor dimension changed")

    full_words = tuple(xor_selected(full_rows, state) for state in zero_anchor_basis)
    even_columns = tuple(select_nodes(word, EVEN_NODES) for word in full_words)
    odd_columns = tuple(select_nodes(word, ODD_NODES) for word in full_words)
    even_basis, quotient_state_basis = image_basis(even_columns, zero_anchor_basis)
    if len(even_basis) != EVEN_DIMENSION:
        raise RuntimeError(f"even trace: rank {len(even_basis)} is not 24")
    even_kernel_coordinates = kernel_basis(even_columns, 64 * len(EVEN_NODES))
    if len(even_kernel_coordinates) != KERNEL_DIMENSION:
        raise RuntimeError("even trace: kernel dimension is not 40")
    even_kernel_state_basis = tuple(
        xor_selected(zero_anchor_basis, coordinates)
        for coordinates in even_kernel_coordinates
    )
    if any(
        select_nodes(xor_selected(full_rows, state), EVEN_NODES)
        for state in even_kernel_state_basis
    ):
        raise RuntimeError("even trace: reconstructed kernel emits an even trace")

    histogram = [0] * (64 * len(EVEN_NODES) + 1)
    low_count = 0
    minimum_nonzero = 1 << 30
    minimum_state = 0
    previous_gray = 0
    word = 0
    state = 0
    histogram[0] = 1
    for index in range(1, 1 << EVEN_DIMENSION):
        gray = index ^ (index >> 1)
        changed = gray ^ previous_gray
        bit = changed.bit_length() - 1
        word ^= even_basis[bit]
        state ^= quotient_state_basis[bit]
        previous_gray = gray
        weight = word.bit_count()
        histogram[weight] += 1
        if weight <= LOW_TOTAL_MAXIMUM:
            low_count += 1
        if weight < minimum_nonzero:
            minimum_nonzero = weight
            minimum_state = state

    if sum(histogram) != 1 << EVEN_DIMENSION:
        raise RuntimeError("even trace: histogram mass changed")
    minimum_word = select_nodes(xor_selected(full_rows, minimum_state), EVEN_NODES)
    minimum_node_weights = [
        ((minimum_word >> (64 * node)) & MASK64).bit_count()
        for node in range(len(EVEN_NODES))
    ]
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal05-zero-anchor-even-trace-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_COMPLETE_EVEN_TRACE_ENUMERATION",
        "source_sha256": digest(Path(__file__).resolve()),
        "anchor_node": ANCHOR_NODE,
        "even_nodes": list(EVEN_NODES),
        "odd_nodes": list(ODD_NODES),
        "zero_anchor_dimension": len(zero_anchor_basis),
        "even_trace_dimension": len(even_basis),
        "even_trace_kernel_dimension": len(even_kernel_state_basis),
        "enumerated_even_traces": 1 << EVEN_DIMENSION,
        "minimum_nonzero_even_trace_weight": minimum_nonzero,
        "minimum_witness": {
            "initial_a_hex": f"0x{minimum_state & MASK64:016x}",
            "initial_b_hex": f"0x{minimum_state >> 64:016x}",
            "even_node_weights": minimum_node_weights,
        },
        "even_trace_weight_at_most_105_count_including_zero": low_count + 1,
        "weight_histogram": {
            str(weight): count for weight, count in enumerate(histogram) if count
        },
        "elapsed_seconds": time.perf_counter() - started,
        "remaining_problem": (
            "For each even trace of weight e at most 105, exclude an odd-trace "
            "word of weight at most 105-e in its 40-dimensional affine coset."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"even_trace_dimension={payload['even_trace_dimension']}")
    print(f"kernel_dimension={payload['even_trace_kernel_dimension']}")
    print(f"minimum_nonzero_even_weight={minimum_nonzero}")
    print(
        "low_even_trace_count="
        f"{payload['even_trace_weight_at_most_105_count_including_zero']}"
    )
    print(f"output={OUTPUT}")
    print(f"elapsed_seconds={payload['elapsed_seconds']:.3f}")
    print("status=EXACT_COMPLETE_EVEN_TRACE_ENUMERATION")


if __name__ == "__main__":
    main()
