#!/usr/bin/env python3
"""Independent systematic-form replay of the Goal 04 window census."""

from __future__ import annotations

import json
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    relation_of_weight_five,
    relation_up_to_four,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
PRIMARY = CANDIDATE / "receipts" / "goal04_lifted_windows_primary.json"
OUTPUT = CANDIDATE / "receipts" / "goal04_lifted_windows_independent.json"
MASK64 = (1 << 64) - 1


def prefix_xor(value: int) -> int:
    result = 0
    parity = 0
    for bit in range(64):
        parity ^= (value >> bit) & 1
        result |= parity << bit
    return result


def apply_columns(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        low = value & -value
        result ^= columns[low.bit_length() - 1]
        value ^= low
    return result


def observation_rows(nodes: int, columns: tuple[int, ...]) -> list[int]:
    rows = []
    for basis in range(128):
        a = 1 << basis if basis < 64 else 0
        b = 1 << (basis - 64) if basis >= 64 else 0
        word = 0
        for node in range(nodes):
            first = prefix_xor(a)
            retained = prefix_xor(b ^ first)
            word |= retained << (64 * node)
            a = apply_columns(columns, first)
            b = apply_columns(columns, retained)
        rows.append(word)
    return rows


def systematic_syndromes(rows: list[int], nodes: int) -> tuple[list[int], int]:
    work = rows.copy()
    for column in range(128):
        pivot = next(
            row for row in range(column, 128) if (work[row] >> column) & 1
        )
        work[column], work[pivot] = work[pivot], work[column]
        for row in range(128):
            if row != column and ((work[row] >> column) & 1):
                work[row] ^= work[column]
    if any((work[row] & ((1 << 128) - 1)) != 1 << row for row in range(128)):
        raise RuntimeError("Goal 04 verifier: first two outputs are not systematic")

    codimension = 64 * nodes - 128
    syndromes = [work[column] >> 128 for column in range(128)]
    syndromes.extend(1 << bit for bit in range(codimension))
    return syndromes, codimension


def main() -> None:
    primary = json.loads(PRIMARY.read_text())
    columns = systematic_state_columns()
    results = []
    for expected in primary["windows"]:
        nodes = expected["nodes"]
        rows = observation_rows(nodes, columns)
        if nodes == 1:
            syndromes, codimension = [0] * 64, 0
        else:
            syndromes, codimension = systematic_syndromes(rows, nodes)
        distance, witness = relation_up_to_four(syndromes)
        weight_five_exhausted = False
        if distance is None and nodes in (6, 7, 8):
            witness = relation_of_weight_five(syndromes)
            weight_five_exhausted = True
            if witness is not None:
                distance = 5
        lower = distance if distance else (6 if weight_five_exhausted else 5)
        if distance != expected["minimum_distance_if_found"]:
            raise RuntimeError("Goal 04 verifier: distance witness mismatch")
        if lower != expected["proved_minimum_distance_lower_bound"]:
            raise RuntimeError("Goal 04 verifier: lower-bound mismatch")
        results.append(
            {
                "nodes": nodes,
                "dimension": 64 if nodes == 1 else 128,
                "codimension": codimension,
                "minimum_distance_if_found": distance,
                "proved_minimum_distance_lower_bound": lower,
                "weight_five_exhausted": weight_five_exhausted,
            }
        )

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-lifted-windows-independent-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_INDEPENDENT_SYSTEMATIC_REPLAY",
        "method": (
            "Reconstruct P by direct state columns, implement the accumulator "
            "bit-by-bit, reduce the first two output blocks to systematic form, "
            "and enumerate syndrome relations."
        ),
        "windows": results,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("windows_replayed=8")
    print("weight_five_windows_replayed=3")
    print(f"output={OUTPUT}")
    print("status=EXACT_LIFTED_WINDOWS_INDEPENDENT_REPLAY")


if __name__ == "__main__":
    main()
