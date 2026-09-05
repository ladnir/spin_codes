#!/usr/bin/env python3
"""Exact small-weight audit of lifted two-lap observation windows."""

from __future__ import annotations

import json
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "constructions"
    / "riffle_packetmul_wrapmul_2lap_g4"
    / "receipts"
    / "goal04_lifted_windows_primary.json"
)
MASK64 = (1 << 64) - 1
WINDOWS = tuple(range(1, 9))


def accumulate(value: int) -> int:
    for shift in (1, 2, 4, 8, 16, 32):
        value ^= value << shift
    return value & MASK64


def build_apply(columns: tuple[int, ...]):
    tables: list[list[int]] = []
    for byte_index in range(8):
        table: list[int] = []
        for byte in range(256):
            image = 0
            for bit in range(8):
                if (byte >> bit) & 1:
                    image ^= columns[8 * byte_index + bit]
            table.append(image)
        tables.append(table)

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


def lifted_step(a: int, b: int, parity) -> tuple[int, int, int]:
    first_output = accumulate(a)
    retained_output = accumulate(b ^ first_output)
    return parity(first_output), parity(retained_output), retained_output


def observation_rows(nodes: int, parity) -> list[int]:
    rows: list[int] = []
    for basis in range(128):
        a = 1 << basis if basis < 64 else 0
        b = 1 << (basis - 64) if basis >= 64 else 0
        word = 0
        for node in range(nodes):
            a, b, output = lifted_step(a, b, parity)
            word |= output << (64 * node)
        rows.append(word)
    return rows


def rref(rows: list[int], columns: int) -> tuple[list[int], list[int]]:
    work = [row for row in rows if row]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(columns):
        selected = next(
            (row for row in range(pivot_row, len(work)) if (work[row] >> column) & 1),
            None,
        )
        if selected is None:
            continue
        work[pivot_row], work[selected] = work[selected], work[pivot_row]
        pivot = work[pivot_row]
        for row in range(len(work)):
            if row != pivot_row and ((work[row] >> column) & 1):
                work[row] ^= pivot
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(work):
            break
    return work[:pivot_row], pivots


def parity_checks(generator: list[int], columns: int) -> tuple[list[int], list[int]]:
    reduced, pivots = rref(generator, columns)
    pivot_set = set(pivots)
    checks: list[int] = []
    for free in range(columns):
        if free in pivot_set:
            continue
        check = 1 << free
        for row, pivot in zip(reduced, pivots):
            if (row >> free) & 1:
                check |= 1 << pivot
        checks.append(check)
    if any((row & check).bit_count() & 1 for row in generator for check in checks):
        raise RuntimeError("lifted windows: parity-check construction failed")
    return checks, pivots


def coordinate_syndromes(checks: list[int], columns: int) -> list[int]:
    return [
        sum(((check >> column) & 1) << row for row, check in enumerate(checks))
        for column in range(columns)
    ]


def relation_up_to_four(syndromes: list[int]) -> tuple[int | None, tuple[int, ...] | None]:
    zero = next((index for index, value in enumerate(syndromes) if value == 0), None)
    if zero is not None:
        return 1, (zero,)

    singletons: dict[int, int] = {}
    for index, value in enumerate(syndromes):
        old = singletons.get(value)
        if old is not None:
            return 2, (old, index)
        singletons[value] = index

    pair_owner: dict[int, tuple[int, int]] = {}
    for left in range(len(syndromes)):
        for right in range(left + 1, len(syndromes)):
            value = syndromes[left] ^ syndromes[right]
            singleton = singletons.get(value)
            if singleton is not None and singleton not in (left, right):
                return 3, (left, right, singleton)
            old = pair_owner.get(value)
            if old is not None:
                if len({old[0], old[1], left, right}) == 4:
                    return 4, (old[0], old[1], left, right)
            else:
                pair_owner[value] = (left, right)
    return None, None


def relation_of_weight_five(syndromes: list[int]) -> tuple[int, ...] | None:
    pair_owner: dict[int, tuple[int, int]] = {}
    for left in range(len(syndromes)):
        for right in range(left + 1, len(syndromes)):
            value = syndromes[left] ^ syndromes[right]
            if value in pair_owner:
                raise RuntimeError(
                    "lifted windows: weight-four relation escaped the prior search"
                )
            pair_owner[value] = (left, right)

    for first in range(len(syndromes)):
        first_value = syndromes[first]
        for second in range(first + 1, len(syndromes)):
            pair_value = first_value ^ syndromes[second]
            for third in range(second + 1, len(syndromes)):
                owner = pair_owner.get(pair_value ^ syndromes[third])
                if owner is not None and not {
                    first,
                    second,
                    third,
                }.intersection(owner):
                    return (owner[0], owner[1], first, second, third)
    return None


def serialize_witness(indices: tuple[int, ...] | None) -> list[dict[str, int]]:
    if indices is None:
        return []
    return [
        {"coordinate": index, "node": index // 64, "bit": index % 64}
        for index in indices
    ]


def main() -> None:
    parity = build_apply(systematic_state_columns())
    rows_by_window = {nodes: observation_rows(nodes, parity) for nodes in WINDOWS}
    results = []
    for nodes, generator in rows_by_window.items():
        columns = 64 * nodes
        checks, pivots = parity_checks(generator, columns)
        syndromes = coordinate_syndromes(checks, columns)
        distance, witness = relation_up_to_four(syndromes)
        weight_five_exhausted = False
        if distance is None and nodes in (6, 7, 8):
            witness = relation_of_weight_five(syndromes)
            weight_five_exhausted = True
            if witness is not None:
                distance = 5
        if witness is not None:
            word = sum(1 << coordinate for coordinate in witness)
            if any((check & word).bit_count() & 1 for check in checks):
                raise RuntimeError("lifted windows: reported witness is not a codeword")
        results.append(
            {
                "nodes": nodes,
                "length": columns,
                "dimension": len(pivots),
                "codimension": len(checks),
                "minimum_distance_if_found": distance,
                "proved_minimum_distance_lower_bound": (
                    distance if distance else (6 if weight_five_exhausted else 5)
                ),
                "weight_five_exhausted": weight_five_exhausted,
                "witness": serialize_witness(witness),
            }
        )

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-lifted-windows-primary-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_SMALL_WEIGHT_ENUMERATION",
        "windows": results,
        "search_scope": "All codewords of weights one through four.",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for result in results:
        print(
            f"nodes={result['nodes']} dimension={result['dimension']} "
            f"distance_found={result['minimum_distance_if_found']}"
        )
    print(f"output={OUTPUT}")
    print("status=EXACT_LIFTED_WINDOW_SMALL_WEIGHT_ENUMERATION")


if __name__ == "__main__":
    main()
