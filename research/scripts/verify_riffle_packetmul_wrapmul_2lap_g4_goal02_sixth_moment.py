#!/usr/bin/env python3
"""Independently replay the Goal 02 sixth-moment reduction."""

from __future__ import annotations

import hashlib
import json
import struct
from pathlib import Path

from verify_riffle_packetmul_wrapmul_2lap_g4_goal02_low_components import (
    accumulate_by_bits,
    build_apply,
    systematic_right_columns,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
PRIMARY_SOURCE = (
    ROOT / "scripts" / "certify_riffle_packetmul_wrapmul_2lap_g4_goal02_sixth_moment.py"
)
RECONSTRUCTION_SOURCE = (
    ROOT / "scripts" / "verify_riffle_packetmul_wrapmul_2lap_g4_goal02_low_components.py"
)
PRIMARY = CANDIDATE / "receipts" / "goal02_sixth_moment_primary.json"
OUTPUT = CANDIDATE / "receipts" / "goal02_sixth_moment_independent.json"
NODES = 32_772
LENGTH = 64 * NODES
THRESHOLD = 377_530
NONZERO_TARGET = 316_605


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def transpose_by_rows(columns: tuple[int, ...]) -> tuple[int, ...]:
    rows = []
    for output in range(64):
        row = 0
        for input_bit, column in enumerate(columns):
            row |= ((column >> output) & 1) << input_bit
        rows.append(row)
    return tuple(rows)


def main() -> None:
    primary = json.loads(PRIMARY.read_text())
    if primary["source_sha256"] != digest(PRIMARY_SOURCE):
        raise RuntimeError("sixth moment independent: primary source changed")
    parity = systematic_right_columns()
    apply_parity = build_apply(parity)
    output_step_columns = tuple(
        accumulate_by_bits(apply_parity(1 << bit)) for bit in range(64)
    )
    transpose_step = build_apply(transpose_by_rows(output_step_columns))
    current = [1 << bit for bit in range(64)]
    seen = set()
    sequence_hash = hashlib.sha256()
    witness_indices = tuple(
        primary["authenticated_weight4_dual_witness"]["flat_coordinate_indices"]
    )
    witness = {}
    flat_index = 0
    for _node in range(NODES):
        for value in current:
            if value == 0 or value in seen:
                raise RuntimeError("sixth moment independent: zero or repeated coordinate")
            seen.add(value)
            sequence_hash.update(struct.pack("<Q", value))
            if flat_index in witness_indices:
                witness[flat_index] = value
            flat_index += 1
        current = [transpose_step(value) for value in current]
    witness_xor = 0
    for index in witness_indices:
        witness_xor ^= witness[index]
    if witness_xor:
        raise RuntimeError("sixth moment independent: dual witness failed")
    census = primary["coordinate_census"]
    if (
        len(seen) != LENGTH
        or census["distinct_coordinate_count"] != len(seen)
        or census["ordered_coordinate_sequence_sha256"] != sequence_hash.hexdigest()
    ):
        raise RuntimeError("sixth moment independent: coordinate census differs")

    deviation = LENGTH // 2 - THRESHOLD
    baseline = LENGTH + 30 * (LENGTH * (LENGTH - 1) // 2) + 90 * (
        LENGTH * (LENGTH - 1) * (LENGTH - 2) // 6
    )
    coefficient4 = 4 * (720 // 6) + (LENGTH - 4) * (720 // 2)
    sparse_budget = (
        (NONZERO_TARGET + 1) * 64 * deviation**6 - (1 << 64) * baseline
    )
    moment = primary["sixth_moment_identity"]
    reduction = primary["markov_reduction"]
    if (
        moment["baseline_even_multiplicity_term"] != baseline
        or moment["weight4_dual_coefficient"] != coefficient4
        or reduction["scaled_joint_sparse_budget"] != sparse_budget
    ):
        raise RuntimeError("sixth moment independent: arithmetic differs")
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-sixth-moment-independent-v1",
        "candidate": primary["candidate"],
        "evidence_label": "INDEPENDENT_EXACT_RECONSTRUCTION",
        "source_sha256": digest(Path(__file__).resolve()),
        "reconstruction_source_sha256": digest(RECONSTRUCTION_SOURCE),
        "primary_source_sha256": digest(PRIMARY_SOURCE),
        "primary_receipt_sha256": digest(PRIMARY),
        "systematic_map_method": "binary Gaussian elimination",
        "accumulation_method": "bit-by-bit prefix parity",
        "distinct_coordinate_count": len(seen),
        "ordered_coordinate_sequence_sha256": sequence_hash.hexdigest(),
        "weight4_witness_xor_hex": hex(witness_xor),
        "baseline_even_multiplicity_term": baseline,
        "weight4_dual_coefficient": coefficient4,
        "weight6_dual_coefficient": 720,
        "scaled_joint_sparse_budget": sparse_budget,
        "result": "PASS; SPARSE_DUAL_COUNTS_REMAIN",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print(f"distinct_coordinates={len(seen)}")
    print(f"weight4_witness_xor={hex(witness_xor)}")
    print("status=INDEPENDENT_EXACT_RECONSTRUCTION")


if __name__ == "__main__":
    main()
