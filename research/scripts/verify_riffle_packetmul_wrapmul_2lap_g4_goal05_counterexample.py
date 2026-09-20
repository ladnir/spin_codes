#!/usr/bin/env python3
"""Independently verify the Goal 05 paired-return counterexample."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
RAW = CANDIDATE / "receipts" / "goal05_zero_anchor_certificate_raw.json"
OUTPUT = CANDIDATE / "receipts" / "goal05_return_counterexample_verified.json"
GENERATOR = 0xF4845518B9582A1F
MASK64 = (1 << 64) - 1
EXPECTED_A = 0x89CF7ABE4FBC2365
EXPECTED_B = 0x87452995C56BE123
EXPECTED_WEIGHTS = (0, 7, 0, 20, 0, 10, 0, 10, 0, 5, 0, 18, 0, 9, 0, 11, 0, 7)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def apply_columns(columns: tuple[int, ...], value: int) -> int:
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
        raise RuntimeError("counterexample verifier: inverse reconstruction failed")
    return result


def systematic_right_columns() -> tuple[int, ...]:
    generator_rows = tuple((GENERATOR << row) | (1 << 127) for row in range(64))
    left = tuple(row & MASK64 for row in generator_rows)
    right = tuple((row >> 64) & MASK64 for row in generator_rows)
    inverse_left = inverse_columns(left)
    result = tuple(apply_columns(right, inverse_left[bit]) for bit in range(64))
    if any(apply_columns(left, inverse_left[bit]) != 1 << bit for bit in range(64)):
        raise RuntimeError("counterexample verifier: systematic map failed")
    return result


def accumulate_by_bits(value: int) -> int:
    running = 0
    result = 0
    for bit in range(64):
        running ^= (value >> bit) & 1
        result |= running << bit
    return result


def lifted_step(a: int, b: int, parity_columns: tuple[int, ...]) -> tuple[int, int, int]:
    first_output = accumulate_by_bits(a)
    retained_output = accumulate_by_bits(b ^ first_output)
    next_a = apply_columns(parity_columns, first_output)
    next_b = apply_columns(parity_columns, retained_output)
    return next_a, next_b, retained_output


def main() -> None:
    raw = json.loads(RAW.read_text())
    if raw["result"] != "COUNTEREXAMPLE":
        raise RuntimeError("counterexample verifier: raw result changed")
    counterexample = raw["counterexample"]
    a = int(counterexample["initial_a_hex"], 16)
    b = int(counterexample["initial_b_hex"], 16)
    if (a, b) != (EXPECTED_A, EXPECTED_B):
        raise RuntimeError("counterexample verifier: raw state changed")

    parity_columns = systematic_right_columns()
    initial_a = a
    initial_b = b
    outputs = []
    states = [(a, b)]
    for _ in range(18):
        a, b, output = lifted_step(a, b, parity_columns)
        outputs.append(output)
        states.append((a, b))
    weights = tuple(output.bit_count() for output in outputs)
    if weights != EXPECTED_WEIGHTS:
        raise RuntimeError(f"counterexample verifier: weights changed: {weights}")
    first_weight = sum(weights[:9])
    second_weight = sum(weights[9:])
    total_weight = first_weight + second_weight
    if (first_weight, second_weight, total_weight) != (47, 50, 97):
        raise RuntimeError("counterexample verifier: paired weights changed")
    if initial_a == 0 and initial_b == 0:
        raise RuntimeError("counterexample verifier: state is zero")
    if total_weight >= 106:
        raise RuntimeError("counterexample verifier: state does not violate the bound")

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal05-counterexample-verified-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "INDEPENDENT_EXACT_COUNTEREXAMPLE_REPLAY",
        "source_sha256": digest(Path(__file__).resolve()),
        "raw_counterexample_sha256": digest(RAW),
        "construction_constants": {
            "systematic_generator_hex": hex(GENERATOR),
            "state_bits": 128,
            "node_output_bits": 64,
        },
        "initial_state": {
            "a_hex": f"0x{initial_a:016x}",
            "b_hex": f"0x{initial_b:016x}",
        },
        "node_outputs_hex": [f"0x{output:016x}" for output in outputs],
        "node_weights": list(weights),
        "state_after_9_nodes": {
            "a_hex": f"0x{states[9][0]:016x}",
            "b_hex": f"0x{states[9][1]:016x}",
        },
        "first_9_weight": first_weight,
        "second_9_weight": second_weight,
        "paired_18_weight": total_weight,
        "claimed_lower_bound": 106,
        "strict_deficit": 106 - total_weight,
        "result": "PAIRED_RETURN_BOUND_REFUTED",
        "verification_scope": (
            "The verifier reconstructs the systematic parity map from the fixed "
            "generator polynomial and evaluates accumulation bit by bit. It does "
            "not import the search recurrence or its matrix tables."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"node_weights={list(weights)}")
    print(f"first_9_weight={first_weight}")
    print(f"second_9_weight={second_weight}")
    print(f"paired_18_weight={total_weight}")
    print(f"output={OUTPUT}")
    print("status=PAIRED_RETURN_BOUND_REFUTED")


if __name__ == "__main__":
    main()
