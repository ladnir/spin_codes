#!/usr/bin/env python3
"""Reduce the Goal 02 low spectrum to sparse dual enumerator counts."""

from __future__ import annotations

import hashlib
import json
import struct
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
OUTPUT = CANDIDATE / "receipts" / "goal02_sixth_moment_primary.json"
TAIL = CANDIDATE / "receipts" / "goal02_tail_list_primary.json"
COMPONENTS = CANDIDATE / "receipts" / "goal02_low_components_primary.json"
GENERATOR = 0xF4845518B9582A1F
MASK64 = (1 << 64) - 1
NODES = 32_772
LENGTH = 64 * NODES
THRESHOLD = 377_530
NONZERO_TARGET = 316_605


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def carryless_multiply_low(left: int, right: int) -> int:
    result = 0
    while right:
        bit = right & -right
        result ^= left << (bit.bit_length() - 1)
        right ^= bit
    return result & MASK64


def systematic_right_columns() -> tuple[int, ...]:
    inverse = 1
    for bit in range(1, 64):
        if (carryless_multiply_low(GENERATOR, inverse) >> bit) & 1:
            inverse |= 1 << bit
    if carryless_multiply_low(GENERATOR, inverse) != 1:
        raise RuntimeError("sixth moment: generator inverse failed")
    result = []
    for column in range(64):
        message = carryless_multiply_low(1 << column, inverse)
        low = 0
        high = 0
        while message:
            bit = message & -message
            shift = bit.bit_length() - 1
            low ^= (GENERATOR << shift) & MASK64
            high ^= (0 if shift == 0 else GENERATOR >> (64 - shift)) | (1 << 63)
            message ^= bit
        if low != 1 << column:
            raise RuntimeError("sixth moment: systematic reconstruction failed")
        result.append(high)
    return tuple(result)


def accumulate(value: int) -> int:
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
            value = 0
            for bit in range(8):
                if (byte >> bit) & 1:
                    value ^= columns[8 * byte_index + bit]
            table.append(value)
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


def transpose(columns: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(((columns[input_bit] >> output) & 1) << input_bit for input_bit in range(64))
        for output in range(64)
    )


def main() -> None:
    parity = systematic_right_columns()
    apply_parity = build_apply(parity)
    output_step_columns = tuple(accumulate(apply_parity(1 << bit)) for bit in range(64))
    transpose_step = build_apply(transpose(output_step_columns))
    current = [1 << bit for bit in range(64)]
    seen: set[int] = set()
    sequence_hash = hashlib.sha256()
    witness_values = {}
    witness_indices = (2, 126, 127, 129)
    coordinate_index = 0
    for _node in range(NODES):
        for value in current:
            if value == 0 or value in seen:
                raise RuntimeError("sixth moment: output coordinate is zero or repeated")
            seen.add(value)
            sequence_hash.update(struct.pack("<Q", value))
            if coordinate_index in witness_indices:
                witness_values[coordinate_index] = value
            coordinate_index += 1
        current = [transpose_step(value) for value in current]
    if len(seen) != LENGTH:
        raise RuntimeError("sixth moment: coordinate census is incomplete")
    if not all(index in witness_values for index in witness_indices):
        raise RuntimeError("sixth moment: local dual witness was not captured")
    witness_xor = 0
    for index in witness_indices:
        witness_xor ^= witness_values[index]
    if witness_xor != 0:
        raise RuntimeError("sixth moment: weight-four dual witness changed")

    deviation = LENGTH // 2 - THRESHOLD
    baseline = 15 * LENGTH**3 - 30 * LENGTH**2 + 16 * LENGTH
    coefficient4 = 360 * LENGTH - 960
    allowed_total = NONZERO_TARGET + 1
    scaled_sparse_budget = (
        allowed_total * 64 * deviation**6 - (1 << 64) * baseline
    )
    if scaled_sparse_budget <= 0:
        raise RuntimeError("sixth moment: sparse-relation budget is not positive")
    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-sixth-moment-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "EXACT_SPARSE_DUAL_REDUCTION",
        "source_sha256": digest(Path(__file__).resolve()),
        "tail_primary_sha256": digest(TAIL),
        "low_components_primary_sha256": digest(COMPONENTS),
        "response_length": LENGTH,
        "low_weight_threshold": THRESHOLD,
        "deviation_below_half": deviation,
        "coordinate_census": {
            "nonzero_coordinate_count": len(seen),
            "distinct_coordinate_count": len(seen),
            "ordered_coordinate_sequence_sha256": sequence_hash.hexdigest(),
        },
        "authenticated_weight4_dual_witness": {
            "flat_coordinate_indices": list(witness_indices),
            "node_coordinate_pairs": [divmod(index, 64) for index in witness_indices],
            "column_hex": [hex(witness_values[index]) for index in witness_indices],
            "xor_hex": hex(witness_xor),
        },
        "sixth_moment_identity": {
            "notation": "M6=E[(n-2X)^6] over uniform 64-bit states",
            "baseline_even_multiplicity_term": baseline,
            "weight4_dual_coefficient": coefficient4,
            "weight6_dual_coefficient": 720,
            "formula": "M6=baseline+(360*n-960)*A4_dual+720*A6_dual",
        },
        "markov_reduction": {
            "maximum_total_low_codewords": allowed_total,
            "maximum_nonzero_low_states": NONZERO_TARGET,
            "scaled_joint_sparse_budget": scaled_sparse_budget,
            "sufficient_inequality": (
                "2^64*((360*n-960)*A4_dual+720*A6_dual) "
                "<= scaled_joint_sparse_budget"
            ),
            "a4_cap_if_a6_zero": scaled_sparse_budget // ((1 << 64) * coefficient4),
            "a6_cap_if_a4_zero": scaled_sparse_budget // ((1 << 64) * 720),
            "baseline_total_low_codeword_bound_floor": (
                (1 << 64) * baseline // (64 * deviation**6)
            ),
        },
        "result": "EXACT_REDUCTION; SPARSE_DUAL_COUNTS_REQUIRED",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"output={OUTPUT}")
    print(f"distinct_coordinates={len(seen)}")
    print(f"weight4_witness={list(witness_indices)}")
    print(f"baseline_total_bound={payload['markov_reduction']['baseline_total_low_codeword_bound_floor']}")
    print("status=EXACT_SPARSE_DUAL_REDUCTION")


if __name__ == "__main__":
    main()
