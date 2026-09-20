#!/usr/bin/env python3
"""Independent replay of the PacketMul Goal 01 algebra and witnesses."""

from __future__ import annotations

import collections
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_2lap_g4"
MANIFEST = CANDIDATE / "manifest.json"
PRIMARY = CANDIDATE / "receipts" / "goal01_zero_symbol_primary.json"
SEARCH = CANDIDATE / "receipts" / "goal01_mixed_character_search.json"
OUTPUT = CANDIDATE / "receipts" / "goal01_zero_symbol_independent.json"
GENERATOR = 0xF4845518B9582A1F
MASK64 = (1 << 64) - 1
FIELD_MODULUS = 0x13
INNER_NODES = 32_772
PACKET_SLOTS = 16
SAMPLE_SPACE = INNER_NODES * PACKET_SLOTS


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def apply_sparse(columns: tuple[int, ...], value: int) -> int:
    result = 0
    while value:
        bit = value & -value
        result ^= columns[bit.bit_length() - 1]
        value ^= bit
    return result


def invert_columns(columns: tuple[int, ...]) -> tuple[int, ...]:
    """Invert a 64-by-64 column matrix through an independent row reduction."""
    augmented = []
    for row_index in range(64):
        left_row = sum(
            ((columns[column_index] >> row_index) & 1) << column_index
            for column_index in range(64)
        )
        augmented.append(left_row | (1 << (64 + row_index)))
    for pivot_index in range(64):
        pivot_row = next(
            row_index
            for row_index in range(pivot_index, 64)
            if (augmented[row_index] >> pivot_index) & 1
        )
        augmented[pivot_index], augmented[pivot_row] = (
            augmented[pivot_row],
            augmented[pivot_index],
        )
        for row_index in range(64):
            if row_index != pivot_index and (
                (augmented[row_index] >> pivot_index) & 1
            ):
                augmented[row_index] ^= augmented[pivot_index]
    inverse_rows = tuple(row >> 64 for row in augmented)
    return tuple(
        sum(
            ((inverse_rows[input_index] >> output_index) & 1) << input_index
            for input_index in range(64)
        )
        for output_index in range(64)
    )


def reconstruct_state_columns() -> tuple[int, ...]:
    generator_rows = tuple((GENERATOR << row) | (1 << 127) for row in range(64))
    left_columns = tuple(row & MASK64 for row in generator_rows)
    right_columns = tuple((row >> 64) & MASK64 for row in generator_rows)
    inverse_left = invert_columns(left_columns)
    columns = tuple(
        apply_sparse(right_columns, inverse_left[column]) for column in range(64)
    )
    if any(
        apply_sparse(left_columns, inverse_left[column]) != 1 << column
        for column in range(64)
    ):
        raise RuntimeError("independent systematic reconstruction failed")
    return columns


def byte_tables(columns: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    result = []
    for byte_index in range(8):
        block = columns[8 * byte_index : 8 * byte_index + 8]
        result.append(
            tuple(apply_sparse(block, value) for value in range(256))
        )
    return tuple(result)


def apply_tables(tables: tuple[tuple[int, ...], ...], value: int) -> int:
    result = 0
    for byte_index in range(8):
        result ^= tables[byte_index][(value >> (8 * byte_index)) & 255]
    return result


def accumulate_by_bits(value: int) -> int:
    running = 0
    result = 0
    for bit_index in range(64):
        running ^= (value >> bit_index) & 1
        result |= running << bit_index
    return result


def reconstruct_basis_contributions() -> np.ndarray:
    state_tables = byte_tables(reconstruct_state_columns())

    def step(state: int) -> int:
        return apply_tables(state_tables, accumulate_by_bits(state))

    result = np.empty((4, SAMPLE_SPACE), dtype=np.uint64)
    for packet_bit in range(4):
        cursor = 0
        for slot in range(PACKET_SLOTS):
            state = (1 << packet_bit) << (4 * slot)
            for _ in range(INNER_NODES):
                state = step(state)
                result[packet_bit, cursor] = state
                cursor += 1
    return result


def value_states(basis_states: np.ndarray, value: int) -> np.ndarray:
    result = np.zeros(SAMPLE_SPACE, dtype=np.uint64)
    for bit_index in range(4):
        if (value >> bit_index) & 1:
            result ^= basis_states[bit_index]
    return result


def gf16_multiply_schoolbook(left: int, right: int) -> int:
    product = 0
    for left_bit in range(4):
        for right_bit in range(4):
            if ((left >> left_bit) & 1) and ((right >> right_bit) & 1):
                product ^= 1 << (left_bit + right_bit)
    for degree in range(6, 3, -1):
        if (product >> degree) & 1:
            product ^= FIELD_MODULUS << (degree - 4)
    return product


def gf16_trace(value: int) -> int:
    conjugates = [value]
    for _ in range(3):
        conjugates.append(gf16_multiply_schoolbook(conjugates[-1], conjugates[-1]))
    result = conjugates[0] ^ conjugates[1] ^ conjugates[2] ^ conjugates[3]
    if result not in (0, 1):
        raise RuntimeError("independent trace escaped GF(2)")
    return result


def independent_dual_table() -> tuple[int, ...]:
    result = [-1] * 16
    for coefficient in range(16):
        functional = 0
        for input_bit in range(4):
            functional |= gf16_trace(
                gf16_multiply_schoolbook(coefficient, 1 << input_bit)
            ) << input_bit
        if result[functional] != -1:
            raise RuntimeError("independent trace pairing is degenerate")
        result[functional] = coefficient
    return tuple(result)


def coefficient_histogram(basis_states: np.ndarray, character: int) -> list[int]:
    raw = np.zeros(SAMPLE_SPACE, dtype=np.uint8)
    for input_bit in range(4):
        parity = np.fromiter(
            (
                (character & int(state)).bit_count() & 1
                for state in basis_states[input_bit]
            ),
            dtype=np.uint8,
            count=SAMPLE_SPACE,
        )
        raw |= parity << input_bit
    dual = np.asarray(independent_dual_table(), dtype=np.uint8)
    return np.bincount(dual[raw], minlength=16).astype(int).tolist()


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    if primary["active_manifest_sha256"] != digest(MANIFEST):
        raise RuntimeError("primary receipt does not bind the active manifest")

    basis_states = reconstruct_basis_contributions()
    basis_digest = hashlib.sha256(
        basis_states.astype("<u8", copy=False).tobytes()
    ).hexdigest()
    if basis_digest != primary["basis_contribution_sha256_little_endian_u64"]:
        raise RuntimeError("independent basis contribution map differs")
    dual = independent_dual_table()
    if list(dual) != primary["trace_dual_table_raw_functional_to_field_coefficient"]:
        raise RuntimeError("independent trace-dual table differs")

    values = []
    value_digests = {}
    for value in range(1, 16):
        states = value_states(basis_states, value)
        values.append(states)
        value_digests[str(value)] = hashlib.sha256(
            states.astype("<u8", copy=False).tobytes()
        ).hexdigest()
    if value_digests != primary["value_contribution_sha256_little_endian_u64"]:
        raise RuntimeError("independent packet-value maps differ")

    all_states = np.concatenate(values)
    _unique, multiplicities = np.unique(all_states, return_counts=True)
    multiplicity_histogram = collections.Counter(map(int, multiplicities))
    collision_sum = sum(
        multiplicity * multiplicity * frequency
        for multiplicity, frequency in multiplicity_histogram.items()
    )
    primary_collisions = primary["cross_value_collisions"]
    if {
        str(key): value for key, value in sorted(multiplicity_histogram.items())
    } != primary_collisions["multiplicity_histogram"]:
        raise RuntimeError("independent collision histogram differs")
    if str(collision_sum) != primary_collisions["ordered_collision_count"]:
        raise RuntimeError("independent collision count differs")

    expected_replays: dict[str, dict] = {
        row["character_hex"]: row for row in primary["primary_maximizer_replays"]
    }
    requested_characters = set(expected_replays)
    for component in primary["pure_component_rows"]:
        requested_characters.add(component["maximum_zero_symbol"]["character_hex"])
        requested_characters.add(component["maximum_absolute_beta"]["character_hex"])
    replay_rows = []
    for character_hex in sorted(requested_characters, key=lambda value: int(value, 16)):
        character = int(character_hex, 16)
        histogram = coefficient_histogram(basis_states, character)
        expected = expected_replays[character_hex]
        if histogram != expected["trace_coefficient_histogram_0_through_15"]:
            raise RuntimeError(f"independent coefficient replay differs for {character_hex}")
        zero_count = histogram[0]
        if zero_count != expected["zero_symbol_count"]:
            raise RuntimeError(f"independent zero count differs for {character_hex}")
        replay_rows.append(
            {
                "character_hex": character_hex,
                "trace_coefficient_histogram_0_through_15": histogram,
                "zero_symbol_count": zero_count,
                "signed_beta_numerator": 16 * zero_count - SAMPLE_SPACE,
                "beta_denominator": 15 * SAMPLE_SPACE,
            }
        )

    mixed_replay = None
    search_sha256 = None
    if SEARCH.exists():
        search = json.loads(SEARCH.read_text(encoding="utf-8"))
        search_sha256 = digest(SEARCH)
        reported = search["best_exact_character_found"]
        character = int(reported["character_hex"], 16)
        histogram = coefficient_histogram(basis_states, character)
        zero_count = histogram[0]
        signed_numerator = 16 * zero_count - SAMPLE_SPACE
        if zero_count != reported["zero_symbol_count"]:
            raise RuntimeError("independent mixed-witness zero count differs")
        if signed_numerator != reported["signed_beta_numerator"]:
            raise RuntimeError("independent mixed-witness beta differs")
        mixed_replay = {
            "character_hex": hex(character),
            "trace_coefficient_histogram_0_through_15": histogram,
            "zero_symbol_count": zero_count,
            "signed_beta_numerator": signed_numerator,
            "absolute_beta_numerator": abs(signed_numerator),
            "beta_denominator": 15 * SAMPLE_SPACE,
        }

    payload = {
        "schema": "riffle-packetmul-2lap-g4-goal01-zero-symbol-independent-v1",
        "candidate": "Riffle PacketMul-2Lap g=4",
        "candidate_id": manifest["candidate_id"],
        "evidence_label": "EXACT_INDEPENDENT_REPLAY",
        "active_manifest_sha256": digest(MANIFEST),
        "primary_receipt_sha256": digest(PRIMARY),
        "mixed_search_receipt_sha256": search_sha256,
        "independence_boundary": (
            "This verifier imports no construction-map helper. It reconstructs the "
            "systematic BCH half by separate row reduction, evaluates the accumulator "
            "bit by bit, implements GF(16) multiplication by schoolbook reduction, and "
            "recomputes the cross-value collision histogram."
        ),
        "basis_contribution_sha256_little_endian_u64": basis_digest,
        "trace_dual_table_raw_functional_to_field_coefficient": list(dual),
        "value_contribution_sha256_little_endian_u64": value_digests,
        "collision_histogram": {
            str(key): value for key, value in sorted(multiplicity_histogram.items())
        },
        "ordered_collision_count": str(collision_sum),
        "maximizer_replays": replay_rows,
        "best_mixed_witness_replay": mixed_replay,
        "all_primary_checks_passed": True,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"basis_contribution_sha256={basis_digest}")
    print(f"ordered_collision_count={collision_sum}")
    print(f"maximizers_replayed={len(replay_rows)}")
    print(f"mixed_witness_replayed={mixed_replay is not None}")
    print(f"output={OUTPUT}")
    print("status=EXACT_INDEPENDENT_ZERO_SYMBOL_REPLAY")


if __name__ == "__main__":
    main()
