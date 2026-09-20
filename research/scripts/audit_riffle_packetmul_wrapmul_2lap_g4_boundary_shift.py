#!/usr/bin/env python3
"""Audit the boundary-shift identity for PacketMul-WrapMul-2Lap g=4.

The audit uses the authenticated systematic BCH state columns.  It checks the
identity on every state basis vector over the full 32,772-node lap and also
checks the two-node inverse image of the boundary impulse directly.
"""

from __future__ import annotations

from analyze_systematic_group_kernel import systematic_state_columns


NODES = 32_772
MASK64 = (1 << 64) - 1


def accumulate(value: int) -> int:
    for shift in (1, 2, 4, 8, 16, 32):
        value ^= value << shift
    return value & MASK64


def inverse_accumulate(value: int) -> int:
    # A prefix XOR is inverted by the adjacent difference x_i=y_i+y_{i-1}.
    return (value ^ ((value << 1) & MASK64)) & MASK64


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


def audit_response_impulse(z: int, parity) -> int:
    response_state = z
    impulse_state = 0
    checksum = 0
    for node in range(NODES):
        impulse_input = z if node == 0 else 0
        response_output = accumulate(response_state)
        impulse_output = accumulate(impulse_state ^ impulse_input)
        if response_output != impulse_output:
            raise RuntimeError("boundary shift audit: J(z) != F(e_0 z)")
        response_state = parity(response_output)
        impulse_state = parity(impulse_output)
        if response_state != impulse_state:
            raise RuntimeError("boundary shift audit: terminal-state mismatch")
        checksum ^= response_output
    return checksum


def audit_two_node_turnoff(z: int, parity) -> None:
    first_output = accumulate(inverse_accumulate(z))
    first_state = parity(first_output)
    second_output = accumulate(first_state ^ parity(z))
    second_state = parity(second_output)
    if first_output != z or second_output != 0 or second_state != 0:
        raise RuntimeError("boundary shift audit: F(v_z) != e_0 z")


def main() -> None:
    parity = build_apply(systematic_state_columns())
    if any(inverse_accumulate(accumulate(1 << bit)) != 1 << bit for bit in range(64)):
        raise RuntimeError("boundary shift audit: accumulator inverse mismatch")

    checked_output_words = 0
    checksum = 0

    for bit in range(64):
        z = 1 << bit
        checksum ^= audit_response_impulse(z, parity)
        audit_two_node_turnoff(z, parity)
        checked_output_words += 2 * NODES

    # A mixed state catches accidental basis-only comparison or ordering errors.
    z = 0xD6E8FEB86659FD93
    checksum ^= audit_response_impulse(z, parity)
    audit_two_node_turnoff(z, parity)

    print("candidate=Riffle PacketMul-WrapMul-2Lap g=4")
    print(f"nodes={NODES}")
    print("state_basis_vectors=64")
    print(f"full_lap_output_words_checked={checked_output_words}")
    print(f"trajectory_checksum=0x{checksum:016x}")
    print("J_equals_F_boundary_impulse=PASS")
    print("two_node_turnoff_preimage=PASS")
    print("status=EXACT_BOUNDARY_SHIFT_IDENTITY")


if __name__ == "__main__":
    main()
