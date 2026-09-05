#!/usr/bin/env python3
"""Verify the finite-field moduli and one-vector ranks used by the inner maps."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path


MODULUS_LOW = {
    31: 0x0F,
    32: 0x1000B,
    46: 0x807,
    47: 0x8B,
    48: 0x20007,
    64: 0x807,
}
DEFAULT_OUTPUT = Path(
    "constructions/riffle_bchperm_transpose_bitshuffle_structuredstepconv/"
    "receipts/algebraic_family_checks.json"
)


def polynomial_mod(value: int, modulus: int) -> int:
    degree = modulus.bit_length() - 1
    while value.bit_length() - 1 >= degree:
        value ^= modulus << (value.bit_length() - 1 - degree)
    return value


def polynomial_gcd(left: int, right: int) -> int:
    while right:
        left, right = right, polynomial_mod(left, right)
    return left


def field_multiply(left: int, right: int, modulus: int) -> int:
    product = 0
    while right:
        if right & 1:
            product ^= left
        right >>= 1
        left <<= 1
    return polynomial_mod(product, modulus)


def prime_factors(value: int) -> list[int]:
    factors = []
    candidate = 2
    while candidate * candidate <= value:
        if value % candidate == 0:
            factors.append(candidate)
            while value % candidate == 0:
                value //= candidate
        candidate += 1
    if value > 1:
        factors.append(value)
    return factors


def is_irreducible(modulus: int, degree: int) -> bool:
    x = 2
    power = x
    for _ in range(degree):
        power = field_multiply(power, power, modulus)
    if power != x:
        return False
    for divisor in prime_factors(degree):
        power = x
        for _ in range(degree // divisor):
            power = field_multiply(power, power, modulus)
        if polynomial_gcd(power ^ x, modulus) != 1:
            return False
    return True


def binary_rank(columns: list[int], rows: int) -> int:
    basis = [0] * rows
    rank = 0
    for value in columns:
        while value:
            pivot = value.bit_length() - 1
            if basis[pivot]:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                rank += 1
                break
    return rank


def toeplitz_image(value: int, seed: int, degree: int) -> int:
    output = 0
    for column in range(degree):
        parity = 0
        for row in range(degree):
            parity ^= ((value >> row) & 1) & (
                (seed >> (degree - 1 + column - row)) & 1
            )
        output |= parity << column
    return output


def check_degree(degree: int, modulus_low: int) -> dict[str, object]:
    modulus = (1 << degree) | modulus_low
    if not is_irreducible(modulus, degree):
        raise AssertionError(f"degree-{degree} modulus is reducible")
    generator = random.Random(0x72696666 ^ degree)
    vectors = [1 << bit for bit in range(degree)]
    vectors.extend(generator.randrange(1, 1 << degree) for _ in range(16))
    minimum_field_rank = degree
    minimum_toeplitz_rank = degree
    for value in vectors:
        field_columns = [
            field_multiply(value, 1 << bit, modulus) for bit in range(degree)
        ]
        minimum_field_rank = min(
            minimum_field_rank, binary_rank(field_columns, degree)
        )
        toeplitz_columns = [
            toeplitz_image(value, 1 << bit, degree)
            for bit in range(2 * degree - 1)
        ]
        minimum_toeplitz_rank = min(
            minimum_toeplitz_rank, binary_rank(toeplitz_columns, degree)
        )
    if min(minimum_field_rank, minimum_toeplitz_rank) != degree:
        raise AssertionError(f"degree-{degree} one-vector map lost rank")
    return {
        "degree": degree,
        "modulus_hex": hex(modulus),
        "irreducible": True,
        "tested_nonzero_vectors": len(vectors),
        "minimum_field_coefficient_map_rank": minimum_field_rank,
        "minimum_toeplitz_seed_map_rank": minimum_toeplitz_rank,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    rows = [check_degree(degree, low) for degree, low in MODULUS_LOW.items()]
    payload = {
        "schema": "riffle-structured-step-algebraic-checks-v1",
        "rows": rows,
        "scope": (
            "Rabin irreducibility checks are exact. Rank checks cover all basis "
            "vectors and sixteen deterministic random nonzero vectors per degree."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
