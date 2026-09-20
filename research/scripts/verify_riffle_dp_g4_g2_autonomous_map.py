#!/usr/bin/env python3
"""Independent replay of the exact Riffle DP g=4 autonomous-map receipt."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns


ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "explorations" / "riffle_dp_g4_g2_autonomous_map.json"
MASK64 = (1 << 64) - 1
EXPECTED_FACTORS = (0x3, 0x7, 0x13, 0x373, 0x519, 0x7C9C3, 0x1E1FFF)


def acc(value: int) -> int:
    for shift in (1, 2, 4, 8, 16, 32):
        value ^= value << shift
    return value & MASK64


def make_linear_map(columns: tuple[int, ...]):
    tables = []
    for byte_index in range(8):
        row = []
        for byte in range(256):
            result = 0
            for bit in range(8):
                if (byte >> bit) & 1:
                    result ^= columns[8 * byte_index + bit]
            row.append(result)
        tables.append(tuple(row))

    def apply(value: int) -> int:
        result = 0
        for byte_index, table in enumerate(tables):
            result ^= table[(value >> (8 * byte_index)) & 0xFF]
        return result

    return apply


def poly_mul(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        left <<= 1
        right >>= 1
    return result


def poly_divmod(dividend: int, divisor: int) -> tuple[int, int]:
    quotient = 0
    degree = divisor.bit_length() - 1
    while dividend and dividend.bit_length() - 1 >= degree:
        shift = dividend.bit_length() - 1 - degree
        quotient ^= 1 << shift
        dividend ^= divisor << shift
    return quotient, dividend


def poly_mod(value: int, modulus: int) -> int:
    return poly_divmod(value, modulus)[1]


def poly_gcd(left: int, right: int) -> int:
    while right:
        left, right = right, poly_mod(left, right)
    return left


def poly_mul_mod(left: int, right: int, modulus: int) -> int:
    return poly_mod(poly_mul(left, right), modulus)


def poly_pow_mod(value: int, exponent: int, modulus: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = poly_mul_mod(result, value, modulus)
        value = poly_mul_mod(value, value, modulus)
        exponent >>= 1
    return result


def integer_prime_factors(value: int) -> tuple[int, ...]:
    factors = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        factors.append(value)
    return tuple(factors)


def is_irreducible(polynomial: int) -> bool:
    degree = polynomial.bit_length() - 1
    x_value = poly_mod(2, polynomial)
    if poly_pow_mod(x_value, 1 << degree, polynomial) != x_value:
        return False
    return all(
        poly_gcd(
            poly_pow_mod(x_value, 1 << (degree // prime), polynomial) ^ x_value,
            polynomial,
        )
        == 1
        for prime in integer_prime_factors(degree)
    )


def root_order(polynomial: int) -> int:
    degree = polynomial.bit_length() - 1
    order = (1 << degree) - 1
    for prime in integer_prime_factors(order):
        while order % prime == 0 and poly_pow_mod(2, order // prime, polynomial) == 1:
            order //= prime
    return order


def apply_polynomial(step, polynomial: int, value: int) -> int:
    result = 0
    current = value
    while polynomial:
        if polynomial & 1:
            result ^= current
        current = step(current)
        polynomial >>= 1
    return result


def nullspace_basis(columns: tuple[int, ...]) -> tuple[int, ...]:
    pivots = [0] * 64
    witnesses = [0] * 64
    kernel = []
    for index, column in enumerate(columns):
        value = column
        witness = 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivots[pivot]:
                value ^= pivots[pivot]
                witness ^= witnesses[pivot]
            else:
                pivots[pivot] = value
                witnesses[pivot] = witness
                break
        if value == 0:
            kernel.append(witness)
    return tuple(kernel)


def replay_subspace(step, row: dict, factor_by_degree: dict[int, int]) -> None:
    degrees = (
        tuple(row["component_degrees"])
        if "component_degrees" in row
        else (row["degree"],)
    )
    factors = tuple(factor_by_degree[degree] for degree in degrees)
    annihilator = 1
    for factor in factors:
        annihilator = poly_mul(annihilator, factor)
    if hex(annihilator) != row.get("annihilator_hex", row.get("factor_hex")):
        raise SystemExit("G2 autonomous verifier: annihilator mismatch")
    dimension = sum(degrees)
    columns = tuple(
        apply_polynomial(step, annihilator, 1 << bit) for bit in range(64)
    )
    basis = nullspace_basis(columns)
    if len(basis) != dimension:
        raise SystemExit("G2 autonomous verifier: kernel dimension mismatch")
    states = [0]
    for vector in basis:
        states += [state ^ vector for state in states]
    complements = []
    for factor in factors:
        quotient, remainder = poly_divmod(annihilator, factor)
        if remainder:
            raise SystemExit("G2 autonomous verifier: factor division failed")
        complements.append(quotient)
    unseen = set(states[1:])
    averages = []
    full_averages = []
    minimum_weight = 65
    maximum_weight = 0
    while unseen:
        start = next(iter(unseen))
        full = all(apply_polynomial(step, complement, start) for complement in complements)
        current = start
        total = 0
        period = 0
        while True:
            unseen.remove(current)
            weight = acc(current).bit_count()
            minimum_weight = min(minimum_weight, weight)
            maximum_weight = max(maximum_weight, weight)
            total += weight
            period += 1
            current = step(current)
            if current == start:
                break
        average = Fraction(total, period)
        averages.append(average)
        if full:
            full_averages.append(average)
    if len(averages) != row["cycle_count"]:
        raise SystemExit("G2 autonomous verifier: cycle count mismatch")
    if minimum_weight != row["minimum_emitted_weight"] or maximum_weight != row["maximum_emitted_weight"]:
        raise SystemExit("G2 autonomous verifier: emitted-weight extrema mismatch")
    minimum = row["minimum_cycle_average"]
    maximum = row["maximum_cycle_average"]
    if min(averages) != Fraction(minimum["numerator"], minimum["denominator"]):
        raise SystemExit("G2 autonomous verifier: minimum cycle average mismatch")
    if max(averages) != Fraction(maximum["numerator"], maximum["denominator"]):
        raise SystemExit("G2 autonomous verifier: maximum cycle average mismatch")
    if full_averages:
        if len(full_averages) != row["full_component_support_cycle_count"]:
            raise SystemExit("G2 autonomous verifier: full-support cycle count mismatch")
        minimum = row["minimum_full_component_support_cycle_average"]
        maximum = row["maximum_full_component_support_cycle_average"]
        if min(full_averages) != Fraction(minimum["numerator"], minimum["denominator"]):
            raise SystemExit("G2 autonomous verifier: full-support minimum mismatch")
        if max(full_averages) != Fraction(maximum["numerator"], maximum["denominator"]):
            raise SystemExit("G2 autonomous verifier: full-support maximum mismatch")


def main() -> None:
    receipt = json.loads(RECEIPT.read_text())
    p_columns = systematic_state_columns()
    apply_p = make_linear_map(p_columns)

    def step(state: int) -> int:
        return apply_p(acc(state))

    columns = tuple(step(1 << bit) for bit in range(64))
    digest = hashlib.sha256(
        b"".join(column.to_bytes(8, "little") for column in columns)
    ).hexdigest()
    if digest != receipt["matrix_sha256"]:
        raise SystemExit("G2 autonomous verifier: matrix hash mismatch")
    factors = EXPECTED_FACTORS
    if not all(is_irreducible(factor) for factor in factors):
        raise SystemExit("G2 autonomous verifier: reducible committed factor")
    polynomial = 1
    for factor in factors:
        polynomial = poly_mul(polynomial, factor)
    if hex(polynomial) != receipt["minimal_polynomial_hex"]:
        raise SystemExit("G2 autonomous verifier: minimal polynomial mismatch")
    if any(apply_polynomial(step, polynomial, 1 << bit) for bit in range(64)):
        raise SystemExit("G2 autonomous verifier: polynomial does not annihilate map")
    factor_by_degree = {factor.bit_length() - 1: factor for factor in factors}
    component_rows = receipt["irreducible_components"]
    if [row["factor_hex"] for row in component_rows] != [hex(factor) for factor in factors]:
        raise SystemExit("G2 autonomous verifier: component list mismatch")
    orders = tuple(root_order(factor) for factor in factors)
    if orders != tuple(row["root_order"] for row in component_rows):
        raise SystemExit("G2 autonomous verifier: component order mismatch")
    if str(math.lcm(*orders)) != receipt["matrix_order"]:
        raise SystemExit("G2 autonomous verifier: matrix order mismatch")
    for row in component_rows:
        replay_subspace(step, row, factor_by_degree)
    for row in receipt["combined_subspaces_dimension_at_most_20"]:
        replay_subspace(step, row, factor_by_degree)
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"matrix_sha256={digest}")
    print("factor_degrees=" + ",".join(str(factor.bit_length() - 1) for factor in factors))
    print(f"matrix_order={math.lcm(*orders)}")
    print("status=EXACT_G2_AUTONOMOUS_MAP_VERIFIED")


if __name__ == "__main__":
    main()
