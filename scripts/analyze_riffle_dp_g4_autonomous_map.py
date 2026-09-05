#!/usr/bin/env python3
"""Exact algebraic audit of the Riffle DP g=4 zero-input state map."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns


REPOSITORY = Path(__file__).resolve().parents[1]
OUTPUT = EXPLORATIONS = (
    REPOSITORY / "explorations" / "riffle_dp_g4_g2_autonomous_map.json"
)
MASK64 = (1 << 64) - 1
ORDER_FACTORS = (3, 5, 17, 257, 641, 65537, 6700417)


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


def polynomial_mod(value: int, modulus: int) -> int:
    degree = modulus.bit_length() - 1
    while value.bit_length() - 1 >= degree:
        value ^= modulus << (value.bit_length() - 1 - degree)
    return value


def polynomial_multiply_mod(left: int, right: int, modulus: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left.bit_length() - 1 == modulus.bit_length() - 1:
            left ^= modulus
    return result


def polynomial_power_mod(value: int, exponent: int, modulus: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = polynomial_multiply_mod(result, value, modulus)
        value = polynomial_multiply_mod(value, value, modulus)
        exponent >>= 1
    return result


def polynomial_gcd(left: int, right: int) -> int:
    while right:
        left, right = right, polynomial_mod(left, right)
    return left


def polynomial_multiply(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        left <<= 1
        right >>= 1
    return result


def polynomial_divide_exact(dividend: int, divisor: int) -> int:
    quotient = 0
    divisor_degree = divisor.bit_length() - 1
    while dividend:
        shift = dividend.bit_length() - 1 - divisor_degree
        if shift < 0:
            break
        quotient ^= 1 << shift
        dividend ^= divisor << shift
    if dividend:
        raise RuntimeError("polynomial division left a remainder")
    return quotient


def polynomial_lcm(left: int, right: int) -> int:
    common = polynomial_gcd(left, right)
    return polynomial_multiply(polynomial_divide_exact(left, common), right)


def binary_rank(columns: tuple[int, ...]) -> int:
    pivots = [0] * 64
    rank = 0
    for column in columns:
        value = column
        while value:
            pivot = value.bit_length() - 1
            if pivots[pivot]:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                rank += 1
                break
    return rank


def apply_polynomial(step, polynomial: int, value: int) -> int:
    result = 0
    current = value
    for exponent in range(polynomial.bit_length()):
        if (polynomial >> exponent) & 1:
            result ^= current
        current = step(current)
    return result


def kernel_basis(columns: tuple[int, ...]) -> tuple[int, ...]:
    pivot_values = [0] * 64
    pivot_representations = [0] * 64
    result = []
    for column_index, column in enumerate(columns):
        value = column
        representation = 1 << column_index
        while value:
            pivot = value.bit_length() - 1
            if pivot_values[pivot]:
                value ^= pivot_values[pivot]
                representation ^= pivot_representations[pivot]
            else:
                pivot_values[pivot] = value
                pivot_representations[pivot] = representation
                break
        if value == 0:
            result.append(representation)
    return tuple(result)


def distinct_degree_factors(polynomial: int) -> tuple[tuple[int, int], ...]:
    remainder = polynomial
    x_value = 2
    x_power = x_value
    degree = 0
    result = []
    while remainder.bit_length() - 1 >= 2:
        degree += 1
        x_power = polynomial_multiply_mod(x_power, x_power, remainder)
        factor = polynomial_gcd(x_power ^ x_value, remainder)
        if factor != 1:
            result.append((degree, factor))
            remainder = polynomial_divide_exact(remainder, factor)
            if remainder == 1:
                break
            x_power = polynomial_mod(x_power, remainder)
        if 2 * degree > remainder.bit_length() - 1:
            break
    if remainder != 1:
        result.append((remainder.bit_length() - 1, remainder))
    if polynomial_multiply_all(factor for _, factor in result) != polynomial:
        raise RuntimeError("factorization does not reconstruct the minimal polynomial")
    return tuple(result)


def polynomial_multiply_all(values) -> int:
    result = 1
    for value in values:
        result = polynomial_multiply(result, value)
    return result


def prime_factors(value: int) -> tuple[int, ...]:
    result = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            result.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        result.append(value)
    return tuple(result)


def root_order(factor: int, degree: int) -> int:
    order = (1 << degree) - 1
    for prime in prime_factors(order):
        while order % prime == 0 and polynomial_power_mod(2, order // prime, factor) == 1:
            order //= prime
    return order


def subspace_cycles(
    step,
    annihilator: int,
    dimension: int,
    orbit_order: int,
    component_factors: tuple[int, ...] = (),
) -> dict:
    columns = tuple(
        apply_polynomial(step, annihilator, 1 << bit) for bit in range(64)
    )
    basis = kernel_basis(columns)
    if len(basis) != dimension:
        raise RuntimeError("invariant subspace has the wrong kernel dimension")
    states = [0]
    for basis_value in basis:
        states += [state ^ basis_value for state in states]
    unseen = set(states[1:])
    cycle_rows = []
    complementary_polynomials = tuple(
        polynomial_divide_exact(annihilator, factor)
        for factor in component_factors
    )
    while unseen:
        start = next(iter(unseen))
        full_component_support = all(
            apply_polynomial(step, complement, start) != 0
            for complement in complementary_polynomials
        )
        current = start
        total_weight = 0
        minimum_weight = 65
        maximum_weight = 0
        period = 0
        while True:
            unseen.remove(current)
            emitted_weight = accumulate(current).bit_count()
            total_weight += emitted_weight
            minimum_weight = min(minimum_weight, emitted_weight)
            maximum_weight = max(maximum_weight, emitted_weight)
            period += 1
            current = step(current)
            if current == start:
                break
        cycle_rows.append(
            {
                "period": period,
                "total_emitted_weight": total_weight,
                "average_emitted_weight": Fraction(total_weight, period),
                "minimum_emitted_weight": minimum_weight,
                "maximum_emitted_weight": maximum_weight,
                "full_component_support": full_component_support,
            }
        )
    minimum_average = min(row["average_emitted_weight"] for row in cycle_rows)
    maximum_average = max(row["average_emitted_weight"] for row in cycle_rows)
    result = {
        "annihilator_hex": hex(annihilator),
        "dimension": dimension,
        "orbit_order": orbit_order,
        "cycle_count": len(cycle_rows),
        "minimum_cycle_average": {
            "numerator": minimum_average.numerator,
            "denominator": minimum_average.denominator,
            "decimal": float(minimum_average),
        },
        "maximum_cycle_average": {
            "numerator": maximum_average.numerator,
            "denominator": maximum_average.denominator,
            "decimal": float(maximum_average),
        },
        "minimum_emitted_weight": min(
            row["minimum_emitted_weight"] for row in cycle_rows
        ),
        "maximum_emitted_weight": max(
            row["maximum_emitted_weight"] for row in cycle_rows
        ),
    }
    full_support_rows = [
        row for row in cycle_rows if row["full_component_support"]
    ]
    if full_support_rows:
        minimum_full_average = min(
            row["average_emitted_weight"] for row in full_support_rows
        )
        maximum_full_average = max(
            row["average_emitted_weight"] for row in full_support_rows
        )
        result["full_component_support_cycle_count"] = len(full_support_rows)
        result["minimum_full_component_support_cycle_average"] = {
            "numerator": minimum_full_average.numerator,
            "denominator": minimum_full_average.denominator,
            "decimal": float(minimum_full_average),
        }
        result["maximum_full_component_support_cycle_average"] = {
            "numerator": maximum_full_average.numerator,
            "denominator": maximum_full_average.denominator,
            "decimal": float(maximum_full_average),
        }
    return result


def component_cycles(step, factor: int, degree: int) -> dict:
    order = root_order(factor, degree)
    result = subspace_cycles(step, factor, degree, order, (factor,))
    result["factor_hex"] = result.pop("annihilator_hex")
    result["degree"] = result.pop("dimension")
    result["root_order"] = result.pop("orbit_order")
    return result


def krylov_polynomial(step, start: int) -> int:
    basis_values = [0] * 64
    basis_representations = [0] * 64
    value = start
    for exponent in range(65):
        reduced = value
        representation = 1 << exponent
        while reduced:
            pivot = reduced.bit_length() - 1
            if basis_values[pivot]:
                reduced ^= basis_values[pivot]
                representation ^= basis_representations[pivot]
            else:
                basis_values[pivot] = reduced
                basis_representations[pivot] = representation
                break
        if reduced == 0:
            return representation
        value = step(value)
    raise RuntimeError("Krylov sequence did not produce a degree-64 relation")


def annihilates(step, polynomial: int, value: int) -> bool:
    result = 0
    current = value
    for exponent in range(polynomial.bit_length()):
        if (polynomial >> exponent) & 1:
            result ^= current
        current = step(current)
    return result == 0


def main() -> None:
    state_columns = systematic_state_columns()
    apply_state = build_apply(state_columns)

    def autonomous(state: int) -> int:
        return apply_state(accumulate(state))

    columns = tuple(autonomous(1 << bit) for bit in range(64))
    if len(set(columns)) != 64 or 0 in columns:
        raise RuntimeError("autonomous map has a repeated or zero column")
    basis_polynomials = tuple(krylov_polynomial(autonomous, 1 << bit) for bit in range(64))
    polynomial = 1
    for basis_polynomial in basis_polynomials:
        polynomial = polynomial_lcm(polynomial, basis_polynomial)
    degree = polynomial.bit_length() - 1
    if any(not annihilates(autonomous, polynomial, 1 << bit) for bit in range(64)):
        raise RuntimeError("minimal polynomial does not annihilate the matrix")

    x_value = 2
    irreducible = degree == 64 and (
        polynomial_power_mod(x_value, 1 << degree, polynomial) == x_value
        and polynomial_gcd(
            polynomial_power_mod(x_value, 1 << (degree // 2), polynomial) ^ x_value,
            polynomial,
        )
        == 1
    )
    order = (1 << 64) - 1
    primitive = irreducible and all(
        polynomial_power_mod(x_value, order // factor, polynomial) != 1
        for factor in ORDER_FACTORS
    )
    rank = binary_rank(columns)
    factors = distinct_degree_factors(polynomial)
    factor_rows = [
        component_cycles(autonomous, factor, factor_degree)
        for factor_degree, factor in factors
    ]
    matrix_order = math.lcm(*(row["root_order"] for row in factor_rows))
    factor_by_degree = {
        factor_degree: factor for factor_degree, factor in factors
    }
    order_by_degree = {
        row["degree"]: row["root_order"] for row in factor_rows
    }
    combined_degree_sets = (
        (20,),
        (18, 2),
        (18, 1),
        (10, 9, 1),
        (10, 4, 2, 1),
        (9, 4, 2, 1),
    )
    combined_rows = []
    for degree_set in combined_degree_sets:
        annihilator = polynomial_multiply_all(
            factor_by_degree[component_degree]
            for component_degree in degree_set
        )
        combined = subspace_cycles(
            autonomous,
            annihilator,
            sum(degree_set),
            math.lcm(*(order_by_degree[component_degree] for component_degree in degree_set)),
            tuple(factor_by_degree[component_degree] for component_degree in degree_set),
        )
        combined["component_degrees"] = list(degree_set)
        combined_rows.append(combined)

    digest = hashlib.sha256(
        b"".join(column.to_bytes(8, "little") for column in columns)
    ).hexdigest()
    payload = {
        "schema": "riffle-dp-g4-g2-autonomous-map-v1",
        "candidate": "Riffle DP g=4@g0-v1",
        "evidence_label": "EXACT",
        "matrix_sha256": digest,
        "minimal_polynomial_hex": hex(polynomial),
        "degree": degree,
        "matrix_rank": rank,
        "basis_krylov_degree_histogram": {
            str(basis_degree): sum(
                basis_polynomial.bit_length() - 1 == basis_degree
                for basis_polynomial in basis_polynomials
            )
            for basis_degree in sorted(
                {basis_polynomial.bit_length() - 1 for basis_polynomial in basis_polynomials}
            )
        },
        "irreducible": irreducible,
        "primitive": primitive,
        "irreducible_components": factor_rows,
        "combined_subspaces_dimension_at_most_20": combined_rows,
        "matrix_order": str(matrix_order),
        "interpretation": (
            "The map is invertible but decomposes into the listed irreducible components. "
            "The combined rows exhaust selected maximal component subsets of dimension at most 20."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("candidate=Riffle DP g=4@g0-v1")
    print(f"minimal_polynomial_hex={hex(polynomial)}")
    print(f"minimal_polynomial_degree={degree}")
    print(f"matrix_rank={rank}")
    print(f"irreducible={str(irreducible).lower()}")
    print(f"primitive={str(primitive).lower()}")
    print(
        "factor_degrees="
        + ",".join(str(row["degree"]) for row in factor_rows)
    )
    print(f"matrix_order={matrix_order}")
    print(
        "component_minimum_cycle_averages="
        + ",".join(f"{row['minimum_cycle_average']['decimal']:.9f}" for row in factor_rows)
    )
    print(
        "combined_minimum_cycle_averages="
        + ",".join(f"{row['minimum_cycle_average']['decimal']:.9f}" for row in combined_rows)
    )
    print(f"output={OUTPUT}")
    print("status=EXACT_G2_AUTONOMOUS_MAP")


if __name__ == "__main__":
    main()
