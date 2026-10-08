"""Exact outer-parity and state-translation cost references.

This file verifies algebra and operation counts only.  It performs no
native benchmark and supplies no full-code distance certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb
from random import Random

from grs32_24_parity_reference import Cost, inverse, multiply, parity_transpose
from region64_reference import phase_geometry, matmul, matvec, transpose
from rate_design_coset import power, xor_sum


EXTRA_ROW = (225, 153, 25, 168, 71, 121, 240, 55, 199, 48, 100, 145, 128, 251, 212, 207, 151, 144, 69, 220, 104, 24, 113, 204)
SCALAR_BASIS = (1, 225, 153, 25, 168, 71, 240, 100)
SCALAR_XORS = (
    (105, 153, 240),
    (112, 25, 105),
    (204, 168, 100),
    (224, 1, 225),
    (55, 71, 112),
    (128, 153, 25),
    (49, 153, 168),
    (24, 1, 25),
    (166, 225, 71),
    (44, 204, 224),
    (121, 224, 153),
    (199, 128, 71),
    (48, 49, 1),
    (145, 112, 225),
    (251, 55, 204),
    (212, 24, 204),
    (207, 166, 105),
    (151, 166, 49),
    (144, 224, 112),
    (69, 44, 105),
    (220, 44, 240),
    (104, 105, 1),
    (113, 112, 1),
)


def interpolation_row(point, dimension=24):
    row = []
    for source in range(dimension):
        numerator = denominator = 1
        for other in range(dimension):
            if other != source:
                numerator = multiply(numerator, point ^ other)
                denominator = multiply(denominator, source ^ other)
        row.append(multiply(numerator, inverse(denominator)))
    return tuple(row)


def extra_transpose(value, offset=0, cost=None):
    """Point 32+offset, offset 0..7: seven products and 23 binary XORs."""
    if not 0 <= offset < 8:
        raise ValueError("extra point must lie in 32+{0,..,7}")
    products = {1: value}
    for coefficient in SCALAR_BASIS[1:]:
        products[coefficient] = multiply(coefficient, value)
        if cost is not None:
            cost.products += 1
    for target, left, right in SCALAR_XORS:
        products[target] = products[left] ^ products[right]
        if cost is not None:
            cost.xor2 += 1
    return tuple(products[EXTRA_ROW[column ^ offset]] for column in range(24))


def extended_parity_transpose(values, cost=None):
    if not 8 <= len(values) <= 16:
        raise ValueError("between 8 and 16 parity inputs required")
    result = list(parity_transpose(values[:8], cost))
    for offset, value in enumerate(values[8:]):
        extra = extra_transpose(value, offset, cost)
        result = [a ^ b for a, b in zip(result, extra)]
        if cost is not None:
            cost.xor2 += 24
    return tuple(result)


def shortened_outer_transpose(values, cost=None):
    """Complete GRS[32,23] transpose using 35 products and 139 booleans."""
    if len(values) != 32:
        raise ValueError("32 encoded symbols required")
    parity = parity_transpose(values[24:], cost)
    missing = values[23] ^ parity[23]
    low = multiply(12, missing)
    middle = low ^ missing
    if cost is not None:
        cost.products += 1
        cost.xor2 += 2
        cost.xor3 += 23
    return tuple(values[j] ^ parity[j] ^ (low if j < 8 else middle if j < 16 else missing)
                 for j in range(23))


def two_missing_correction(t22, t23, cost=None):
    """Seven products and 32 Boolean instructions for both missing values."""
    def xor2(a, b):
        if cost is not None:
            cost.xor2 += 1
        return a ^ b

    def xor3(a, b, c):
        if cost is not None:
            cost.xor3 += 1
        return a ^ b ^ c

    def cube(base, directions, count):
        result = [base]
        for index in range(1, count):
            previous = index & (index-1)
            bit = (index ^ previous).bit_length()-1
            result.append(xor2(result[previous], directions[bit]))
        return result

    s = xor2(t22, t23)
    s2, s4, s8, s16, s32, s64 = (multiply(coefficient, s) for coefficient in (2, 4, 8, 16, 32, 64))
    w = xor2(t22, xor3(s2, s4, s16))
    low = multiply(12, w)
    if cost is not None:
        cost.products += 7
    middle = xor2(low, w)
    s12, s24, s48 = xor2(s4, s8), xor2(s8, s16), xor2(s16, s32)
    s13, s26, s52 = xor2(s12, s), xor2(s24, s2), xor2(s48, s4)
    s104 = xor3(s8, s32, s64)
    return tuple(cube(low, (s12, s24, s48), 8)
                 + cube(xor2(middle, s104), (s13, s26, s52), 8)
                 + cube(xor2(w, s16), (s, s2, s4), 6))


def twice_shortened_outer_transpose(values, cost=None):
    """Complete GRS[32,22] transpose: 41 products and 170 booleans."""
    if len(values) != 32:
        raise ValueError("32 encoded symbols required")
    parity = parity_transpose(values[24:], cost)
    correction = two_missing_correction(values[22] ^ parity[22], values[23] ^ parity[23], cost)
    if cost is not None:
        cost.xor2 += 2
        cost.xor3 += 22
    return tuple(values[j] ^ parity[j] ^ correction[j] for j in range(22))


def translation_matrix(state, delta):
    return tuple(tuple(power(delta, i-j) if j <= i and comb(i, j) & 1 else 0
                       for j in range(state)) for i in range(state))


def translate(values, delta, adjoint=False):
    """Truncated Boolean-cube shears; return state and state-product count."""
    result = list(values)
    bits = [1 << bit for bit in range((len(values)-1).bit_length())]
    if adjoint:
        bits.reverse()
    count = 0
    for bit in bits:
        coefficient = power(delta, bit)  # Setup constant, not a state product.
        for index in range(len(values)):
            if index & bit:
                destination, source = ((index ^ bit, index) if adjoint else (index, index ^ bit))
                result[destination] ^= multiply(coefficient, result[source])
                count += 1
    return tuple(result), count


def monomial_maps(state):
    points = phase_geometry()[0]
    emission = tuple(tuple(power(u, j) for j in range(state)) for u in points)
    return emission, transpose(emission)


def additive_values(u, state):
    points = phase_geometry()[0]
    x2, x3 = multiply(u, u), power(u, 3)
    s2 = power(u, 4) ^ multiply(points[7], x2) ^ multiply(points[6], u)
    return (1, u, x2 ^ u, x3 ^ x2, s2, multiply(u, s2))[:state]


def additive_coefficients(values):
    q = tuple(values) + (0,)*(6-len(values))
    seven = phase_geometry()[0][7]
    return (q[0], q[1] ^ q[2] ^ q[3] ^ q[4] ^ q[5],
            q[2] ^ q[3] ^ multiply(seven, q[4]) ^ q[5],
            q[3] ^ multiply(seven, q[5]), q[4], q[5])[:len(values)]


def normalized_cost(length, dimension, state):
    if state not in (4, 5, 6):
        raise ValueError("cost model covers states 4, 5, and 6")
    if dimension == 24 and 32 <= length <= 40:
        outer_products = 34 + 7*(length-32)
    elif (length, dimension) == (32, 23):
        outer_products = 35
    elif (length, dimension) == (32, 22):
        outer_products = 41
    else:
        raise ValueError("unsupported outer cost model")
    phase_products = sum(i.bit_count() for i in range(state))
    step_zmm = 16 + 2*(state-4)  # Both FFT maps plus four diagonal products.
    step_xmm = 2*(state-4) + phase_products
    zmm = Fraction(step_zmm*length, 16*dimension) + Fraction(outer_products, 4*dimension)
    xmm = Fraction(step_xmm*length, 16*dimension)
    return dict(length=length, dimension=dimension, state=state,
                rate=Fraction(dimension, length), outer_products=outer_products,
                step_zmm=step_zmm, step_xmm=step_xmm,
                zmm_per_source=zmm, xmm_per_source=xmm,
                gfni_byte_equivalent=zmm+xmm/4,
                record_bytes_scatter=16*(3*Fraction(length, dimension)+1),
                record_bytes_region=16*(5*Fraction(length, dimension)+1))


def check():
    assert interpolation_row(32) == EXTRA_ROW
    span = {0}
    for coefficient in SCALAR_BASIS:
        span |= {value ^ coefficient for value in tuple(span)}
    assert len(span) == 256
    known = set(SCALAR_BASIS)
    for target, left, right in SCALAR_XORS:
        assert left in known and right in known and target == left ^ right
        assert target not in known
        known.add(target)
    for offset in range(8):
        row = interpolation_row(32+offset)
        assert row == tuple(EXTRA_ROW[column ^ offset] for column in range(24))
        for value in range(256):
            cost = Cost()
            assert extra_transpose(value, offset, cost) == tuple(multiply(c, value) for c in row)
            assert (cost.products, cost.xor2, cost.xor3) == (7, 23, 0)
    # GRS32/23 extends 23 source values to point 23, then uses GRS32/24.
    # This extra map has only one nontrivial product, shared by all outputs.
    short_row = interpolation_row(23, 23)
    span = {0, 1}
    for coefficient in short_row:
        span |= {value ^ coefficient for value in tuple(span)}
    assert short_row == (12,)*8 + (13,)*8 + (1,)*7
    assert len(span) == 4
    two_rows = tuple(interpolation_row(point, 22) for point in (22, 23))
    for row, point in zip(two_rows, (22, 23)):
        other_missing = 23 if point == 22 else 22
        assert row == tuple(multiply(12 if j < 8 else 13 if j < 16 else 1, j ^ other_missing)
                            for j in range(22))
    # Exact GF2 dimension bound: free input signals span two dimensions;
    # all required output forms together with them span nine dimensions.
    # Each fixed-scalar product introduces at most one new signal dimension.
    pairs = [1, 256] + [a | (b << 8) for a, b in zip(*two_rows)]
    pivots = {}
    for value in pairs:
        while value:
            pivot = value.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    assert len(pivots) == 9
    tables = tuple(tuple(tuple(multiply(coefficient, value) for value in range(256))
                         for coefficient in row) for row in two_rows)
    for t22 in range(256):
        for t23 in range(256):
            expected = tuple(tables[0][j][t22] ^ tables[1][j][t23] for j in range(22))
            assert two_missing_correction(t22, t23) == expected
    cost = Cost()
    two_missing_correction(137, 249, cost)
    assert cost.products == 7 and cost.boolean_instructions == 32

    rng = Random(20261005)
    for length in (33, 34, 36):
        rows = tuple(interpolation_row(point) for point in range(24, length))
        for _ in range(32):
            values = tuple(rng.randrange(256) for _ in rows)
            expected = tuple(xor_sum(multiply(row[column], value) for row, value in zip(rows, values))
                             for column in range(24))
            cost = Cost()
            assert extended_parity_transpose(values, cost) == expected
            assert cost.products == 34+7*(length-32)
    short_rows = tuple(interpolation_row(point, 23) for point in range(23, 32))
    for _ in range(64):
        values = tuple(rng.randrange(256) for _ in range(32))
        expected = tuple(values[j] ^ xor_sum(multiply(row[j], value)
                                            for row, value in zip(short_rows, values[23:]))
                         for j in range(23))
        cost = Cost()
        assert shortened_outer_transpose(values, cost) == expected
        assert cost.products == 35 and cost.boolean_instructions == 139
    two_short_rows = tuple(interpolation_row(point, 22) for point in range(22, 32))
    for _ in range(64):
        values = tuple(rng.randrange(256) for _ in range(32))
        expected = tuple(values[j] ^ xor_sum(multiply(row[j], value)
                                            for row, value in zip(two_short_rows, values[22:]))
                         for j in range(22))
        cost = Cost()
        assert twice_shortened_outer_transpose(values, cost) == expected
        assert cost.products == 41 and cost.boolean_instructions == 170

    for state, expected_products in ((4, 4), (5, 5), (6, 7)):
        emission, injection = monomial_maps(state)
        assert matmul(injection, emission) == tuple((0,)*state for _ in range(state))
        for delta in range(256):
            matrix = translation_matrix(state, delta)
            assert matmul(matrix, matrix) == tuple(tuple(int(i == j) for j in range(state)) for i in range(state))
            for coordinate in range(state):
                values = tuple(255 if i == coordinate else 0 for i in range(state))
                actual, count = translate(values, delta)
                transposed, dual_count = translate(values, delta, True)
                assert actual == matvec(matrix, values)
                assert transposed == matvec(transpose(matrix), values)
                assert count == dual_count == expected_products
            gamma = rng.randrange(256)
            assert matmul(matrix, translation_matrix(state, gamma)) == translation_matrix(state, delta ^ gamma)
        for _ in range(64):
            values = tuple(rng.randrange(256) for _ in range(state))
            basis = additive_coefficients(values)
            assert matvec(emission, values) == tuple(
                xor_sum(multiply(a, b) for a, b in zip(additive_values(u, state), basis))
                for u in phase_geometry()[0])
    assert normalized_cost(32, 24, 5)["gfni_byte_equivalent"] == 2
    assert normalized_cost(32, 24, 6)["gfni_byte_equivalent"] == Fraction(9, 4)
    assert normalized_cost(32, 22, 4)["gfni_byte_equivalent"] == Fraction(177, 88)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        check()
        print("extra parity, GRS32/23 and GRS32/22, state 4/5/6 translations/adjoints, CA=0, FFT basis and normalized costs pass")
    else:
        for params in ((32,24,4),(33,24,4),(34,24,4),(36,24,4),(32,24,5),(32,24,6),(32,23,4),(32,22,4)):
            row = normalized_cost(*params)
            print(" ".join(f"{key}={value}" for key, value in row.items()))
