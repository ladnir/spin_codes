"""Exact GRS[32,24] parity transpose and a fixed-coefficient cost circuit.

The field is AES GF256, with modulus X^8+X^4+X^3+X+1.  Source values
occupy points 0..23; parity values occupy points 24..31.  These are four
cosets of the additive subspace {0,..,7}, not cosets of a GF8 subfield.

This scalar reference uses 34 fixed GF256 products and 114 Boolean
instructions when an XOR of three operands costs one instruction.
Systematic contributions require 24 additional XORs.  No native timing
or minimum-distance certificate is claimed.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from random import Random


def multiply(a: int, b: int) -> int:
    result = 0
    for _ in range(8):
        if b & 1:
            result ^= a
        a = (a << 1) ^ (0x11B if a & 128 else 0)
        b >>= 1
    return result


def inverse(a: int) -> int:
    if not a:
        raise ZeroDivisionError("zero GF256 inverse")
    result, power, exponent = 1, a, 254
    while exponent:
        if exponent & 1:
            result = multiply(result, power)
        power = multiply(power, power)
        exponent >>= 1
    return result


def lagrange_parity() -> tuple[tuple[int, ...], ...]:
    """Independent interpolation matrix, with no circuit constants."""
    rows = []
    for point in range(24, 32):
        row = []
        for source in range(24):
            numerator = denominator = 1
            for other in range(24):
                if other != source:
                    numerator = multiply(numerator, point ^ other)
                    denominator = multiply(denominator, source ^ other)
            row.append(multiply(numerator, inverse(denominator)))
        rows.append(tuple(row))
    return tuple(rows)


# First rows of the three symmetric XOR-circulant blocks.
CIRCULANT = (
    (208, 230, 192, 134, 42, 68, 78, 81),
    (162, 127, 147, 40, 38, 203, 227, 105),
    (52, 172, 158, 17, 209, 32, 250, 29),
)

# Superset zeta transforms of CIRCULANT.  The zero coefficient is one.
GAMMA = (
    (1, 117, 89, 215, 113, 21, 31, 81),
    (1, 245, 49, 65, 103, 162, 138, 105),
    (1, 128, 104, 12, 22, 61, 231, 29),
)

# For each transformed input Z_m, first compute the fixed products in
# "basis".  Then each (target,a,b) entry gives target*Z_m=a*Z_m+b*Z_m.
# Coefficient one is always available without a multiplication.
SCALAR_PROGRAMS = (
    ((61, 117, 12, 22, 89, 103, 215), ((49, 61, 12), (113, 22, 103), (48, 1, 49), (40, 89, 113), (231, 215, 48), (104, 89, 49), (121, 117, 12), (128, 103, 231), (162, 117, 215), (21, 40, 61), (29, 104, 117), (30, 121, 103), (31, 30, 1), (65, 48, 113), (81, 121, 40), (105, 48, 89), (138, 162, 40), (245, 128, 117))),
    ((22, 104, 49, 113, 138, 31, 231), ((89, 49, 104), (103, 113, 22))),
    ((128, 22, 21, 162, 61, 117, 113), ((103, 113, 22), (245, 117, 128))),
    ((22, 113), ((103, 113, 22),)),
    ((128, 49, 65, 215, 245, 12, 89), ((104, 89, 49), (117, 245, 128))),
    ((49, 89), ((104, 89, 49),)),
    ((117, 245), ((128, 245, 117),)),
    ((), ()),
)


@dataclass
class Cost:
    products: int = 0
    xor2: int = 0
    xor3: int = 0

    @property
    def boolean_instructions(self):
        return self.xor2 + self.xor3


def mobius(values, cost=None):
    """The same characteristic-two superset transform in both directions."""
    result = list(values)
    for bit in (1, 2, 4):
        for mask in range(8):
            if not mask & bit:
                result[mask] ^= result[mask | bit]
                if cost is not None:
                    cost.xor2 += 1
    return result


def xor_sum(terms, cost=None):
    """Accumulate using three-input XOR whenever three terms remain."""
    result = terms[0]
    index = 1
    while index+1 < len(terms):
        result ^= terms[index] ^ terms[index+1]
        if cost is not None:
            cost.xor3 += 1
        index += 2
    if index < len(terms):
        result ^= terms[index]
        if cost is not None:
            cost.xor2 += 1
    return result


def parity_transpose(values, cost=None):
    """Return 24 parity contributions from eight GF256 parity inputs."""
    if len(values) != 8:
        raise ValueError("eight parity inputs required")
    transformed = mobius(values, cost)
    products = []
    for value, (basis, operations) in zip(transformed, SCALAR_PROGRAMS):
        available = {1: value}
        for coefficient in basis:
            available[coefficient] = multiply(coefficient, value)
            if cost is not None:
                cost.products += 1
        for target, left, right in operations:
            available[target] = available[left] ^ available[right]
            if cost is not None:
                cost.xor2 += 1
        products.append(available)
    output = []
    for gamma in GAMMA:
        result = []
        for mask in range(8):
            terms = [products[mask ^ subset][gamma[subset]]
                     for subset in range(8) if subset & mask == subset]
            result.append(xor_sum(terms, cost))
        output.extend(mobius(result, cost))
    return tuple(output)


def dense_transpose(matrix, values):
    return tuple(xor_sum([multiply(matrix[row][column], values[row])
                          for row in range(8)])
                 for column in range(24))


def check():
    matrix = lagrange_parity()
    for block in range(3):
        assert tuple(matrix[0][8*block:8*block+8]) == CIRCULANT[block]
        assert tuple(mobius(CIRCULANT[block])) == GAMMA[block]
        assert GAMMA[block][0] == 1
        for row in range(8):
            for column in range(8):
                assert matrix[row][8*block+column] == CIRCULANT[block][row ^ column]
    for basis, operations in SCALAR_PROGRAMS:
        known = {1, *basis}
        for target, left, right in operations:
            assert left in known and right in known and target == left ^ right
            assert target not in known
            known.add(target)

    assert parity_transpose((0,)*8) == (0,)*24
    checks = 1
    # All nonzero scalar multiples of each coordinate basis vector.
    # The circuit and matrix are GF256-linear, so these include a basis.
    for coordinate in range(8):
        for scalar in range(1, 256):
            values = [0]*8
            values[coordinate] = scalar
            expected = tuple(multiply(matrix[coordinate][column], scalar)
                             for column in range(24))
            assert parity_transpose(values) == expected
            checks += 1
    rng = Random(20261004)
    for _ in range(256):
        values = tuple(rng.randrange(256) for _ in range(8))
        cost = Cost()
        assert parity_transpose(values, cost) == dense_transpose(matrix, values)
        assert (cost.products, cost.xor2, cost.xor3) == (34, 96, 18)
        checks += 1
    return checks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        print(f"GRS[32,24]: {check()} exact matrix/circuit checks pass")
    else:
        cost = Cost()
        parity_transpose(tuple(range(8)), cost)
        print(f"{cost.products} fixed GF256 products + {cost.boolean_instructions} XOR/ternlog instructions")
        print(f"Boolean breakdown: {cost.xor2} binary XOR + {cost.xor3} ternary XOR")
        print("Add 24 XOR instructions for the systematic contribution")
