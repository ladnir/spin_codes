"""Exact systematic GRS[16,12] parity-transpose circuit over GF(16).

This is a scalar reference and cost derivation, not a native implementation
or a distance certificate.  Source values occupy points 0..11 and parity
values occupy points 12..15 in GF(2)[X]/(X^4+X+1).  The optimized function
uses seven field multiplications, all by fixed subfield constants.  With
three-input XOR available, its other work is 26 Boolean instructions.
Adding the twelve systematic inputs costs twelve more XOR instructions.

Each four-by-four block of the parity transpose is XOR-circulant with
first row (a,b,c,d) and a+b+c+d=1.  In the group-algebra basis u=g+1,
v=h+1, where u^2=v^2=0, define

    s=x0+x1+x2+x3, z1=x1+x3, z2=x2+x3,
    r0=s,
    r1=z1+(b+d)*s,
    r2=z2+(c+d)*s,
    r3=x3+(b+d)*z2+(c+d)*z1+d*s.

Its output is (r0+r1+r2+r3, r1+r3, r2+r3, r3).  The three
blocks have (b+d,c+d,d)=(13,9,7),(5,15,12),(8,6,7).
The implementation below shares products across these blocks.

Because the identities use only field operations and fixed GF16 scalars,
they extend to compatible GF(2^m) symbol representations, m=8,32,64,128.
All coefficients must be embedded into the representation used by GFNI;
their abstract nibbles are not AES-field byte representations.
"""

from __future__ import annotations

import argparse


def multiply(a: int, b: int) -> int:
    """Multiply abstract GF16 nibbles, independently of the circuit."""
    result = 0
    for _ in range(4):
        if b & 1:
            result ^= a
        a = (a << 1) ^ (0x13 if a & 8 else 0)
        b >>= 1
    return result


def inverse(a: int) -> int:
    if not a:
        raise ZeroDivisionError("zero GF16 inverse")
    result = 1
    for _ in range(14):
        result = multiply(result, a)
    return result


def lagrange_parity() -> tuple[tuple[int, ...], ...]:
    """Construct the four-by-twelve parity matrix from interpolation."""
    rows = []
    for point in range(12, 16):
        row = []
        for source in range(12):
            numerator = denominator = 1
            for other in range(12):
                if other != source:
                    numerator = multiply(numerator, point ^ other)
                    denominator = multiply(denominator, source ^ other)
            row.append(multiply(numerator, inverse(denominator)))
        rows.append(tuple(row))
    return tuple(rows)


PARITY = (
    (2, 10, 14, 7, 7, 9, 3, 12, 8, 15, 1, 7),
    (10, 2, 7, 14, 9, 7, 12, 3, 15, 8, 7, 1),
    (14, 7, 2, 10, 3, 12, 7, 9, 1, 7, 8, 15),
    (7, 14, 10, 2, 12, 3, 9, 7, 7, 1, 15, 8),
)


def parity_transpose(x: tuple[int, int, int, int]) -> tuple[int, ...]:
    """Seven fixed-scalar products; 26 XOR/three-input-XOR operations."""
    x0, x1, x2, x3 = x
    # Three binary XORs and one ternary XOR.
    t, z1, z2 = x0 ^ x3, x1 ^ x3, x2 ^ x3
    s = t ^ z1 ^ z2

    # Three products and five binary XORs give every needed multiple of s.
    s12, s8, s6 = multiply(12, s), multiply(8, s), multiply(6, s)
    s13, s9, s7 = s ^ s12, s ^ s8, s ^ s6
    s5, s15 = s12 ^ s9, s8 ^ s7

    # Four products.  The third block's mixed term is recovered by XOR.
    a13, a5 = multiply(13, z2), multiply(5, z2)
    b9, b15 = multiply(9, z1), multiply(15, z1)
    # Five Boolean instructions: two for r3a, one each for mix/r3b/r3c.
    r3a = (a13 ^ b9 ^ s7) ^ x3
    mix = a5 ^ b15
    r3b = mix ^ x3 ^ s12
    r3c = r3a ^ mix

    # Four Boolean instructions per block: two for y0, one each for y1/y2.
    return (
        (t ^ s13 ^ s9) ^ r3a, z1 ^ s13 ^ r3a, z2 ^ s9 ^ r3a, r3a,
        (t ^ s5 ^ s15) ^ r3b, z1 ^ s5 ^ r3b, z2 ^ s15 ^ r3b, r3b,
        (t ^ s8 ^ s6) ^ r3c, z1 ^ s8 ^ r3c, z2 ^ s6 ^ r3c, r3c,
    )


def check_exhaustive() -> int:
    """Compare all 16^4 inputs with an independently interpolated matrix."""
    expected_matrix = lagrange_parity()
    assert expected_matrix == PARITY
    # Precomputation only accelerates the exhaustive test's dense oracle.
    products = tuple(tuple(multiply(a, b) for b in range(16)) for a in range(16))
    for packed in range(1 << 16):
        x = tuple((packed >> (4 * j)) & 15 for j in range(4))
        expected = tuple(
            products[expected_matrix[0][j]][x[0]]
            ^ products[expected_matrix[1][j]][x[1]]
            ^ products[expected_matrix[2][j]][x[2]]
            ^ products[expected_matrix[3][j]][x[3]]
            for j in range(12)
        )
        assert parity_transpose(x) == expected, (x, expected)
    return 1 << 16


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="check all 65536 parity inputs")
    args = parser.parse_args()
    if args.check:
        print(f"GRS[16,12]: {check_exhaustive()} exhaustive parity-transpose checks pass")
    else:
        print("7 fixed GF16 products + 26 XOR/ternlog operations per four parity inputs")
        print("Add 12 XOR operations for the systematic contribution")
