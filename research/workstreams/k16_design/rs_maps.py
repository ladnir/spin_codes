"""Explicit systematic RS[8,4] map for the K16 candidate, not a SPIN certificate.

Bytes are GF(256) elements modulo x^8+x^4+x^3+x+1 (0x11b), matching GFNI
byte multiplication. A row contains four message bytes. The unique polynomial
of degree below four taking those values at 0,1,2,3 is evaluated at 0,...,7.
Four independent rows form a 128-to-256-bit binary outer group. Its eight
four-byte symbol positions are the proposed independent GL32 mixing blocks.
Arithmetic here operates on code coordinates, not on the 128-bit payloads.
"""
from itertools import combinations


def multiply(a, b):
    if not 0 <= a < 256 or not 0 <= b < 256:
        raise ValueError('GF256 byte operands required')
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        if a & 256:
            a ^= 0x11b
        b >>= 1
    return result


def inverse(a):
    if not 1 <= a < 256:
        raise ValueError('nonzero GF256 byte required')
    result, exponent = 1, 254
    while exponent:
        if exponent & 1:
            result = multiply(result, a)
        a = multiply(a, a)
        exponent >>= 1
    return result


def systematic_generator():
    """Return eight output rows, each a four-coefficient GF256 tuple."""
    generator = []
    for x in range(8):
        row = []
        for i in range(4):
            numerator, denominator = 1, 1
            for j in range(4):
                if i != j:
                    numerator = multiply(numerator, x ^ j)
                    denominator = multiply(denominator, i ^ j)
            row.append(multiply(numerator, inverse(denominator)))
        generator.append(tuple(row))
    return tuple(generator)


def encode_row(message):
    if len(message) != 4 or any(not 0 <= x < 256 for x in message):
        raise ValueError('four GF256 message bytes required')
    output = []
    for coefficients in systematic_generator():
        value = 0
        for coefficient, symbol in zip(coefficients, message):
            value ^= multiply(coefficient, symbol)
        output.append(value)
    return tuple(output)


def adjoint_multiply(coefficient, value):
    """Binary-coordinate adjoint of multiplication, not field multiplication.

    Polynomial-basis byte coordinates use the ordinary binary dot product.
    A GFNI transposed kernel uses the transposed 8-by-8 affine matrix here.
    """
    if not 0 <= coefficient < 256 or not 0 <= value < 256:
        raise ValueError('GF256 byte operands required')
    return sum(((multiply(coefficient, 1 << bit) & value).bit_count() & 1) << bit
               for bit in range(8))


def parity_factored(x, *, transpose=False):
    """Apply the symmetric four-by-four parity map with five byte products.

    This is an exact factorization for evaluation points 0,...,7 and 0x11b,
    not a replacement ensemble. The binary transpose substitutes adjoint
    byte maps: symmetry over GF256 does not imply binary symmetry.
    """
    if len(x) != 4 or any(not 0 <= value < 256 for value in x):
        raise ValueError('four GF256 bytes required')
    action = adjoint_multiply if transpose else multiply
    total = x[0] ^ x[1] ^ x[2] ^ x[3]
    u, v = action(6, total), action(8, total)
    common = action(26, total) ^ action(6, x[1] ^ x[3]) ^ action(8, x[2] ^ x[3])
    return (x[0] ^ common, x[1] ^ common ^ u,
            x[2] ^ common ^ v, x[3] ^ common ^ u ^ v)


def transpose_row(values):
    if len(values) != 8 or any(not 0 <= value < 256 for value in values):
        raise ValueError('eight GF256 coordinate bytes required')
    parity = parity_factored(values[4:], transpose=True)
    return tuple(a ^ b for a, b in zip(values[:4], parity))


def rank(rows):
    work = [list(row) for row in rows]
    pivots = 0
    for column in range(len(work[0])):
        pivot = next((r for r in range(pivots, len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        work[pivots], work[pivot] = work[pivot], work[pivots]
        scale = inverse(work[pivots][column])
        work[pivots] = [multiply(value, scale) for value in work[pivots]]
        for r in range(len(work)):
            if r != pivots:
                scale = work[r][column]
                work[r] = [a ^ multiply(scale, b) for a, b in zip(work[r], work[pivots])]
        pivots += 1
    return pivots


def verify():
    """Check systematic form, all 70 erasure patterns, and binary linearity."""
    generator = systematic_generator()
    assert generator[:4] == tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
    assert all(rank([generator[i] for i in subset]) == 4 for subset in combinations(range(8), 4))
    assert all(generator[i + 4][j] == generator[j + 4][i] for i in range(4) for j in range(4))
    for a in range(1, 256):
        assert multiply(a, inverse(a)) == 1
    encoded_basis = []
    for i in range(32):
        basis = [0] * 4
        basis[i // 8] = 1 << (i % 8)
        image = encode_row(basis)
        encoded_basis.append(image)
        assert image[:4] == tuple(basis)
        assert parity_factored(basis) == image[4:]
        assert parity_factored(parity_factored(basis)) == tuple(basis)
        for j in range(32):
            other = [0] * 4
            other[j // 8] = 1 << (j % 8)
            assert encode_row([a ^ b for a, b in zip(basis, other)]) == tuple(
                a ^ b for a, b in zip(image, encode_row(other)))
            combined = [a ^ b for a, b in zip(basis, other)]
            assert parity_factored(combined) == encode_row(combined)[4:]
    for bit in range(64):
        value = [0] * 8
        value[bit // 8] = 1 << (bit % 8)
        transposed = transpose_row(value)
        for input_bit, image in enumerate(encoded_basis):
            assert ((image[bit // 8] >> (bit % 8)) & 1) == (
                (transposed[input_bit // 8] >> (input_bit % 8)) & 1)
    return generator


if __name__ == '__main__':
    print('RS[8,4] GF256 generator:', verify())
    print('PASS: systematic map, 70 nonsingular four-symbol subsets, binary linearity, five-product factor and full binary adjoint')
