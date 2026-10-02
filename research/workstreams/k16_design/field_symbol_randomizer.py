"""Exact transitive symbol randomizers and their nine-map binary adjoints.

The byte field is GF(2)[z]/(z^8+z^4+z^3+z+1).  Set mu=0x20 and adjoin
u with u^2+u+mu=0, then v with v^2+v+mu*u=0.  Coordinates are the four
bytes (1,u,v,uv), least significant first.  Both quadratic polynomials are
irreducible: Tr_GF256/GF2(mu)=1 and Tr_GF65536/GF2(mu*u)=1.

Sample an independent uniform nonzero 32-bit scalar for every outer symbol.
Multiplication fixes zero and sends each fixed nonzero symbol uniformly to
all nonzero symbols.  Thus it preserves the fixed-message symbol law used by
the first-moment outer proof.  It is not a uniform GL32 matrix or seed-compatible
replacement; no joint law for different messages is asserted.

Arithmetic acts on code coordinates, not on the 128-bit payload elements.
"""
import rs_maps

MU = 0x20
NU = MU << 8
THETA = NU << 16


def mul8(a, b):
    return rs_maps.multiply(a, b)


def mul16(a, b):
    if not 0 <= a < 1 << 16 or not 0 <= b < 1 << 16:
        raise ValueError("16-bit tower coordinates required")
    a0, a1, b0, b1 = a & 255, a >> 8, b & 255, b >> 8
    lo = mul8(a0, b0) ^ mul8(MU, mul8(a1, b1))
    hi = mul8(a0, b1) ^ mul8(a1, b0) ^ mul8(a1, b1)
    return lo | (hi << 8)


def mul32(a, b):
    if not 0 <= a < 1 << 32 or not 0 <= b < 1 << 32:
        raise ValueError("32-bit tower coordinates required")
    a0, a1, b0, b1 = a & 65535, a >> 16, b & 65535, b >> 16
    lo = mul16(a0, b0) ^ mul16(NU, mul16(a1, b1))
    hi = mul16(a0, b1) ^ mul16(a1, b0) ^ mul16(a1, b1)
    return lo | (hi << 16)


def mul64(a, b):
    if not 0 <= a < 1 << 64 or not 0 <= b < 1 << 64:
        raise ValueError("64-bit tower coordinates required")
    a0, a1, b0, b1 = a & 0xffffffff, a >> 32, b & 0xffffffff, b >> 32
    lo = mul32(a0, b0) ^ mul32(THETA, mul32(a1, b1))
    hi = mul32(a0, b1) ^ mul32(a1, b0) ^ mul32(a1, b1)
    return lo | (hi << 32)


def power(a, exponent, multiply=mul32):
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, a)
        a = multiply(a, a)
        exponent >>= 1
    return result


def trace(a, degree, multiply):
    result = 0
    for _ in range(degree):
        result ^= a
        a = multiply(a, a)
    return result


def coefficients16(scalar):
    c0, c1 = scalar & 255, scalar >> 8
    return c0, mul8(MU, c1), c0 ^ c1


def _coefficients32(scalar):
    c0, c1 = scalar & 65535, scalar >> 16
    return (coefficients16(c0) + coefficients16(mul16(NU, c1))
            + coefficients16(c0 ^ c1))


def coefficients32(scalar):
    """Nine byte multipliers, with tower reductions folded into constants."""
    if not 0 < scalar < 1 << 32:
        raise ValueError("uniform nonzero scalar required")
    return _coefficients32(scalar)


def coefficients64(scalar):
    """Twenty-seven byte multipliers in the third quadratic tower."""
    if not 0 < scalar < 1 << 64:
        raise ValueError("uniform nonzero scalar required")
    c0, c1 = scalar & 0xffffffff, scalar >> 32
    return (_coefficients32(c0) + _coefficients32(mul32(THETA, c1))
            + _coefficients32(c0 ^ c1))


def affine(value, matrix):
    return sum((((matrix >> (8 * (7 - bit))) & value).bit_count() & 1) << bit
               for bit in range(8))


def adjoint_matrix8(scalar):
    return sum(mul8(scalar, 1 << bit) << (8 * (7 - bit)) for bit in range(8))


def adjoint_coefficients(scalar):
    return tuple(adjoint_matrix8(c) for c in coefficients32(scalar))


def transpose16(y0, y1, coeff):
    p0 = affine(y0 ^ y1, coeff[0])
    p1 = affine(y0, coeff[1])
    p2 = affine(y1, coeff[2])
    return p0 ^ p2, p1 ^ p2


def _transpose32(value, coeff):
    y0, y1, y2, y3 = ((value >> (8 * j)) & 255 for j in range(4))
    p0, p1 = transpose16(y0 ^ y2, y1 ^ y3, coeff[:3])
    q0, q1 = transpose16(y0, y1, coeff[3:6])
    r0, r1 = transpose16(y2, y3, coeff[6:])
    return (p0 ^ r0) | ((p1 ^ r1) << 8) | ((q0 ^ r0) << 16) | ((q1 ^ r1) << 24)


def transpose32(value, scalar):
    """Binary-coordinate adjoint, not ordinary field multiplication."""
    return _transpose32(value, adjoint_coefficients(scalar))


def adjoint_coefficients64(scalar):
    return tuple(adjoint_matrix8(c) for c in coefficients64(scalar))


def transpose64(value, scalar):
    coeff = adjoint_coefficients64(scalar)
    y0, y1 = value & 0xffffffff, value >> 32
    p = _transpose32(y0 ^ y1, coeff[:9])
    q = _transpose32(y0, coeff[9:18])
    r = _transpose32(y1, coeff[18:])
    return (p ^ r) | ((q ^ r) << 32)


def adjoint_rows(scalar):
    """Row masks of the binary transpose matrix: forward basis images."""
    return tuple(mul32(scalar, 1 << bit) for bit in range(32))


def adjoint_rows64(scalar):
    return tuple(mul64(scalar, 1 << bit) for bit in range(64))


def multiply16_factored(x0, x1, coeff):
    p0 = mul8(x0, coeff[0])
    p1 = mul8(x1, coeff[1])
    p2 = mul8(x0 ^ x1, coeff[2])
    return p0 ^ p1, p0 ^ p2


def multiply32_factored(value, scalar):
    coeff = coefficients32(scalar)
    x0, x1, x2, x3 = ((value >> (8 * j)) & 255 for j in range(4))
    p0, p1 = multiply16_factored(x0, x1, coeff[:3])
    q0, q1 = multiply16_factored(x2, x3, coeff[3:6])
    r0, r1 = multiply16_factored(x0 ^ x2, x1 ^ x3, coeff[6:])
    return (p0 ^ q0) | ((p1 ^ q1) << 8) | ((p0 ^ r0) << 16) | ((p1 ^ r1) << 24)


def multiply_rows32(scalar):
    """Ordinary multiplication rows for an adjoint-multiplier forward family."""
    columns = adjoint_rows(scalar)
    return tuple(sum(((column >> row) & 1) << j for j, column in enumerate(columns))
                 for row in range(32))


def _multiply32_coefficients(value, coeff):
    x0, x1, x2, x3 = ((value >> (8 * j)) & 255 for j in range(4))
    p0, p1 = multiply16_factored(x0, x1, coeff[:3])
    q0, q1 = multiply16_factored(x2, x3, coeff[3:6])
    r0, r1 = multiply16_factored(x0 ^ x2, x1 ^ x3, coeff[6:])
    return (p0 ^ q0) | ((p1 ^ q1) << 8) | ((p0 ^ r0) << 16) | ((p1 ^ r1) << 24)


def multiply64_factored(value, scalar):
    coeff = coefficients64(scalar)
    x0, x1 = value & 0xffffffff, value >> 32
    p = _multiply32_coefficients(x0, coeff[:9])
    q = _multiply32_coefficients(x1, coeff[9:18])
    r = _multiply32_coefficients(x0 ^ x1, coeff[18:])
    return (p ^ q) | ((p ^ r) << 32)


def multiply_rows64(scalar):
    columns = adjoint_rows64(scalar)
    return tuple(sum(((column >> row) & 1) << j for j, column in enumerate(columns))
                 for row in range(64))
