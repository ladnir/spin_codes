"""Exact polynomial-basis field arithmetic and small-coefficient shears.

This is a research oracle, not an optimized encoder. Bit i represents X**i.
"""

MODULI = {
    8: (1 << 8) | 0x1B,
    32: (1 << 32) | 0x8D,
    64: (1 << 64) | 0x1B,
    128: (1 << 128) | 0x87,
}


def carryless_product(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


def polynomial_remainder(a, modulus):
    degree = modulus.bit_length() - 1
    while a.bit_length() - 1 >= degree:
        a ^= modulus << (a.bit_length() - 1 - degree)
    return a


def polynomial_gcd(a, b):
    while b:
        a, b = b, polynomial_remainder(a, b)
    return a


def multiply(a, b, modulus):
    return polynomial_remainder(carryless_product(a, b), modulus)


def frobenius(a, exponent, modulus):
    """Return a**(2**exponent), using only exact polynomial operations."""
    for _ in range(exponent):
        a = multiply(a, a, modulus)
    return a


def irreducible(modulus):
    """Rabin's criterion over GF(2), including every prime divisor of degree."""
    degree = modulus.bit_length() - 1
    if degree < 1:
        return False
    if degree == 1:
        return True
    if not modulus & 1:
        return False
    factors = set()
    remaining = degree
    divisor = 2
    while divisor * divisor <= remaining:
        while remaining % divisor == 0:
            factors.add(divisor)
            remaining //= divisor
        divisor += 1
    if remaining > 1:
        factors.add(remaining)
    if frobenius(2, degree, modulus) != 2:
        return False
    return all(polynomial_gcd(frobenius(2, degree // p, modulus) ^ 2,
                              modulus) == 1 for p in factors)


def multiply_byte_horner(value, coefficient, modulus):
    """At most seven reduced shifts and eight XORs; coefficient is a byte."""
    if not 0 <= coefficient < 256:
        raise ValueError("coefficient must be a byte polynomial")
    degree = modulus.bit_length() - 1
    mask = (1 << degree) - 1
    if not 0 <= value <= mask:
        raise ValueError("value does not fit field representation")
    result = value if coefficient & 128 else 0
    for bit in range(6, -1, -1):
        carry = result >> (degree - 1)
        result = (result << 1) & mask
        if carry:
            result ^= modulus ^ (1 << degree)
        if coefficient & (1 << bit):
            result ^= value
    return result


def multiply_byte_onefold(value, coefficient, modulus):
    """Exact single-fold reduction for the selected m>=32 sparse moduli."""
    degree = modulus.bit_length() - 1
    tail = modulus ^ (1 << degree)
    if degree < 8 or not 0 <= coefficient < 256:
        raise ValueError("invalid byte product geometry")
    if not 0 <= value < (1 << degree):
        raise ValueError("value does not fit field representation")
    # The unreduced high part has degree at most six. Its product with the
    # modulus tail must already fit below degree m for one fold to suffice.
    if 6 + tail.bit_length() - 1 >= degree:
        raise ValueError("this modulus needs additional reduction")
    product = carryless_product(value, coefficient)
    return (product & ((1 << degree) - 1)) ^ carryless_product(product >> degree, tail)


def two_shears(state, a, b, modulus):
    x, y = state
    u = x ^ multiply(a, y, modulus)
    v = y ^ multiply(b, u, modulus)
    return u, v


def three_shears(state, a, b, c, modulus):
    u, v = two_shears(state, a, b, modulus)
    return u ^ multiply(c, v, modulus), v


def shear_matrix(a, b, c, modulus):
    """Row-major U_c L_b U_a; det=1 in characteristic two."""
    ab = multiply(a, b, modulus)
    bc = multiply(b, c, modulus)
    return 1 ^ bc, a ^ c ^ multiply(ab, c, modulus), b, 1 ^ ab
