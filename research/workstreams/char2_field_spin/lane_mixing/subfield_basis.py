"""Explicit GF256-basis extension fields for the subfield encoder.

Byte i is the coefficient of Z**i. GF256 bytes use the AES polynomial.
Only the reference/tests need general extension multiplication. Encoding by
GF256 coefficients multiplies each byte independently in this declared basis.
"""

from functools import lru_cache


# Low coefficient first; all coefficients are AES-field bytes, not integers.
MODULI = {
    8: (0, 1),
    32: (0x3b, 0xd1, 0, 0xd4, 1),
    64: (0x57, 0x48, 0, 0, 0, 0, 0, 0x0f, 1),
    128: (0x63, 0xa5) + (0,) * 13 + (0x41, 1),
}


@lru_cache(maxsize=65536)
def gf_mul(a, b):
    if not (0 <= a < 256 and 0 <= b < 256):
        raise ValueError("GF256 coefficient")
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        if a & 256:
            a ^= 0x11b
        b >>= 1
    return result


def gf_pow(a, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = gf_mul(result, a)
        a = gf_mul(a, a)
        exponent >>= 1
    return result


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return tuple(a)


def add(a, b):
    out = list(a) + [0] * max(0, len(b) - len(a))
    for i, value in enumerate(b):
        out[i] ^= value
    return trim(out)


def remainder(a, modulus):
    a = list(trim(a))
    modulus = trim(modulus)
    if not modulus:
        raise ZeroDivisionError("zero polynomial")
    inv_lead = gf_pow(modulus[-1], 254)
    while len(a) >= len(modulus):
        factor = gf_mul(a[-1], inv_lead)
        shift = len(a) - len(modulus)
        for i, value in enumerate(modulus):
            a[shift + i] ^= gf_mul(factor, value)
        while a and a[-1] == 0:
            a.pop()
    return tuple(a)


def multiply(a, b, modulus):
    if not a or not b:
        return ()
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] ^= gf_mul(x, y)
    return remainder(out, modulus)


def square(a, modulus):
    if not a:
        return ()
    out = [0] * (2 * len(a) - 1)
    for i, x in enumerate(a):
        out[2 * i] = gf_mul(x, x)
    return remainder(out, modulus)


def frobenius(a, modulus):
    for _ in range(8):
        a = square(a, modulus)
    return a


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    if not a:
        return ()
    inv_lead = gf_pow(a[-1], 254)
    return tuple(gf_mul(x, inv_lead) for x in a)


def irreducible(modulus):
    """Exact Rabin test over GF256, for the power-of-two degrees used here."""
    modulus = trim(modulus)
    degree = len(modulus) - 1
    if degree < 1 or degree & (degree - 1) or modulus[-1] != 1:
        raise ValueError("a monic polynomial of power-of-two positive degree is required")
    if degree == 1:
        return True
    x = (0, 1)
    value = x
    for i in range(1, degree + 1):
        value = frobenius(value, modulus)
        if i == degree // 2 and gcd(add(value, x), modulus) != (1,):
            return False
    return value == x


def search_sparse(degree):
    for power in range(degree - 1, 0, -2):
        for constant in range(1, 256):
            candidate = [0] * (degree + 1)
            candidate[0], candidate[power], candidate[-1] = constant, 1, 1
            if irreducible(candidate):
                return tuple(candidate)
    # The base field contains many small subfields; restricting every
    # nonconstant coefficient to one can exclude all irreducible choices.
    from random import Random
    rng = Random(20261003 + degree)
    for _ in range(4096):
        candidate = [rng.randrange(1, 256)] + [0] * (degree - 1) + [1]
        candidate[1] = rng.randrange(1, 256)
        candidate[degree - 1] = rng.randrange(1, 256)
        if irreducible(candidate):
            return tuple(candidate)
    raise ValueError("no sparse polynomial in this search")


if __name__ == "__main__":
    for degree in (4, 8, 16):
        print(degree, search_sparse(degree), flush=True)
