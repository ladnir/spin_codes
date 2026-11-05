"""Generic GF256-polynomial oracle for the sparse-byte GF128 primitive."""

from field128_oracle import aes_multiply, polynomial_remainder

MODULUS = (0x32, 0, 0, 1, 0, 1) + (0,) * 10 + (1,)
MASK = (1 << 128) - 1


def trim(p):
    while p and not p[-1]:
        p.pop()
    return p


def byte_power(a, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = aes_multiply(result, a)
        a = aes_multiply(a, a)
        exponent >>= 1
    return result


def remainder(a, modulus):
    """Generic polynomial long division over AES GF256."""
    a = trim(list(a))
    inverse = byte_power(modulus[-1], 254)
    while len(a) >= len(modulus):
        offset = len(a) - len(modulus)
        scale = aes_multiply(a[-1], inverse)
        for j, coefficient in enumerate(modulus):
            if coefficient:
                a[offset + j] ^= aes_multiply(scale, coefficient)
        trim(a)
    return a


def product(a, b):
    if not a or not b:
        return []
    result = [0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        if left:
            for j, right in enumerate(b):
                if right:
                    result[i + j] ^= aes_multiply(left, right)
    return trim(result)


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a


def square_mod(a, modulus):
    result = [0] * (2 * len(a) - 1)
    for i, coefficient in enumerate(a):
        result[2 * i] = aes_multiply(coefficient, coefficient)
    return remainder(result, modulus)


def irreducible():
    # Rabin's criterion over q=256 for degree16. Its only prime divisor is
    # 2, so Y^(q^16)=Y and gcd(Y^(q^8)-Y,P)=1 are sufficient and necessary.
    value = [0, 1]
    for squaring in range(1, 129):
        value = square_mod(value, MODULUS)
        if squaring == 64:
            middle = value + [0] * max(0, 2 - len(value))
            middle[1] ^= 1
            if len(gcd(list(MODULUS), trim(middle))) != 1:
                return False
    return value == [0, 1]


def coefficients(value):
    return [(value >> (8 * i)) & 255 for i in range(16)]


def multiply(a, b):
    p = remainder(product(coefficients(a), coefficients(b)), MODULUS)
    return sum(c << (8 * i) for i, c in enumerate(p))


def bytewise(value, coefficient):
    return sum(aes_multiply(c, coefficient) << (8 * i)
               for i, c in enumerate(coefficients(value)))


def multiply_spaced(a, b):
    """Model the independent spaced CLMUL and two-fold SIMD algorithm."""
    left = sum(c << (16 * i) for i, c in enumerate(coefficients(a)))
    right = sum(c << (16 * i) for i, c in enumerate(coefficients(b)))
    raw = 0
    while right:
        bit = (right & -right).bit_length() - 1
        raw ^= left << bit
        right &= right - 1
    reduced = sum(polynomial_remainder((raw >> (16 * i)) & 65535, 0x11B) << (8 * i)
                  for i in range(32))
    low, high = reduced & MASK, reduced >> 128
    overflow = (high >> (8 * 13)) ^ (high >> (8 * 11))
    first = bytewise(high, 0x32) ^ ((high << 24) & MASK) ^ ((high << 40) & MASK)
    second = bytewise(overflow, 0x32) ^ (overflow << 24) ^ (overflow << 40)
    return low ^ first ^ second
