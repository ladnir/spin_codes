"""Small independent arithmetic oracle for the literal GF(2^128) lift.

There is no encoder or timing code here. ``--constants`` prints the pinned
AES embedding and the rank of the proposed sixteen-byte coordinate basis.
"""

import argparse

BITS = 128
MASK = (1 << BITS) - 1
MODULUS = (1 << BITS) | 0x87


def multiply(a, b):
    """Shift-and-reduce native polynomial multiplication."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        a = ((a << 1) & MASK) ^ (0x87 if a >> 127 else 0)
        b >>= 1
    return result


def polynomial_remainder(a, modulus):
    while a.bit_length() >= modulus.bit_length():
        a ^= modulus << (a.bit_length() - modulus.bit_length())
    return a


def multiply_oracle(a, b):
    """Unreduced schoolbook product, then independent polynomial division."""
    product = 0
    while b:
        bit = (b & -b).bit_length() - 1
        product ^= a << bit
        b &= b - 1
    return polynomial_remainder(product, MODULUS)


def power(a, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply_oracle(result, a)
        a = multiply_oracle(a, a)
        exponent >>= 1
    return result


def aes_multiply(a, b):
    product = 0
    for bit in range(8):
        if b >> bit & 1:
            product ^= a << bit
    return polynomial_remainder(product, 0x11B)


def aes_polynomial(a):
    a2 = multiply_oracle(a, a)
    a3 = multiply_oracle(a2, a)
    a4 = multiply_oracle(a2, a2)
    a8 = multiply_oracle(a4, a4)
    return a8 ^ a4 ^ a3 ^ a ^ 1


def derive_aes_root():
    """First AES root among successive powers of X^((2^128-1)/255)."""
    generator = power(2, MASK // 255)
    value = 1
    for exponent in range(255):
        if aes_polynomial(value) == 0:
            return value, generator, exponent
        value = multiply_oracle(value, generator)
    raise ValueError("the pinned native X projection contains no AES root")


def embedding(root):
    columns = tuple(power(root, bit) for bit in range(8))
    return tuple(transform(value, columns) for value in range(256))


def transform(value, columns):
    result = 0
    while value:
        bit = (value & -value).bit_length() - 1
        result ^= columns[bit]
        value &= value - 1
    return result


def byte_basis_columns(embedded):
    # Byte j is its AES coefficient times native X^j, for j=0,...,15.
    return tuple(multiply_oracle(embedded[1 << bit], 1 << byte)
                 for byte in range(16) for bit in range(8))


def inverse_columns(columns):
    pivots = {}
    for index, value in enumerate(columns):
        source = 1 << index
        while value:
            bit = value.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = value, source
                break
            pivot, preimage = pivots[bit]
            value ^= pivot
            source ^= preimage
        if not value:
            raise ValueError("dependent byte-coordinate basis")
    inverse = []
    for bit in range(len(columns)):
        value, source = 1 << bit, 0
        while value:
            pivot, preimage = pivots[value.bit_length() - 1]
            value ^= pivot
            source ^= preimage
        inverse.append(source)
    return tuple(inverse)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constants", action="store_true", required=True)
    parser.parse_args()
    root, generator, exponent = derive_aes_root()
    columns = byte_basis_columns(embedding(root))
    inverse = inverse_columns(columns)
    print(f"native_modulus=0x{MODULUS:x}")
    print(f"projected_generator=0x{generator:032x}")
    print(f"root_power={exponent}")
    print(f"aes_root=0x{root:032x}")
    print(f"byte_basis_rank={len(inverse)}")


if __name__ == "__main__":
    main()
