#!/usr/bin/env python3
"""Exact b64/b128 causal-minor obstructions in the native AES tower.

No external packages or generated files. Run serially with python -B.
This checks an impulse-rank counterexample, not a code-distance failure.
"""


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def byte_mul(a, b):
    result = 0
    for _ in range(8):
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 256:
            a ^= 0x11b
    return result


def mul16(a, b):
    # Native u^2+u=0x20, direct four-product expansion.
    a0, a1 = a & 255, a >> 8
    b0, b1 = b & 255, b >> 8
    bd = byte_mul(a1, b1)
    low = byte_mul(a0, b0) ^ byte_mul(0x20, bd)
    high = byte_mul(a0, b1) ^ byte_mul(a1, b0) ^ bd
    return low | (high << 8)


def mul(a, b):
    # Native v^2+v=0x20*u, byte order [1,u,v,uv].
    a0, a1 = a & 65535, a >> 16
    b0, b1 = b & 65535, b >> 16
    bd = mul16(a1, b1)
    low = mul16(a0, b0) ^ mul16(0x2000, bd)
    high = mul16(a0, b1) ^ mul16(a1, b0) ^ bd
    return low | (high << 16)


def power(a, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = mul(result, a)
        a = mul(a, a)
        exponent >>= 1
    return result


def delta(a, count=1):
    for _ in range(count):
        a = mul(a, a) ^ a
    return a


def divide(a, b):
    require(b != 0, "division by zero")
    return mul(a, power(b, (1 << 32) - 2))


def xor_all(values):
    result = 0
    for value in values:
        result ^= value
    return result


def cantor_basis():
    # Match Transform256.h's byte roots, then AdditivePersistent32.h's
    # root with input bit0 fixed to zero. Gaussian rows are F2 equations.
    basis = [1]
    for _ in range(7):
        basis.append(next(x for x in range(256) if delta(x) == basis[-1]))
    equations = [0] * 32
    for column in range(1, 32):
        image = delta(1 << column)
        for row in range(32):
            if (image >> row) & 1:
                equations[row] |= 1 << column

    def root(target):
        rows = [equations[i] | (((target >> i) & 1) << 32)
                for i in range(32)]
        used = 0
        pivot_rows = {}
        for column in range(1, 32):
            pivot = next((i for i in range(used, 32)
                          if (rows[i] >> column) & 1), None)
            require(pivot is not None, "Delta columns lost independence")
            rows[used], rows[pivot] = rows[pivot], rows[used]
            for row in range(32):
                if row != used and ((rows[row] >> column) & 1):
                    rows[row] ^= rows[used]
            pivot_rows[column] = used
            used += 1
        require(all((row & 0xffffffff) or not (row >> 32) for row in rows),
                "Artin-Schreier target has no root")
        return sum(((rows[pivot_rows[column]] >> 32) & 1) << column
                   for column in range(1, 32))

    while len(basis) < 32:
        basis.append(root(basis[-1]))
    require(all(delta(basis[i]) == basis[i - 1] for i in range(1, 32)),
            "Cantor chain identity failed")
    return tuple(basis)


def coordinate_map(basis):
    pivots = {}
    for index, value in enumerate(basis):
        coordinates = 1 << index
        while value:
            bit = value.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = value, coordinates
                break
            pv, pc = pivots[bit]
            value ^= pv
            coordinates ^= pc
        require(value != 0, "Cantor basis dependence")

    def convert(value):
        result = 0
        while value:
            bit = value.bit_length() - 1
            require(bit in pivots, "point outside basis")
            pv, pc = pivots[bit]
            value ^= pv
            result ^= pc
        return result

    return convert


def inverse_gray(value):
    result = 0
    while value:
        result ^= value
        value >>= 1
    return result


def check_example(h, d, e, offset, eigenvalue, basis, coordinates):
    require(d and e and d != e and eigenvalue, "degenerate eigenspace")
    require(all(delta(x, h) == mul(eigenvalue, x)
                for x in (d, e, d ^ e)), "eigenvalue identity failed")
    points = [offset, offset ^ d, offset ^ e, offset ^ d ^ e]
    entries = []
    for point in points:
        coordinate = coordinates(point)
        require(coordinate < (1 << 21), "point is outside V21")
        quotient, lane = coordinate >> h, coordinate & ((1 << h) - 1)
        require(quotient != 0, "point is in the omitted predecessor coset")
        gray_index = inverse_gray(quotient)
        batch = gray_index - 1
        require(0 <= batch < (1 << (21 - h)), "batch outside K20 schedule")
        require((gray_index ^ (gray_index >> 1)) == quotient,
                "Gray inverse mismatch")
        representative = xor_all(basis[i + h]
                                 for i in range(32 - h) if (quotient >> i) & 1)
        base_point = xor_all(basis[i] for i in range(h) if (lane >> i) & 1)
        require(point == representative ^ base_point, "native point mismatch")
        entries.append((batch, lane, point))
    entries.sort()
    require(len({entry[0] for entry in entries}) == 4, "cosets are not distinct")
    require(entries[1][0] < entries[2][0], "causal split failed")
    sources, observations = entries[:2], entries[2:]
    matrix = []
    for target_batch, _, x in observations:
        row = []
        for source_batch, _, z in sources:
            require(source_batch < target_batch, "unexpected causal zero")
            difference = x ^ z
            require(delta(difference, h) != 0, "same-coset difference")
            row.append(divide(delta(difference, h), difference))
        matrix.append(row)
    require(all(value == eigenvalue for row in matrix for value in row),
            "impulse entries are not the same nonzero scalar")
    determinant = mul(matrix[0][0], matrix[1][1]) ^ mul(matrix[0][1], matrix[1][0])
    require(determinant == 0 and matrix[0][0] != 0, "matrix does not have rank one")
    print(f"b={1 << h}: rank-one noncausal-zero minor, "
          f"lambda=0x{eigenvalue:08x}, batches/lanes="
          f"{[(t, j) for t, j, _ in entries]}")


def main():
    require(delta(0x100) == 0x20, "native u relation failed")
    require(delta(0x10000) == 0x2000, "native v relation failed")
    basis = cantor_basis()
    coordinates = coordinate_map(basis)

    omega = next(x for x in range(2, 256) if byte_mul(x, x) ^ x == 1)
    d64 = next(x for x in range(1, 256) if delta(x, 6) != 0)
    e64 = byte_mul(omega, d64)
    lambda64 = divide(delta(d64, 6), d64)
    require(delta(0x100, 8) == 1 and delta(0x100, 9) == 0,
            "b64 offset prefix membership failed")
    check_example(6, d64, e64, 0x100, lambda64, basis, coordinates)

    zeta, mu = 0x02, 0x8d
    require(byte_mul(zeta, mu) == 1, "AES inverse mismatch")
    conjugates = [mu]
    for _ in range(7):
        conjugates.append(byte_mul(conjugates[-1], conjugates[-1]))
    require(xor_all(conjugates) == 1, "mu trace is not one")
    require(xor_all(byte_mul(conjugates[i], conjugates[j])
                    for i in range(8) for j in range(i)) == 0,
            "mu subtrace is not zero")
    t = next(x for x in range(256) if byte_mul(x, x) ^ x == mu ^ 0x20)
    w = 0x100 ^ t
    require(delta(w) == mu and w >= 256, "native quadratic root failed")
    d128, e128 = mu, mul(mu, w)
    require(delta(0x10000, 16) == 1 and delta(0x10000, 17) == 0,
            "b128 offset prefix membership failed")
    check_example(7, d128, e128, 0x10000, zeta, basis, coordinates)
    print("PASS: exact b64/b128 impulse obstructions in the full-Cantor K20 prefix")


if __name__ == "__main__":
    main()
