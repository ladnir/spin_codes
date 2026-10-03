"""Independent, scalar algebra for the byte-packet candidate; not a certificate.

No experimental C++ code or existing proof evaluator is imported. Integers
encode binary coordinates little-endian; payload dot products are over GF(2).
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations
from random import Random


def multiply(a, b, *, degree=8, modulus=0x11B):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a >> degree:
            a ^= modulus
    return result


@lru_cache(None)
def multiplication_table():
    return tuple(tuple(multiply(a, b) for b in range(256)) for a in range(256))


def dot(x, y):
    return (x & y).bit_count() & 1


def adjoint_multiply(c, x):
    """Binary transpose of multiplication by c in the AES polynomial basis."""
    return sum(dot(multiply(c, 1 << bit), x) << bit for bit in range(8))


def gfni_adjoint_matrix(c):
    """Conventional GFNI matrix operand: high byte is output row zero."""
    return sum(multiply(c, 1 << bit) << (8 * (7 - bit)) for bit in range(8))


def gfni_affine(x, matrix):
    return sum(dot(x, (matrix >> (8 * (7 - bit))) & 255) << bit
               for bit in range(8))


def pack_payload_bytes(raw):
    """Literal result of Fast.cpp's VBMI gather followed by GFNI transpose.

    Eight 128-bit records enter. Within each payload byte, the packed lanes
    enumerate payload bits seven through zero, not zero through seven.
    """
    if len(raw) != 128:
        raise ValueError('eight 128-bit records required')
    result = bytearray(128)
    for half in range(2):
        for byte in range(8):
            matrix = sum(raw[16 * (7 - coordinate) + 8 * half + byte] << (8 * coordinate)
                         for coordinate in range(8))
            for lane in range(8):
                result[64 * half + 8 * byte + lane] = gfni_affine(1 << (7 - lane), matrix)
    return bytes(result)


def unpack_payload_bytes(packed):
    if len(packed) != 128:
        raise ValueError('two 64-byte bitplane vectors required')
    return bytes(sum(((packed[64 * (byte // 8) + 8 * (byte % 8) + 7 - bit] >> coordinate) & 1) << bit
                     for bit in range(8))
                 for coordinate in range(8) for byte in range(16))


def rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def expansion(state, *, points=tuple(range(8)), scales=None):
    a, b = state & 255, state >> 8
    table = multiplication_table()
    if scales is None:
        scales = (1,) * len(points)
    return sum(table[c][a ^ table[alpha][b]] << (8 * i)
               for i, (alpha, c) in enumerate(zip(points, scales)))


def feedback(word, *, points=tuple(range(8)), scales=None):
    table = multiplication_table()
    a = b = 0
    if scales is None:
        scales = (1,) * len(points)
    for i, (alpha, c) in enumerate(zip(points, scales)):
        value = table[inverse(c)][(word >> (8 * i)) & 255]
        a ^= value
        b ^= table[alpha][value]
    return a | (b << 8)


@lru_cache(None)
def inverse(value):
    if not 0 < value < 256:
        raise ValueError('nonzero field byte required')
    return next(candidate for candidate in range(1, 256)
                if multiply(value, candidate) == 1)


def apply_columns(columns, word):
    result = 0
    while word:
        bit = (word & -word).bit_length() - 1
        result ^= columns[bit]
        word &= word - 1
    return result


def apply_transpose(columns, word):
    return sum(dot(image, word) << bit for bit, image in enumerate(columns))


def local_maps(*, points=tuple(range(8)), scales=None):
    a_columns = tuple(expansion(1 << bit, points=points, scales=scales)
                      for bit in range(16))
    c_columns = tuple(feedback(1 << bit, points=points, scales=scales)
                      for bit in range(8 * len(points)))
    return a_columns, c_columns


def matrix_action(matrix, state):
    aa, ab, ba, bb = matrix
    a, b = state & 255, state >> 8
    return (multiply(aa, a) ^ multiply(ab, b)) | (
        (multiply(ba, a) ^ multiply(bb, b)) << 8)


def matrix_adjoint(matrix, state):
    aa, ab, ba, bb = matrix
    a, b = state & 255, state >> 8
    return (adjoint_multiply(aa, a) ^ adjoint_multiply(ba, b)) | (
        (adjoint_multiply(ab, a) ^ adjoint_multiply(bb, b)) << 8)


def sample_matrix(random):
    while True:
        matrix = tuple(random.randrange(256) for _ in range(4))
        aa, ab, ba, bb = matrix
        if multiply(aa, bb) ^ multiply(ab, ba):
            return matrix


def inner_forward(words, matrices):
    if len(words) != len(matrices):
        raise ValueError('one sampled matrix per physical step required')
    state, result = 0, []
    for raw, matrix in zip(words, matrices):
        result.append(raw ^ expansion(state))
        state = matrix_action(matrix, state) ^ feedback(raw)
    return result


def inner_transpose(words, matrices):
    if len(words) != len(matrices):
        raise ValueError('one sampled matrix per physical step required')
    ac, cc = local_maps()
    state, result = 0, [0] * len(words)
    for i in range(len(words) - 1, -1, -1):
        raw = words[i]
        result[i] = raw ^ apply_transpose(cc, state)
        state = matrix_adjoint(matrices[i], state) ^ apply_transpose(ac, raw)
    return result


def byte_route(groups, random, *, padded_stride=260):
    """Return stream-byte index -> padded outer-coordinate base address."""
    per_group = [random.sample(range(32), 32) for _ in range(groups)]
    route = [0] * (32 * groups)
    for region in range(32):
        slots = random.sample(range(groups), groups)
        for group in range(groups):
            route[region * groups + slots[group]] = (
                padded_stride * group + 8 * per_group[group][region])
    return route


def tower16_multiply(a, b):
    a0, a1, b0, b1 = a & 255, a >> 8, b & 255, b >> 8
    low = multiply(a0, b0) ^ multiply(0x20, multiply(a1, b1))
    high = multiply(a0, b1) ^ multiply(a1, b0) ^ multiply(a1, b1)
    return low | (high << 8)


@lru_cache(None)
def rs16_generator():
    def product(items):
        result = 1
        for item in items:
            result = multiply(result, item, degree=4, modulus=0x13)
        return result
    rows = []
    for point in range(16):
        row = []
        for source in range(8):
            numerator = product(point ^ other for other in range(8) if other != source)
            denominator = product(source ^ other for other in range(8) if other != source)
            inv = next(value for value in range(1, 16)
                       if multiply(value, denominator, degree=4, modulus=0x13) == 1)
            row.append(multiply(numerator, inv, degree=4, modulus=0x13))
        rows.append(tuple(row))
    return tuple(rows)


def outer_forward(message, scalars):
    """Native-layout RS16 outer: four rows, sixteen symbol adjoints.

    Four row nibbles form each 16-bit symbol. Absorbing the fixed 4x4
    transpose into the stored symbol matrix leaves its binary adjoint in
    this physical output order, as in the retained native-field outer.
    """
    generator = rs16_generator()
    result = 0
    for symbol, scalar in enumerate(scalars):
        coded = 0
        for lane in range(4):
            value = 0
            for source in range(8):
                nibble = (message >> (32 * lane + 4 * source)) & 15
                value ^= multiply(generator[symbol][source], nibble,
                                  degree=4, modulus=0x13)
            coded |= value << (4 * lane)
        mixed = sum(dot(tower16_multiply(scalar, 1 << bit), coded) << bit
                    for bit in range(16))
        result |= mixed << (16 * symbol)
    return result


def outer_columns(scalars):
    return tuple(outer_forward(1 << bit, scalars) for bit in range(128))


def code_forward(message, outer, route, matrices):
    groups = len(outer)
    routed = [0] * (260 * groups)
    for group, columns in enumerate(outer):
        encoded = apply_columns(columns, (message >> (128 * group)) & ((1 << 128) - 1))
        for bit in range(256):
            routed[260 * group + bit] = (encoded >> bit) & 1
    stream = [sum(routed[base + bit] << bit for bit in range(8)) for base in route]
    words = [sum(stream[i + j] << (8 * j) for j in range(8))
             for i in range(0, len(stream), 8)]
    encoded = inner_forward(words, matrices)
    return sum(word << (64 * step) for step, word in enumerate(encoded))


def code_transpose(word, outer, route, matrices):
    groups = len(outer)
    words = [(word >> (64 * step)) & ((1 << 64) - 1) for step in range(len(matrices))]
    reversed_words = inner_transpose(words, matrices)
    routed = [0] * (260 * groups)
    for i, base in enumerate(route):
        value = (reversed_words[i // 8] >> (8 * (i % 8))) & 255
        for bit in range(8):
            routed[base + bit] = (value >> bit) & 1
    result = 0
    for group, columns in enumerate(outer):
        value = sum(routed[260 * group + bit] << bit for bit in range(256))
        result |= apply_transpose(columns, value) << (128 * group)
    return result


def expansion_spectrum(*, points=tuple(range(8)), scales=None):
    return dict(sorted(Counter(expansion(state, points=points, scales=scales).bit_count()
                               for state in range(65536)).items()))


def local_summary(*, points=tuple(range(8)), scales=None):
    a_columns, c_columns = local_maps(points=points, scales=scales)
    return dict(points=list(points), scales=list(scales or (1,) * len(points)),
                a_rank=rank(a_columns), c_rank=rank(c_columns),
                ca_zero=all(apply_columns(c_columns, value) == 0 for value in a_columns),
                one_packet_ranks=[rank(c_columns[8 * i:8 * i + 8]) for i in range(8)],
                two_packet_ranks=sorted(set(rank(c_columns[8 * i:8 * i + 8] +
                                                c_columns[8 * j:8 * j + 8])
                                            for i, j in combinations(range(8), 2))),
                expansion_spectrum=expansion_spectrum(points=points, scales=scales))


if __name__ == '__main__':
    import json
    print(json.dumps(local_summary(), indent=2))
