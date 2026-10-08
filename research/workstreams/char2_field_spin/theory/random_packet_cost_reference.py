"""Exact algebra and cost checks for fresh Toeplitz packet maps.

This is a scalar reference, not a kernel, benchmark, or distance certificate.
The primary map is Y=X+A R, Q=R+C X with fresh independent Toeplitz A,C.
An optional two-symbol feedthrough uses a uniform nonscalar invertible
trace-zero matrix; its action on every nonzero GF16 vector is uniform.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from random import Random

from grs16_12_parity_reference import multiply


MUL = tuple(tuple(multiply(a, b) for b in range(16)) for a in range(16))


def xor_sum(values):
    result = 0
    for value in values:
        result ^= value
    return result


def dot(left, right):
    return xor_sum(MUL[a][b] for a, b in zip(left, right))


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def transpose(matrix):
    return tuple(zip(*matrix))


def identity(size):
    return tuple(tuple(int(i == j) for j in range(size)) for i in range(size))


def toeplitz(rows, columns, taps):
    if len(taps) != rows + columns - 1:
        raise ValueError("one tap per Toeplitz diagonal required")
    return tuple(tuple(taps[i-j+columns-1] for j in range(columns))
                 for i in range(rows))


def toeplitz2(diagonal, upper, lower, values):
    """Three products; coefficient pre-additions belong to setup."""
    x0, x1 = values
    p = MUL[diagonal][x0 ^ x1]
    q = MUL[diagonal ^ upper][x1]
    r = MUL[diagonal ^ lower][x0]
    return p ^ q, p ^ r


def emission2(taps, state):
    """Even width t by two Toeplitz map: 3t/2 scalar products."""
    width = len(taps)-1
    if width <= 0 or width % 2:
        raise ValueError("positive even emission width required")
    result = []
    for row in range(0, width, 2):
        result.extend(toeplitz2(taps[row+1], taps[row], taps[row+2], state))
    return tuple(result)


def injection2(taps, values):
    """Two by even width t Toeplitz map: 3t/2 scalar products."""
    width = len(values)
    if width <= 0 or width % 2 or len(taps) != width+1:
        raise ValueError("positive even injection width and t+1 taps required")
    result = [0, 0]
    for column in range(0, width, 2):
        index = width-1-column
        block = toeplitz2(taps[index], taps[index-1], taps[index+1], values[column:column+2])
        result[0] ^= block[0]
        result[1] ^= block[1]
    return tuple(result)


def balanced_feedthroughs():
    """Every nonzero GF16^2 input has 15 preimages at each nonzero output."""
    return tuple(((a, b), (c, a))
                 for a in range(16) for b in range(16) for c in range(16)
                 if (b or c) and MUL[a][a] ^ MUL[b][c])


def matrix_inverse(matrix):
    size = len(matrix)
    work = [list(row)+list(unit) for row, unit in zip(matrix, identity(size))]
    for column in range(size):
        pivot = next((i for i in range(column, size) if work[i][column]), None)
        if pivot is None:
            raise ValueError("singular feedthrough")
        work[column], work[pivot] = work[pivot], work[column]
        reciprocal = next(a for a in range(1, 16) if MUL[a][work[column][column]] == 1)
        work[column] = [MUL[reciprocal][x] for x in work[column]]
        for row in range(size):
            if row != column:
                factor = work[row][column]
                work[row] = [x ^ MUL[factor][y] for x, y in zip(work[row], work[column])]
    return tuple(tuple(row[size:]) for row in work)


def coefficient_rank(matrix):
    work = [list(row) for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(rank, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        reciprocal = next(a for a in range(1, 16) if MUL[a][work[rank][column]] == 1)
        work[rank] = [MUL[reciprocal][x] for x in work[rank]]
        for row in range(rank+1, len(work)):
            factor = work[row][column]
            work[row] = [x ^ MUL[factor][y] for x, y in zip(work[row], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def toeplitz_tap_rank(rows, columns, vector):
    images = []
    for index in range(rows+columns-1):
        taps = tuple(int(j == index) for j in range(rows+columns-1))
        images.append(matvec(toeplitz(rows, columns, taps), vector))
    return coefficient_rank(transpose(images))


def forward(packets, maps, initial):
    state = tuple(initial)
    outputs = []
    for values, (feedthrough, emission, injection) in zip(packets, maps):
        output = tuple(a ^ b for a, b in zip(matvec(feedthrough, values), matvec(emission, state)))
        state = tuple(a ^ b for a, b in zip(state, matvec(injection, values)))
        outputs.append(output)
    return tuple(outputs), state


def causal_inverse(outputs, maps, initial):
    state = tuple(initial)
    packets = []
    for output, (feedthrough, emission, injection) in zip(outputs, maps):
        mixed = tuple(a ^ b for a, b in zip(output, matvec(emission, state)))
        values = matvec(matrix_inverse(feedthrough), mixed)
        state = tuple(a ^ b for a, b in zip(state, matvec(injection, values)))
        packets.append(values)
    return tuple(packets), state


def adjoint(dual_outputs, maps, terminal):
    state = tuple(terminal)
    inputs = []
    for output, (feedthrough, emission, injection) in zip(reversed(dual_outputs), reversed(maps)):
        inputs.append(tuple(a ^ b for a, b in zip(matvec(transpose(feedthrough), output),
                                                matvec(transpose(injection), state))))
        state = tuple(a ^ b for a, b in zip(state, matvec(transpose(emission), output)))
    return tuple(reversed(inputs)), state


def check():
    family = balanced_feedthroughs()
    assert len(family) == 3825
    # Exhaust every nonzero vector and every three-tap matrix.  This checks
    # the three-product circuit, the uniform Toeplitz image, and the balanced
    # invertible subfamily independently of a presumed group structure.
    for x0 in range(16):
        for x1 in range(16):
            if not (x0 or x1):
                continue
            full_histogram, balanced_histogram = Counter(), Counter()
            for a in range(16):
                for b in range(16):
                    for c in range(16):
                        expected = MUL[a][x0] ^ MUL[b][x1], MUL[c][x0] ^ MUL[a][x1]
                        assert toeplitz2(a, b, c, (x0, x1)) == expected
                        full_histogram[expected] += 1
                        if (b or c) and MUL[a][a] ^ MUL[b][c]:
                            balanced_histogram[expected] += 1
            assert len(full_histogram) == 256 and set(full_histogram.values()) == {16}
            assert len(balanced_histogram) == 255 and set(balanced_histogram.values()) == {15}
            assert (0, 0) not in balanced_histogram
            assert toeplitz_tap_rank(16, 2, (x0, x1)) == 16

    rng = Random(20261005)
    for width in (2, 4, 8, 16, 32):
        for _ in range(64):
            taps = tuple(rng.randrange(16) for _ in range(width+1))
            values = tuple(rng.randrange(16) for _ in range(width))
            state = tuple(rng.randrange(16) for _ in range(2))
            assert emission2(taps, state) == matvec(toeplitz(width, 2, taps), state)
            assert injection2(taps, values) == matvec(toeplitz(2, width, taps), values)
            if any(values):
                assert toeplitz_tap_rank(2, width, values) == 2
        for length in (1, 2, 7):
            for balanced in (False, True) if width == 2 else (False,):
                maps = []
                for _ in range(length):
                    feedthrough = rng.choice(family) if balanced else identity(width)
                    emission = toeplitz(width, 2, tuple(rng.randrange(16) for _ in range(width+1)))
                    injection = toeplitz(2, width, tuple(rng.randrange(16) for _ in range(width+1)))
                    maps.append((feedthrough, emission, injection))
                packets = tuple(tuple(rng.randrange(16) for _ in range(width)) for _ in range(length))
                initial, terminal = tuple(rng.randrange(16) for _ in range(2)), tuple(rng.randrange(16) for _ in range(2))
                dual = tuple(tuple(rng.randrange(16) for _ in range(width)) for _ in range(length))
                output, final = forward(packets, maps, initial)
                assert causal_inverse(output, maps, initial) == (packets, final)
                transposed, first_dual = adjoint(dual, maps, terminal)
                lhs = xor_sum(dot(x, y) for x, y in zip(output, dual)) ^ dot(final, terminal)
                rhs = xor_sum(dot(x, y) for x, y in zip(packets, transposed)) ^ dot(initial, first_dual)
                assert lhs == rhs

    # Leading GFNI512-equivalent counts, with four 128-bit records per ZMM.
    assert Fraction(48, 4*16) == Fraction(3, 4)
    assert Fraction(3, 4)*Fraction(48, 32)+Fraction(162, 4*32) == Fraction(153, 64)
    assert Fraction(3, 4)*Fraction(32, 24)+Fraction(81, 4*24) == Fraction(59, 32)
    assert Fraction(3, 4)*Fraction(16, 12)+Fraction(7, 4*12) == Fraction(55, 48)
    assert Fraction(9, 4*2)*Fraction(16, 12)+Fraction(7, 4*12) == Fraction(79, 48)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        check()
        print("GF16 uniform Toeplitz images, balanced feedthrough, 3-product circuits, causal inverse, adjoint, and normalized costs pass")
    else:
        print("B=I, s=2, even width t: 3t scalar products = 3/4 GFNI512-equivalent per encoded record")
        print("Optional balanced t=2 feedthrough: total nine scalar products per packet")
