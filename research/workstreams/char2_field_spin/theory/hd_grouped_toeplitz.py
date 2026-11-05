"""Exact interface and toy checks for invertible Toeplitz column mixers.

This changes only the grouped outer, not the recursive inner. See the
Toeplitz section of HD_GROUPED_OUTER.md for the conditional atom argument.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product


def column_cap(q, g):
    if q < 2 or g < 1:
        raise ValueError('Invalid field or column size')
    return Fraction(1, (q-1)*q**(g-1))


def group_cap(q=256, g=4, length=64, dimension=32):
    if not 0 < dimension < length:
        raise ValueError('Invalid RS dimensions')
    return column_cap(q, g)**(length-dimension)


def binary_field(q):
    modulus = {2: 3, 4: 7, 256: 0x11b}[q]
    def mul(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a & q:
                a ^= modulus
        return r
    return [[mul(a, b) for b in range(q)] for a in range(q)]


def matrix(t):
    g = (len(t)+1)//2
    if len(t) != 2*g-1:
        raise ValueError('Odd coefficient count required')
    return [[t[g-1+i-j] for j in range(g)] for i in range(g)]


def nonsingular(t, mul):
    """Independent normalized Gaussian elimination, unlike native setup."""
    a = matrix(t)
    g = len(a)
    for c in range(g):
        p = next((i for i in range(c, g) if a[i][c]), None)
        if p is None:
            return False
        a[c], a[p] = a[p], a[c]
        inv = mul[a[c][c]].index(1)
        a[c] = [mul[x][inv] for x in a[c]]
        for i in range(c+1, g):
            f = a[i][c]
            a[i] = [x ^ mul[f][y] for x, y in zip(a[i], a[c])]
    return True


def circuit(t, x, mul):
    """Three/nine-product recursive circuit; coefficients may be zero."""
    def pair(a, b, c, u, v):
        z = mul[a][u ^ v]
        return z ^ mul[a ^ b][v], z ^ mul[a ^ c][u]
    if len(x) == 2 and len(t) == 3:
        return pair(t[1], t[0], t[2], *x)
    if len(x) != 4 or len(t) != 7:
        raise ValueError('Only the two/four-symbol circuits are implemented')
    u = pair(t[3], t[2], t[4], x[0] ^ x[2], x[1] ^ x[3])
    v = pair(t[5] ^ t[3], t[4] ^ t[2], t[6] ^ t[4], *x[:2])
    w = pair(t[1] ^ t[3], t[0] ^ t[2], t[2] ^ t[4], *x[2:])
    return u[0] ^ w[0], u[1] ^ w[1], u[0] ^ v[0], u[1] ^ v[1]


@lru_cache(None)
def toy_counts(q=4, g=4):
    """Enumerate every accepted matrix and every nonzero input/output atom."""
    if (q, g) not in ((2, 2), (2, 3), (2, 4), (4, 2), (4, 3), (4, 4)):
        raise ValueError('Exhaustive toy sizes only')
    mul = binary_field(q)
    bits = q.bit_length()-1
    words = tuple(tuple((x >> (j*bits)) & (q-1) for j in range(g)) for x in range(q**g))
    counts = [[0]*q**g for _ in words]
    accepted = 0
    for t in product(range(q), repeat=2*g-1):
        if not nonsingular(t, mul):
            continue
        accepted += 1
        m = matrix(t)
        basis = [[sum(mul[m[i][j]][a] << (i*bits) for i in range(g))
                  for a in range(q)] for j in range(g)]
        for row, x in zip(counts[1:], words[1:]):
            y = 0
            for j, a in enumerate(x):
                y ^= basis[j][a]
            row[y] += 1
    counts[0][0] = accepted
    return counts, accepted
