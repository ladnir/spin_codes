"""Exact AES-field rank-seven Toeplitz4 circuit and counted alternatives.

This is an isolated scalar algebra helper, not a native kernel or timing
claim. Seven projective evaluation points lie in GF16, not a nonexistent
GF8 subfield of GF256. The original-map circuit and fixed-congruence
circuit are explicitly distinguished. The optional GF4 transform uses
eight variable products by leaving one quadratic place un-diagonalized.
"""
from functools import lru_cache


OMEGA, OMEGA2, GAMMA = 0xbc, 0xbd, 0x5c
POINTS = (0, 1, OMEGA, OMEGA2, None, GAMMA, GAMMA ^ 1)
COUNTS = {
    'original_toeplitz_rank7': dict(variable_products=7, fixed_maps=8, data_xors=25,
                                  prepared_bytes=7, variable_key_expansions=7),
    'congruent_rank7': dict(variable_products=7, fixed_maps=6, data_xors=19,
                           prepared_bytes=7, variable_key_expansions=7),
    'congruent_quadratic_rank8': dict(variable_products=8, fixed_maps=4, data_xors=18,
                                     prepared_bytes=8, variable_key_expansions=8),
}


def mul(a, b):
    result = 0
    for _ in range(8):
        if b & 1:
            result ^= a
        a = ((a << 1) ^ (0x11b if a & 128 else 0)) & 255
        b >>= 1
    return result


def power(a, n):
    result = 1
    while n:
        if n & 1:
            result = mul(result, a)
        a = mul(a, a)
        n >>= 1
    return result


def dot(a, b):
    result = 0
    for x, y in zip(a, b):
        result ^= mul(x, y)
    return result


def apply(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def transpose(matrix):
    return tuple(zip(*matrix))


def inverse(matrix):
    n = len(matrix)
    rows = [list(row)+[int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            raise ValueError('Singular matrix')
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = power(rows[j][j], 254)
        rows[j] = [mul(x, scale) for x in rows[j]]
        for i in range(n):
            if i != j and rows[i][j]:
                scale = rows[i][j]
                rows[i] = [x ^ mul(scale, y) for x, y in zip(rows[i], rows[j])]
    return tuple(tuple(row[n:]) for row in rows)


def evaluation_row(point, degree):
    # Projective infinity selects the leading coefficient.
    return (0,)*(degree-1)+(1,) if point is None else tuple(power(point, j) for j in range(degree))


@lru_cache(None)
def key_transform():
    return inverse(transpose(tuple(evaluation_row(point, 7) for point in POINTS)))


def prepare7(t):
    """Setup only: solve V^T lambda=t, after sampling/rejecting original T."""
    if len(t) != 7:
        raise ValueError('Seven Toeplitz/Hankel diagonals required')
    return apply(key_transform(), t)


def monomial_to_novel(c):
    u2 = c[2] ^ c[3]
    u1 = c[1] ^ u2
    return c[0], u1, u2, c[3]


def novel_transpose_to_monomial(u):
    c2 = u[1] ^ u[2]
    c3 = c2 ^ u[3]
    return u[0], u[1], c2, c3


def eval7_novel(u):
    """4 fixed maps, 9 XORs; basis(1,x,x^2+x,x^3+x^2)."""
    low = u[1] ^ u[3]
    base = u[0] ^ u[2]
    y1 = u[0] ^ u[1]
    y2 = base ^ mul(OMEGA, low)
    y3 = y2 ^ low
    d = u[1] ^ mul(OMEGA, u[3])
    c = u[0] ^ mul(OMEGA, u[2])
    y5 = c ^ mul(GAMMA, d)
    y6 = y5 ^ d
    return u[0], y1, y2, y3, u[3], y5, y6


def eval7_novel_transpose(r):
    """4 fixed maps, 12 XORs; algebraic transpose, not inverse."""
    t0 = r[5] ^ r[6]
    t1 = r[6] ^ mul(GAMMA, t0)
    t2 = r[2] ^ r[3]
    t3 = r[3] ^ mul(OMEGA, t2)
    u0 = r[0] ^ r[1] ^ t2 ^ t0
    u1 = r[1] ^ t3 ^ t1
    u2 = t2 ^ mul(OMEGA, t0)
    u3 = r[4] ^ t3 ^ mul(OMEGA, t1)
    return u0, u1, u2, u3


def original_from_keys(keys, x):
    """Original T(t)x exactly: input reversal is register naming only."""
    values = eval7_novel(monomial_to_novel(x[::-1]))
    products = tuple(mul(key, value) for key, value in zip(keys, values))
    return novel_transpose_to_monomial(eval7_novel_transpose(products))


@lru_cache(None)
def congruence_matrix():
    # C converts evaluations on GF4 to monomial coefficients.
    return inverse(tuple(evaluation_row(x, 4) for x in POINTS[:4]))


def eval7_congruent(x):
    """3 fixed maps, 8 XORs; E*C, first four output rows are identity."""
    a = x[0] ^ x[1]
    b = x[2] ^ x[3]
    e = a ^ b
    d = mul(OMEGA, e) ^ a
    u = x[0] ^ x[3]
    c = mul(OMEGA2, u) ^ x[2]
    v = mul(GAMMA, d) ^ c
    w = v ^ d
    return x[0], x[1], x[2], x[3], e, v, w


def eval7_congruent_transpose(r):
    """3 fixed maps, 11 XORs."""
    h = r[5] ^ r[6]
    d = r[6] ^ mul(GAMMA, h)
    e = r[4] ^ mul(OMEGA, d)
    a = e ^ d
    u = mul(OMEGA2, h)
    x0 = r[0] ^ a ^ u
    x1 = r[1] ^ a
    x2 = r[2] ^ e ^ h
    x3 = r[3] ^ e ^ u
    return x0, x1, x2, x3


def congruent_from_keys(keys, x):
    """C^T H(t) C x, not the old T(t)x map."""
    values = eval7_congruent(x)
    products = tuple(mul(key, value) for key, value in zip(keys, values))
    return eval7_congruent_transpose(products)


def quadratic_eval(x):
    """2 GF4 fixed maps, 6 XORs; last outputs are the quadratic-place c,d."""
    a = x[0] ^ x[1]
    b = x[2] ^ x[3]
    e = a ^ b
    d = mul(OMEGA, e) ^ a
    u = x[0] ^ x[3]
    c = mul(OMEGA2, u) ^ x[2]
    return x[0], x[1], x[2], x[3], e, c, d


def quadratic_eval_transpose(r):
    """2 GF4 fixed maps, 9 XORs."""
    e = r[4] ^ mul(OMEGA, r[6])
    a = e ^ r[6]
    u = mul(OMEGA2, r[5])
    x0 = r[0] ^ a ^ u
    x1 = r[1] ^ a
    x2 = r[2] ^ e ^ r[5]
    x3 = r[3] ^ e ^ u
    return x0, x1, x2, x3


def prepare8(t):
    keys = prepare7(t)
    k = keys[5] ^ keys[6]
    ell = mul(GAMMA, keys[5]) ^ mul(GAMMA ^ 1, keys[6])
    return keys[:5]+(ell, k ^ ell, mul(OMEGA, k))


def quadratic_from_keys(keys, x):
    values = quadratic_eval(x)
    first = tuple(mul(keys[i], values[i]) for i in range(5))
    c, d = values[5:]
    shared = mul(keys[5], c ^ d)
    out_c = shared ^ mul(keys[6], c)
    out_d = shared ^ mul(keys[7], d)
    return quadratic_eval_transpose(first+(out_c, out_d))


def dense_hankel(t, x):
    return tuple(dot(t[i:i+4], x) for i in range(4))


def dense_original(t, x):
    return dense_hankel(t, x[::-1])


def dense_congruent(t, x):
    c = congruence_matrix()
    return apply(transpose(c), dense_hankel(t, apply(c, x)))
