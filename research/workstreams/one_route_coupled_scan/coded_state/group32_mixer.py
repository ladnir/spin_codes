"""Exact structural checks for cheap intergroup mixing. NOT a certificate.

Use a sparse companion over GF256, G=8 or16. Its action is a cyclic
buffer rotation, one group scaling, and a few group XORs. Uniform
full-field feedback and the outer/route are unchanged. No native timing.
"""
import time
from check_maps import Field
from cheap_reset_maps import BinaryField, novel_values, xor


def trim(p):
    while p and not p[-1]:
        p.pop()
    return p


def remainder(f, p, q):
    p = p[:]
    inv = f.inverse[q[-1]]
    while len(p) >= len(q):
        scale, shift = f.table[p[-1]][inv], len(p)-len(q)
        for i, value in enumerate(q):
            p[i+shift] ^= f.table[scale][value]
        trim(p)
    return p


def square_mod(f, p, modulus):
    out = [0]*(2*len(p)-1)
    for i, value in enumerate(p):
        out[2*i] = f.table[value][value]
    return remainder(f, out, modulus)


def gcd(f, p, q):
    while q:
        p, q = q, remainder(f, p, q)
    return p


def irreducible(f, modulus):
    n = len(modulus)-1
    assert n & (n-1) == 0
    p = [0, 1]
    for i in range(1, n+1):
        for _ in range(8):
            p = square_mod(f, p, modulus)
        if i == n//2:
            difference = p+[0]*max(0, 2-len(p))
            difference[1] ^= 1
            if len(gcd(f, modulus, trim(difference))) != 1:
                return False
    return p == [0, 1]


def companion(f, vector, a, c):
    out = [f.table[a][vector[-1]]] + vector[:-1]
    for index in c:
        out[index] ^= vector[-1]
    return out


def find_companion(f, G, deadline):
    for c in ((1, 3, 4), (1, 2, 5), (1, 3, 5), (1, 2, 3)):
        for a in range(1, 256):
            modulus = [a]+[0]*(G-1)+[1]
            for index in c:
                modulus[index] = 1
            if irreducible(f, modulus):
                return a, c
            assert time.monotonic() < deadline, 'Bounded polynomial search expired'
    raise ValueError('No irreducible sparse polynomial in this family')


def group_columns(f, G, a, c, steps, initial=0):
    vector = [int(j == initial) for j in range(G)]
    columns = []
    for _ in range(steps):
        columns.append(vector)
        vector = companion(f, vector, a, c)
    return columns


def tensor_columns(f, G, a, c, d, times, initial=0):
    # All moment points are distinct in the toy field. Positive powers of
    # P are equivalent to the inverse-companion candidate, with same cost.
    field = BinaryField(8, 0x11b)
    vectors = group_columns(f, G, a, c, max(times)+1, initial)
    return [[f.table[u][w] for u in vectors[t]
             for w in novel_values(field, d.bit_length()-1, t*d+j)]
            for t in times for j in range(d)]


def column_rank(f, columns):
    return f.rank([list(row) for row in zip(*columns)])


def main():
    start = time.monotonic()
    deadline = start+40
    f = Field(8, 0x11b)
    for G in (8, 16):
        a, c = find_companion(f, G, deadline)
        vectors = group_columns(f, G, a, c, 8193)
        assert len(set(map(tuple, vectors))) == 8193
        # Irreducibility implies every nonzero GF256 vector is cyclic.
        # These rank identities extend to F, but arbitrary F-valued vectors
        # need not be cyclic there. Test every physical group and offsets.
        checks = 0
        for initial in range(G):
            columns = group_columns(f, G, a, c, 4*G, initial)
            for offset in (0, 1, G-1, 2*G):
                assert column_rank(f, columns[offset:offset+G]) == G
                checks += 1
        # Consecutive one-group full batches now contribute d*G reset rank,
        # not d. This is a local theorem, not an all-support MDS property.
        for d in (2, 4):
            columns = tensor_columns(f, G, a, c, d, list(range(G)))
            assert column_rank(f, columns) == d*G
        # Compare the old minimal one-group support across d+1 batches.
        d = 4
        columns = tensor_columns(f, G, a, c, d, list(range(d+1)))
        sparse = columns[::d]
        assert column_rank(f, sparse) == d+1  # old grouped rank was d
        # Negative control: P e_0=e_1. Tracking that group shift over two
        # batches leaves the group direction unchanged in the unrolled
        # check. Any d+1 distinct local moment columns then have rank d.
        # Hence this sparse mixer does NOT remove cheap short excursions.
        tracking = (tensor_columns(f, G, a, c, d, [0], initial=1)
                    + tensor_columns(f, G, a, c, d, [1], initial=0)[:1])
        assert column_rank(f, tracking) == d
        assert all(column_rank(f, tracking[:j]+tracking[j+1:]) == d
                   for j in range(d+1))
        assert time.monotonic() < deadline
        print(dict(kind='EXACT_LOCAL_NOT_CERTIFICATE', groups=G, field_modulus='0x11b',
                   polynomial=f'x^{G}+'+'+'.join(f'x^{i}' for i in c)+f'+{a}', period_greater_than=8192,
                   cyclic_rank_checks=checks, one_group_full_batch_rank=4*G,
                   old_weight5_pattern_rank=column_rank(f, sparse),
                   tracking_weight5_rank=column_rank(f, tracking),
                   conclusion='short-excursion weakness remains',
                   scalar_products_at_N21=2**21//G,
                   extra_xors_at_N21=len(c)*2**21//G), flush=True)
    print(dict(seconds=round(time.monotonic()-start, 3)), flush=True)


if __name__ == '__main__':
    main()
