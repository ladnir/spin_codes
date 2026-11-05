"""Exact ranks for the implemented RS[16,4] emission basis.

Run with python -B. No files are written. The tables count row masks, not
probabilities: using them in a proof still requires the actual routing and
conditional support distribution. This is not a distance certificate.
"""

from fractions import Fraction
from itertools import combinations
from math import comb
import json


def mul(a, b):
    """AES GF256 multiplication, independent of encoder lookup tables."""
    value = 0
    for _ in range(8):
        if b & 1:
            value ^= a
        a <<= 1
        if a & 256:
            a ^= 0x11B
        b >>= 1
    return value


def power(a, exponent):
    value = 1
    while exponent:
        if exponent & 1:
            value = mul(value, a)
        a = mul(a, a)
        exponent >>= 1
    return value


INVERSE = [0] + [power(x, 254) for x in range(1, 256)]


def basis(degree, point):
    # Scalar reproduction of hd_grouped_forward.cpp::basis, degrees 0..3.
    result = 1
    for bit in range(2):
        if degree & (1 << bit):
            value = normalizer = 1
            for j in range(1 << bit):
                value = mul(value, point ^ j)
                normalizer = mul(normalizer, (1 << bit) ^ j)
            result = mul(result, mul(value, INVERSE[normalizer]))
    return result


A = [[basis(degree, point) for degree in range(4)] for point in range(16)]


def rank(rows):
    """Row rank over AES GF256; empty matrices have rank zero."""
    rows = [list(row) for row in rows]
    if not rows:
        return 0
    pivot = 0
    for col in range(len(rows[0])):
        selected = next((i for i in range(pivot, len(rows)) if rows[i][col]), None)
        if selected is None:
            continue
        rows[pivot], rows[selected] = rows[selected], rows[pivot]
        inv = INVERSE[rows[pivot][col]]
        rows[pivot] = [mul(x, inv) for x in rows[pivot]]
        for i in range(pivot + 1, len(rows)):
            factor = rows[i][col]
            if factor:
                rows[i] = [x ^ mul(factor, y) for x, y in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def rank_counts():
    """counts[H][z][r] counts z-row subsets Z with rank A[Z,H] = r.

    All four-row minors are checked independently below. Therefore every
    >=4-row restriction has full column rank, including after deleting
    columns. Only the <=3-row restrictions need explicit enumeration.
    """
    counts = [[[0] * 5 for _ in range(17)] for _ in range(16)]
    for mask in range(16):
        columns = [j for j in range(4) if mask >> j & 1]
        for size in range(4):
            for chosen in combinations(range(16), size):
                r = rank([[A[i][j] for j in columns] for i in chosen])
                counts[mask][size][r] += 1
        for size in range(4, 17):
            counts[mask][size][len(columns)] = comb(16, size)
    return counts


def support_sum(counts, mask, tilt, tau):
    """Sum_Z tilt^(16-|Z|) tau^(-rank A[Z,H]), not an average."""
    return sum((number * tilt ** (16 - size) / tau**r
                for size, row in enumerate(counts[mask])
                for r, number in enumerate(row)), Fraction(0))


def self_check(counts):
    assert all(mul(x, INVERSE[x]) == 1 for x in range(1, 256))
    assert A[0] == [1, 0, 0, 0]
    for point, row in enumerate(A):
        quadratic = mul(mul(point, point ^ 1), INVERSE[6])
        assert row == [1, point, quadratic, mul(point, quadratic)]
    assert all(rank([A[i] for i in subset]) == 4
               for subset in combinations(range(16), 4))
    for mask in range(16):
        h = mask.bit_count()
        for size, row in enumerate(counts[mask]):
            assert sum(row) == comb(16, size)
            lower = max(0, h - max(0, 4 - size))
            upper = min(size, h)
            assert all(number == 0 or lower <= r <= upper
                       for r, number in enumerate(row))
        for tilt in (Fraction(1, 4), Fraction(1, 2), Fraction(1)):
            assert support_sum(counts, mask, tilt, Fraction(1)) == (1 + tilt)**16
    # In this basis, a one-dimensional state need not affect the first row.
    assert counts[2][1][0] == 1
    assert counts[4][1][0] == 2
    assert counts[8][1][0] == 2
    # More than one mask type occurs even at a fixed state-support size.
    assert counts[1] != counts[2]


def main():
    counts = rank_counts()
    self_check(counts)
    defects = []
    for mask in range(1, 16):
        h = mask.bit_count()
        rows = []
        for size in range(1, 4):
            bad = {str(r): n for r, n in enumerate(counts[mask][size])
                   if n and r < min(size, h)}
            if bad:
                rows.append({"zero_rows": size, "counts_by_rank": bad})
        if rows:
            defects.append({"state_mask": mask, "state_weight": h, "defects": rows})
    print(json.dumps({
        "status": "exact_local_rank_counts_not_a_whole_code_certificate",
        "field_polynomial": "0x11b", "row_points": list(range(16)),
        "basis": ["1", "x", "x(x+1)/6", "x^2(x+1)/6"],
        "checked_four_row_minors": comb(16, 4),
        "represented_submatrices": 16 * 2**16,
        "cardinality_only_min_formula_counterexamples": defects,
        "self_checks": "passed",
    }, indent=2))


if __name__ == "__main__":
    main()
