"""Exact small-field checks for the source-rank gate; not a certificate.

The outer is four parallel RS[4,2] rows over GF4, extended to GF16.
Each column has two independent invertible two-record mixers. Constraints
may have coefficients outside GF4. No numerical job writes files.
Run with Python -B under a 30-second / 256-MiB process cap.
"""
from fractions import Fraction
from itertools import combinations, product

from test_small_state_bridge import E, FIELD, mul


def power(a, n):
    out = 1
    while n:
        if n & 1:
            out = mul(out, a)
        a = mul(a, a)
        n >>= 1
    return out


def rank(rows, width):
    rows = [list(row) for row in rows]
    assert all(len(row) == width for row in rows)
    pivot = 0
    for j in range(width):
        hit = next((i for i in range(pivot, len(rows)) if rows[i][j]), None)
        if hit is None:
            continue
        rows[pivot], rows[hit] = rows[hit], rows[pivot]
        inverse = power(rows[pivot][j], 14)
        rows[pivot] = [mul(inverse, a) for a in rows[pivot]]
        for i in range(len(rows)):
            if i != pivot and rows[i][j]:
                scale = rows[i][j]
                rows[i] = [a ^ mul(scale, b) for a, b in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def matmul(left, right):
    assert right and all(len(row) == len(right[0]) for row in right)
    assert all(len(row) == len(right) for row in left)
    out = []
    for row in left:
        target = []
        for j in range(len(right[0])):
            total = 0
            for a, other in zip(row, right):
                total ^= mul(a, other[j])
            target.append(total)
        out.append(tuple(target))
    return tuple(out)


def outer_column(point, mixer):
    evaluate = [tuple((1 if j == 2*i else point if j == 2*i+1 else 0)
                      for j in range(8)) for i in range(4)]
    return matmul(mixer, evaluate)


def source_rank_checks():
    identity = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
    mixed = ((1, 1, 0, 0), (1, 6, 0, 0), (0, 0, 1, 6), (0, 0, 1, 7))
    # This invertible row transform uses extension-field entries. Its
    # prefixes give arbitrary-rank constraint fixtures, not only erasures.
    transform = ((1, 2, 3, 4), (0, 1, 5, 8), (0, 0, 1, 9), (0, 0, 0, 1))
    assert rank(mixed, 4) == rank(transform, 4) == 4
    fixtures = 0
    full_rank_cutoffs = 0
    for style in range(3):
        columns = [outer_column(point, mixed if (j+style) % 2 else identity)
                   for j, point in enumerate(E)]
        constraints = []
        for j in range(4):
            # Rotate independent rows, retaining their prefix ranks.
            basis = transform[j:] + transform[:j] if style else identity
            constraints.append([matmul(basis[:r], columns[j]) for r in range(5)])
        for profile in product(range(5), repeat=4):
            equations = [row for j, r in enumerate(profile) for row in constraints[j][r]]
            actual = rank(equations, 8)
            lower = sum(sorted(profile)[-2:])
            assert actual >= lower
            for chosen in combinations(range(4), 2):
                selected = [row for j in chosen for row in constraints[j][profile[j]]]
                assert rank(selected, 8) == sum(profile[j] for j in chosen)
            if lower == 8:
                assert actual == 8
                full_rank_cutoffs += 1
            fixtures += 1
    return fixtures, full_rank_cutoffs


def projective(q, c):
    return sum(q**i for i in range(c))


def diagonal_rank_checks():
    # A 2x2 superregular subfield map stays superregular after extension.
    h = ((1, 1), (1, 6))
    fixtures = assignments = 0
    for field in (E, FIELD):
        q = len(field)
        vectors = (((1,), (0,)), ((0,), (1,)), ((1,), (1,)),
                   ((1,), (field[2],)), ((1, 0), (0, 1)))
        for image in vectors:
            c = len(image[0])
            assert rank(image, c) == c
            for chosen in ((0,), (1,), (0, 1)):
                checks = len(chosen)
                hist = [0]*(c+1)
                for diagonal in product(field[1:], repeat=2):
                    scaled = tuple(tuple(mul(diagonal[i], v) for v in image[i])
                                   for i in range(2))
                    transformed = matmul(tuple(h[j] for j in chosen), scaled)
                    hist[rank(transformed, c)] += 1
                    assignments += 1
                for d in range(c):
                    probability = Fraction(sum(hist[:d+1]), (q-1)**2)
                    cap = min(Fraction(1), Fraction(projective(q, c),
                              (q-1)**checks * projective(q, c-d)))
                    assert probability <= cap
                    fixtures += 1
    return fixtures, assignments


def top_half_identity():
    fixtures = 0
    for values in product(range(5), repeat=4):
        top = sum(sorted(values)[-2:])
        assert top == min(2*v+sum(max(x-v, 0) for x in values) for v in range(5))
        for dimension in range(1, 9):
            retained = any(2*v+sum(max(x-v, 0) for x in values) < dimension
                           for v in range(5))
            assert retained == (top < dimension)
        fixtures += 1
    return fixtures


def projected_rank_checks():
    # Eliminate one incoming-state coordinate from a Cauchy check matrix.
    # Its column span is unchanged by that coordinate's nonzero D1 scalar.
    h = tuple(tuple(power(i ^ (j+4), 14) for j in range(4)) for i in range(4))
    fixtures = assignments = minors = 0
    images = (((1,), (2,), (3,)), ((1, 0), (0, 1), (1, 6)))
    for incoming, row_count in product((0, 1), (2, 3)):
        outside = tuple(j for j in range(4) if j != incoming)
        projected = []
        for i in range(1, row_count):
            scale = mul(h[i][incoming], power(h[0][incoming], 14))
            assert h[i][incoming] ^ mul(scale, h[0][incoming]) == 0
            projected.append(tuple(h[i][j] ^ mul(scale, h[0][j]) for j in outside))
        checks = row_count-1
        for size in range(1, checks+1):
            for chosen in combinations(range(3), size):
                assert rank([[row[j] for j in chosen] for row in projected], size) == size
                minors += 1
        for image in images:
            k = len(image[0])
            assert rank(image, k) == k
            hist = [0]*(k+1)
            for diagonal in product(FIELD[1:], repeat=3):
                scaled = tuple(tuple(mul(diagonal[i], v) for v in image[i])
                               for i in range(3))
                hist[rank(matmul(projected, scaled), k)] += 1
                assignments += 1
            for d in range(k):
                probability = Fraction(sum(hist[:d+1]), 15**3)
                cap = min(Fraction(1), Fraction(projective(16, k),
                          15**checks * projective(16, k-d)))
                assert probability <= cap
                fixtures += 1
    return fixtures, assignments, minors


def main():
    source, cutoffs = source_rank_checks()
    local, assignments = diagonal_rank_checks()
    thresholds = top_half_identity()
    projected, projected_assignments, minors = projected_rank_checks()
    print(dict(exact=True, whole_code_certificate=False,
               source_rank_fixtures=source, full_rank_source_cutoffs=cutoffs,
               diagonal_rank_caps=local, diagonal_assignments=assignments,
               top_half_threshold_fixtures=thresholds,
               projected_rank_caps=projected,
               projected_diagonal_assignments=projected_assignments,
               projected_mds_minors=minors))


if __name__ == "__main__":
    main()
