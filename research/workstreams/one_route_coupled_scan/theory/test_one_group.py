"""Exact small-field and ordered-DP checks for the all-one-group proof."""

from decimal import ROUND_CEILING, localcontext
from fractions import Fraction as F
from itertools import combinations, product

from check_one_group import profile_moments
from check_sparse_full_gate import upward


def mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        if a & 4:
            a ^= 7
        b >>= 1
    return result


def kernel(matrix, columns):
    rows = [list(row) for row in matrix]
    pivots = []
    for column in range(columns):
        selected = next((i for i in range(len(pivots), len(rows)) if rows[i][column]), None)
        if selected is None:
            continue
        pivot = len(pivots)
        rows[pivot], rows[selected] = rows[selected], rows[pivot]
        inverse = (0, 1, 3, 2)[rows[pivot][column]]
        rows[pivot] = [mul(x, inverse) for x in rows[pivot]]
        for i in range(len(rows)):
            if i != pivot and rows[i][column]:
                scalar = rows[i][column]
                rows[i] = [a ^ mul(scalar, b) for a, b in zip(rows[i], rows[pivot])]
        pivots.append(column)
    basis = []
    for free in range(columns):
        if free not in pivots:
            row = [0] * columns
            row[free] = 1
            for i, pivot in enumerate(pivots):
                row[pivot] = rows[i][free]
            basis.append(row)
    return len(pivots), basis


def polynomial_product(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] ^= mul(x, y)
    return result


def line_constraints():
    lines = [(1, x) for x in range(4)] + [(0, 1)]
    fixtures = deficient = 0
    for b in range(1, 4):
        for r in range(5):
            evaluations = []
            for alpha in range(r):
                row = [1]
                for _ in range(1, b):
                    row.append(mul(row[-1], alpha))
                evaluations.append(row)
            for pattern in product(lines, repeat=r):
                constraints = [[mul(a, x) for x in point] + [mul(c, x) for x in point]
                               for (a, c), point in zip(pattern, evaluations)]
                rank, basis = kernel(constraints, 2 * b)
                if rank < r:
                    deficient += 1
                    # Every pair of row-pair kernel vectors has zero polynomial determinant.
                    for v, w in combinations(basis, 2):
                        assert polynomial_product(v[:b], w[b:]) == polynomial_product(v[b:], w[:b])
                fixtures += 1
    assert deficient > 0
    print(f"PASS: {fixtures} exact GF4 line systems; {deficient} deficient kernels have rational rank<=1.")


def pair_bridge():
    spaces = set()
    for generators in product(range(8), repeat=3):
        words = {0}
        for generator in generators:
            words |= {word ^ generator for word in words}
        spaces.add(tuple(sorted(words)))
    fixtures = 0
    for words in spaces:
        for support in range(1, 8):
            shortened = [word for word in words if word & ~support == 0]
            total_support = 0
            for word in shortened:
                total_support |= word
            if total_support != support:
                continue
            pairs = sum(a | b == support for a in shortened for b in shortened)
            dimension_power = len(shortened)**2
            assert pairs >= F(dimension_power * (4 - support.bit_count()), 4)
            assert F(dimension_power, 4) <= F(pairs, 4 - support.bit_count())
            fixtures += 1
    print(f"PASS: {fixtures} exact-pair support bridges across all{len(spaces)} binary length3 subspaces.")


def ordered_profiles():
    matrices = [[[F(1), F(0)], [F(0), F(2, 3)]],
                [[F(0), F(3, 7)], [F(2, 5), F(1, 11)]],
                [[F(1, 5), F(2, 7)], [F(3, 5), F(1, 13)]]]
    terminal = F(7, 5)
    fixtures = 0
    for columns in range(4):
        exact = [[F(0) for _ in range(columns + 1)] for _ in range(columns + 1)]
        for word in product(range(3), repeat=columns):
            row = [F(1), F(0)]
            for label in word:
                matrix = matrices[label]
                row = [sum((row[k] * matrix[k][j] for k in range(2)), F(0)) for j in range(2)]
            h, r = sum(bool(x) for x in word), word.count(1)
            exact[h][r] += row[0] + row[1] * terminal
        with localcontext() as context:
            context.prec = 80
            context.rounding = ROUND_CEILING
            converted = [[[upward(x) for x in row] for row in matrix] for matrix in matrices]
            observed = profile_moments(*converted, upward(terminal), columns=columns)
            for h in range(columns + 1):
                for r in range(h + 1):
                    difference = F(observed[h][r]) - exact[h][r]
                    assert 0 <= difference < F(1, 10**60)
                    fixtures += 1
    print(f"PASS: {fixtures} ordered rank-profile coefficients against exact direct enumeration.")


if __name__ == "__main__":
    line_constraints()
    pair_bridge()
    ordered_profiles()
