"""Tiny exact checks of the partial-column proof interfaces, not K20 evidence."""

from decimal import Decimal as D, ROUND_CEILING, localcontext
from fractions import Fraction as F
from itertools import combinations

from check_partial_columns import coefficient_moment
from check_sparse_full_gate import upward


def rank(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = row
                break
            row ^= basis[pivot]
    return len(basis)


def image(rows, vector):
    return sum(((row & vector).bit_count() & 1) << j for j, row in enumerate(rows))


def toeplitz_checks():
    matrices = [tuple(sum(((diagonals >> (i - j + 3)) & 1) << j for j in range(4))
                      for i in range(4)) for diagonals in range(128)]
    invertible = [matrix for matrix in matrices if rank(matrix) == 4]
    assert len(invertible) == 64  # (q-1) q^6 at q=2.
    subspaces = {(a,): (a,) for a in range(1, 16)}
    for a, b in combinations(range(1, 16), 2):
        subspaces.setdefault(tuple(sorted((a, b, a ^ b))), (a, b))
    assert len(subspaces) == 15 + 35
    fixtures = 0
    for basis in subspaces.values():
        t = len(basis)
        for allowed in range(1, 16):
            k = allowed.bit_count()
            forbidden = 15 ^ allowed
            equations = [sum(((vector >> j) & 1) << (i - j + 3) for j in range(4))
                         for i in range(4) if forbidden & (1 << i) for vector in basis]
            if forbidden:
                assert rank(equations) >= (4 - k) + t - 1
            count = sum(all(image(matrix, vector) & forbidden == 0 for vector in basis)
                        for matrix in invertible)
            if t == 1:
                bound = min(F(1), F(2**k - 1, 8))
            else:
                bound = (F(0), F(0), F(1, 4), F(1, 2), F(1))[k]
            assert F(count, len(invertible)) <= bound
            if t > k:
                assert count == 0
            fixtures += 1
    print(f"PASS: {fixtures} exact GF2 Toeplitz subspace fixtures; all64 invertible matrices.")


def ordered_dp_checks():
    empty = [[F(1), F(0)], [F(0), F(2, 3)]]
    occupied = [[F(0), F(3, 7)], [F(2, 5), F(1, 11)]]
    terminal = F(5, 3)
    fixtures = 0
    for columns in range(5):
        for degree in range(columns + 1):
            exact = F(0)
            for selected in combinations(range(columns), degree):
                row = [F(1), F(0)]
                for j in range(columns):
                    matrix = occupied if j in selected else empty
                    row = [sum((row[k] * matrix[k][i] for k in range(2)), F(0)) for i in range(2)]
                exact += row[0] + row[1] * terminal
            with localcontext() as context:
                context.prec = 80
                context.rounding = ROUND_CEILING
                observed = coefficient_moment([[upward(x) for x in row] for row in empty],
                                              [[upward(x) for x in row] for row in occupied],
                                              upward(terminal), degree, columns)
                assert F(observed) >= exact
                assert F(observed) - exact < F(1, 10**60)
            fixtures += 1
    print(f"PASS: {fixtures} ordered coefficient-DP fixtures, exact reference and upward arithmetic.")


def main():
    toeplitz_checks()
    ordered_dp_checks()


if __name__ == "__main__":
    main()
