"""Exact checks for sparse uniform and independently shifted envelopes."""
import unittest
from fractions import Fraction as F
from itertools import combinations, product
from math import comb

from flint import arb, ctx

from hd_component_rank import arb_fraction, exact_kernels, gap_moment_exact
from hd_component_sparse import class_matrix, minimum_kernel
from hd_component_shift import birth_return_sum, birth_sum, shift_matrix


def region(cuts, rows, q, kernels):
    cuts = sorted(cuts)
    lengths = [b-a for a, b in zip(cuts, cuts[1:])]
    result = [[F(0) for _ in range(4)] for _ in range(4)]
    for old, new in product((0, 1), repeat=2):
        for inside in product((0, 1), repeat=len(lengths)):
            h = old*cuts[0]+new*(rows-cuts[-1])+sum(x*n for x, n in zip(inside, lengths))
            weight = (q-1)**(sum(inside)+new)
            for i in range(2):
                for j in range(2):
                    result[2*old+i][2*new+j] += weight*kernels[h][i][j]
    return result


class SparseTests(unittest.TestCase):
    def setUp(self):
        self.previous = ctx.prec
        ctx.prec = 160

    def tearDown(self):
        ctx.prec = self.previous

    def test_birth_hockey_stick(self):
        for rows in range(2, 16):
            for gap in range(rows):
                actual = sum(F(comb(rows-gap-1, h-1), comb(rows, h))
                             for h in range(1, rows-gap+1))
                self.assertEqual(actual, F(rows+1, (gap+1)*(gap+2)))

    def test_first_return_vandermonde(self):
        for rows in range(3, 16):
            for gap in range(rows-2):
                actual = sum(F((h-2)*(gap+1)*comb(rows-gap-2, h-2), comb(rows, h))
                             for h in range(3, rows-gap+1))
                self.assertEqual(actual, F(6*(rows+1)*(rows-gap-2),
                                          (gap+2)*(gap+3)*(gap+4)))

    def test_birth_sums(self):
        for rows, q, mark in product(range(3, 10), (4, 16), (F(1, 2), F(3, 4))):
            base, p = mark+(1-mark)/F(q-1), F(1, q-1)
            expected = sum(mark*base**(rows-1)*gap_moment_exact(rows, h, 1, base)
                           for h in range(1, rows+1))
            returned = sum((h-2)*mark**3*p*p*base**(rows-3)
                           *gap_moment_exact(rows, h, 2, base) for h in range(3, rows+1))
            self.assertTrue(birth_sum(rows, q, mark).overlaps(arb_fraction(expected)))
            self.assertTrue(birth_return_sum(rows, q, mark).overlaps(arb_fraction(returned)))

    def test_minimum_kernels(self):
        for rows, q, mark in product(range(3, 10), (16, 256), (F(1, 2), F(3, 4))):
            exact = exact_kernels(rows, q, mark, collapsed=True)
            for minimum in range(1, rows+1):
                upper = minimum_kernel(rows, q, mark, minimum)
                for h in range(minimum, rows+1):
                    for i, j in product(range(2), repeat=2):
                        self.assertTrue(upper[i][j] >= arb_fraction(exact[h][i][j]))

    def test_every_small_three_cut_matrix(self):
        for rows, mark in product((4, 8, 12), (F(1, 2), F(3, 4))):
            q = 16
            kernels = exact_kernels(rows, q, mark, collapsed=True)
            fixed = class_matrix(rows, q, mark, 3)
            averaged = shift_matrix(rows, q, mark, 3)
            for cuts in combinations(range(rows), 3):
                actual = region(cuts, rows, q, kernels)
                shifted = [[[x for x in row] for row in region(
                    [(x+s) % rows for x in cuts], rows, q, kernels)]
                    for s in range(0, rows, 4)]
                for i, j in product(range(4), repeat=2):
                    self.assertTrue(fixed[i][j] >= arb_fraction(actual[i][j]))
                    mean = sum(matrix[i][j] for matrix in shifted)/len(shifted)
                    self.assertTrue(averaged[i][j] >= arb_fraction(mean))


if __name__ == '__main__':
    unittest.main()
