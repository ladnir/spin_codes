from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest

import hd_lambda_renyi_support as kernel


def mul(a, b):
    value = (a if b & 1 else 0) ^ ((a << 1) if b & 2 else 0)
    return value ^ (7 if value & 4 else 0)


def rank4(matrix, width):
    rows, used = [list(row) for row in matrix], 0
    for col in range(width):
        found = next((i for i in range(used, len(rows)) if rows[i][col]), None)
        if found is None:
            continue
        rows[used], rows[found] = rows[found], rows[used]
        inverse = (0, 1, 3, 2)[rows[used][col]]
        rows[used] = [mul(x, inverse) for x in rows[used]]
        for i in range(len(rows)):
            if i != used:
                factor = rows[i][col]
                rows[i] = [x ^ mul(factor, y) for x, y in zip(rows[i], rows[used])]
        used += 1
    return used


def brute(rows, width):
    counts = [0]*(1 << width)
    matrix = tuple(tuple(2 ^ (row >> j & 1) for j in range(width)) for row in rows)
    for x in product(range(4), repeat=width):
        syndrome = []
        for row in matrix:
            total = 0
            for c, value in zip(row, x):
                total ^= mul(c, value)
            syndrome.append(total)
        if not any(syndrome):
            counts[sum(bool(value) << j for j, value in enumerate(x))] += 1
    return tuple(counts)


class RenyiSupportTests(unittest.TestCase):
    def test_all_small_binary_matrices_against_GF4(self):
        for rows in product(range(8), repeat=2):
            data = kernel.support_counts_from_binary(rows, 3, q=4)
            self.assertEqual(data['counts'], brute(rows, 3))
            for subset, rank in enumerate(data['subset_ranks']):
                cols = tuple(i for i in range(3) if subset >> i & 1)
                actual = tuple(tuple(2 ^ (row >> j & 1) for j in cols) for row in rows)
                self.assertEqual(rank, rank4(actual, len(cols)))

    def test_full_diagonal_is_not_cross_only_matrix(self):
        complement = (0, 1, 2, 3)
        rows = kernel.parity_binary_rows(complement)
        for o, row in zip(complement, rows):
            self.assertEqual(row, sum((int(o == i) ^ ((o & i).bit_count() & 1)) << i for i in range(16)))
        # A six-coordinate restriction permits a tiny full GF4 enumeration.
        local = tuple(row & 63 for row in rows)
        data = kernel.support_counts_from_binary(local, 6, q=4)
        self.assertEqual(data['counts'], brute(local, 6))
        cross_only = tuple(row ^ (1 << o) for o, row in zip(complement, local))
        self.assertNotEqual(data['counts'], kernel.support_counts_from_binary(cross_only, 6, q=4)['counts'])

    def test_exact_density_moment_direct_uniform_space(self):
        rows, width = (0b1011, 0b0101), 4
        data = kernel.support_counts_from_binary(rows, width, q=4)
        for order in (1, 2, 3, 4):
            direct = F(0)
            for y in product(range(4), repeat=width):
                support = sum(bool(v) << i for i, v in enumerate(y))
                likelihood = F(4**data['codimension']*data['counts'][support], 3**support.bit_count())
                direct += likelihood**order / 4**width
            self.assertEqual(kernel.density_moment_from_counts(data['counts'], data['codimension'], q=4, order=order), direct)

    def test_uniform_and_zero_code_endpoints(self):
        full = kernel.support_counts_from_binary((), 4, q=256)
        self.assertEqual(full['counts'], tuple(255**s.bit_count() for s in range(16)))
        self.assertEqual(kernel.density_moment_from_counts(full['counts'], 0), 1)
        zero = kernel.support_counts_from_binary((1, 2, 4, 8), 4, q=4)
        self.assertEqual(zero['counts'], (1,)+(0,)*15)
        self.assertEqual(kernel.density_moment_from_counts(zero['counts'], 4, q=4), 4**8)

    def test_permuted_parity_rows_preserve_labels_and_counts(self):
        a = kernel.support_counts_from_binary((0b00101, 0b11100), 5, q=4)
        b = kernel.support_counts_from_binary((0b11100, 0b00101), 5, q=4)
        self.assertEqual(a, b)
        permutation = (3, 0, 4, 1, 2)
        rows = tuple(sum(((row >> i) & 1) << permutation[i] for i in range(5)) for row in (0b00101, 0b11100))
        moved = kernel.support_counts_from_binary(rows, 5, q=4)
        for s, count in enumerate(a['counts']):
            target = sum(((s >> i) & 1) << permutation[i] for i in range(5))
            self.assertEqual(moved['counts'][target], count)
        self.assertEqual(kernel.density_moment_from_counts(a['counts'], a['codimension'], q=4),
                         kernel.density_moment_from_counts(moved['counts'], moved['codimension'], q=4))

    def test_guards(self):
        for rows, width in (((1,)*5, 4), ((8,), 3), ((True,), 3), ((), 17)):
            with self.assertRaises(ValueError):
                kernel.subset_ranks(rows, width)
        for complement in ((0, 0), (16,), tuple(range(5)), (True,)):
            with self.assertRaises(ValueError):
                kernel.parity_binary_rows(complement)
        with self.assertRaises(ValueError):
            kernel.support_counts_from_binary((0,), 2, q=8)
        with self.assertRaises(ValueError):
            kernel.density_moment_from_counts((1, 2), 0, q=4)


if __name__ == '__main__':
    unittest.main()
