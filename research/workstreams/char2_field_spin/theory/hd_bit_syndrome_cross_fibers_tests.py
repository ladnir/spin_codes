from itertools import combinations
import unittest

from gf8_packet_maps import matvec, rank
from hd_bit_syndrome_outer import affine_matrix, matrix, quiet_support_counts
from hd_bit_syndrome_singleton_cosets import counts_for
from hd_bit_syndrome_cross_fibers import (
    cross_fiber, information_sets, inside_support_counts, inverse_matrix,
    outside_fiber, supported_count,
)


class CrossFiberTests(unittest.TestCase):
    CASES = (
        (tuple(range(7)), 4, 3, 3),
        ((0, 1, 2, 3, 4, 5, 8), 5, 2, 2),
        (tuple(range(7, 16)), 4, 5, 4),
        ((6, 7, 9, 10, 11, 12, 13, 14, 15), 5, 4, 4),
    )

    def test_representative_generic_and_exceptional_arities(self):
        for a, crossrank, zero_dim, quiet_dim in self.CASES:
            f = cross_fiber(a)
            self.assertEqual(f['cross_rank'], crossrank)
            self.assertEqual(f['inside_zero_fiber_dimension'], zero_dim)
            self.assertEqual(f['quiet_dimension'], quiet_dim)
            self.assertEqual(rank(f['outside_image_basis']), crossrank)
            self.assertEqual(rank(f['inside_zero_fiber_basis']), zero_dim)
            self.assertFalse(f['metadata']['complete_outer_kernel'])

    def test_actual_a7_image_not_full_outside_affine_code(self):
        f = cross_fiber(tuple(range(7)))
        outside_v = [affine_matrix()[i] for i in f['outside']]
        self.assertEqual(rank(outside_v), 5)
        self.assertEqual(f['outside_image_dimension'], 4)
        self.assertTrue(f['zero_outside_equals_quiet'])
        self.assertTrue(all(outside_fiber(f, z) is not None for z in f['outside_image_basis']))
        ambient_columns = [tuple(row[j] for row in outside_v) for j in range(5)]
        missing = next(z for z in ambient_columns if outside_fiber(f, z) is None)
        self.assertEqual(inside_support_counts(f, missing), {})

    def test_a9_nonquiet_zero_outside_branch(self):
        f = cross_fiber(tuple(range(7, 16)))
        self.assertEqual(f['zero_outside_excess_dimension'], 1)
        self.assertFalse(f['zero_outside_equals_quiet'])
        witness = next(v for v in f['input_kernel_basis'] if any(matvec(f['quiet_constraints'], v)))
        x = [0] * 16
        for i, value in zip(f['source'], witness):
            x[i] = value
        y = matvec(matrix(), x)
        self.assertTrue(all(y[i] == 0 for i in f['outside']))
        self.assertNotEqual(tuple(x), y)
        counts = inside_support_counts(f)
        self.assertEqual(sum(counts.values()), 256 ** 5)
        quiet = quiet_support_counts()
        self.assertEqual(sum(quiet[mask] for mask in range(65536)
                             if not mask & ~sum(1 << i for i in f['source'])), 256 ** 4)

    def test_independent_forward_constraints(self):
        h = matrix()
        for a, _, _, _ in self.CASES:
            f = cross_fiber(a)
            x = tuple((37 * i + 11) % 256 for i in range(len(a)))
            z = matvec(f['cross_matrix'], x)
            local = outside_fiber(f, z)
            self.assertIsNotNone(local)
            for w in range(4):
                for s in combinations(a, w):
                    rows = f['outside'] + tuple(i for i in a if i not in s)
                    m = [[h[i][j] for j in a] for i in rows]
                    rhs = list(z) + [0] * (len(a)-len(s))
                    r = rank(m)
                    expected = 256 ** (len(a)-r) if rank([row+[v] for row,v in zip(m,rhs)]) == r else 0
                    self.assertEqual(supported_count(f, s, z), expected)
            self.assertEqual(supported_count(f, a, z), 256 ** f['input_kernel_dimension'])

    def test_information_sets_and_shifted_singleton_recovery(self):
        for a, _, _, _ in self.CASES:
            f = cross_fiber(a)
            info = information_sets(f)
            self.assertTrue(info)
            self.assertTrue(all(len(s) == f['inside_zero_fiber_dimension'] for s in info))
            for s in info:
                cols = tuple(a.index(i) for i in s)
                self.assertEqual(rank([[row[i] for i in cols] for row in f['inside_zero_fiber_basis']]), len(s))
            counts = inside_support_counts(f)
            self.assertEqual(sum(counts.values()), 256 ** f['inside_zero_fiber_dimension'])
        a, j = (0, 1, 4, 6, 8, 10, 12, 14), 2
        f = cross_fiber(a)
        z = tuple(int(i == j) for i in f['outside'])
        self.assertEqual(inside_support_counts(f, z), counts_for(a, j))

    def test_all_arity_algebra_edges_and_resource_guard(self):
        self.assertEqual(inverse_matrix(), matrix())
        one = cross_fiber((0,))
        self.assertEqual(one['inside_zero_fiber_dimension'], 0)
        self.assertEqual(information_sets(one), ((),))
        self.assertEqual(inside_support_counts(one), {0: 1})
        all_rows = cross_fiber(tuple(range(16)))
        self.assertEqual(all_rows['outside_image_dimension'], 0)
        self.assertEqual(all_rows['inside_zero_fiber_dimension'], 16)
        self.assertEqual(supported_count(all_rows, (0, 3, 9)), 256 ** 3)
        with self.assertRaises(ValueError):
            inside_support_counts(all_rows)
        with self.assertRaises(ValueError):
            inside_support_counts(one, max_subsets=1)
        with self.assertRaises(ValueError):
            cross_fiber((2, 1))
        with self.assertRaises(ValueError):
            outside_fiber(one, (0,))


if __name__ == '__main__':
    unittest.main()
