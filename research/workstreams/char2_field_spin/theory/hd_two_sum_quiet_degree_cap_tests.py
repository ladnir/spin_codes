from fractions import Fraction as F
from itertools import product
from math import comb, log2
from unittest import TestCase, main

import numpy as np

import hd_two_sum_quiet_degree_cap as gate


class QuietDegreeCapTests(TestCase):
    def test_cap_bounds_finite_bad_count_per_degree(self):
        z, limit = F(1, 2), 2
        # Each degree class has its own output-weight distribution.
        classes = [[(F(1, 5), 0), (F(2, 5), 5)],
                   [(F(1, 7), 2), (F(2, 7), 3)],
                   [(F(1, 9), 7), (F(3, 9), 8)]]
        outer = [sum(mass for mass, _ in row) for row in classes]
        marked = [sum(mass*z**weight for mass, weight in row) for row in classes]
        bad = [sum(mass for mass, weight in row if weight <= limit) for row in classes]
        capped, info = gate.degree_cap(np.log2(list(map(float, outer))),
            np.log2(list(map(float, marked))), -limit*log2(z), 0)
        self.assertTrue(np.all(2**capped+1e-15 >= list(map(float, bad))))
        self.assertLessEqual(info['capped_tail_log2'], info['marked_tail_log2'])
        self.assertLessEqual(info['capped_tail_log2'], info['unweighted_tail_log2'])

    def test_degreewise_minimum_is_stronger_than_minimum_of_tails(self):
        _, info = gate.degree_cap(np.log2([1., 100.]), np.log2([100., 1.]), 0., 0)
        self.assertAlmostEqual(info['capped_tail_log2'], 1.)
        self.assertGreater(min(info['unweighted_tail_log2'], info['marked_tail_log2']), 6.)

    def test_outer_word_sums_fixed_rowsets_once_not_by_column(self):
        size, arity, columns = 4, 3, 2
        exact = [0]*(columns*arity+1)
        for rowset in range(1 << size):
            if rowset.bit_count() != arity:
                continue
            subsets = [mask for mask in range(1 << size) if mask & ~rowset == 0]
            for masks in product(subsets, repeat=columns):
                mass = 1
                for mask in masks:
                    mass *= gate.zero_two_moments(mask.bit_count())
                exact[sum(mask.bit_count() for mask in masks)] += mass
        expected = [log2(x) if x else -np.inf for x in exact]
        actual = gate.outer_word(1, columns, arity, size=size)
        np.testing.assert_allclose(actual, expected, atol=2e-14)
        self.assertAlmostEqual(actual[0], log2(comb(size, arity)))
        # Independent row-set choices for each column would give 16 here.
        self.assertNotEqual(2**actual[0], comb(size, arity)**columns)

    def test_guards(self):
        with self.assertRaises(ValueError):
            gate.outer_word(2, 2, 5, size=4)
        with self.assertRaises(ValueError):
            gate.degree_cap([0.], [0., 1.], 0., 0)
        with self.assertRaises(ValueError):
            gate.run(max_seconds=31)


if __name__ == '__main__':
    main()
