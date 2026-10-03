"""Geometry and exact outer-count tests for the wider byte-packet experiment."""
from fractions import Fraction
from itertools import product
from math import comb
import unittest

import numpy as np

import wider_screen as screen


class WiderTests(unittest.TestCase):
    def test_geometry(self):
        wide = screen.geometry()
        self.assertEqual((wide['group_dimension'], wide['group_output_bits']), (256, 512))
        self.assertEqual(wide['parallel_RS_rows'], 8)
        self.assertEqual((wide['outer_groups'], wide['regions']), (256, 64))
        self.assertEqual(wide['physical_steps_per_region'], 32)
        self.assertEqual(wide['regions']*wide['slots_per_region']*8, wide['N'])
        self.assertEqual(wide['physical_steps_total']*64, wide['N'])
        small = screen.geometry(symbol_bits=16)
        self.assertEqual((small['outer_groups'], small['regions']), (512, 32))
        self.assertEqual(small['physical_steps_total'], wide['physical_steps_total'])
        self.assertEqual(small['cutoff'], 13107)

    def test_exact_outer_counts_and_envelope(self):
        envelope, counts = screen.outer(screen.geometry())
        self.assertEqual(len(counts), 65)
        self.assertEqual(counts[:9], (0,)*9)
        self.assertEqual(sum(counts), (1 << 256)-1)
        self.assertEqual(envelope.beta, Fraction(1 << 512, ((1 << 32)-1)**8))
        self.assertEqual(envelope.activity, Fraction(255, 256))
        self.assertTrue(all(a <= b for a, b in zip(counts, envelope.shell_caps())))

    def test_exact_symbol_shells(self):
        envelope, _ = screen.outer(screen.geometry())
        symbols = screen.local_screen.rs_outer.mds_symbol_weight_counts(1 << 32, 16, 8)
        self.assertEqual(symbols[:9], (1,)+(0,)*8)
        self.assertEqual(sum(symbols), 1 << 256)
        for h, density in enumerate(envelope.exact_pointwise_symbol_counts()):
            if h:
                self.assertEqual(density*comb(16, h)*((1 << 32)-1)**h, symbols[h])

    def test_ordered_geometry_not_independent_replacement(self):
        # The adapter calls the same generic without-replacement recurrence
        # with actual region dimensions. Verify its order on a tiny fixture.
        local = np.array([[[1., .2], [0., .7]], [[.4, .3], [.1, .8]], [[.2, .6], [.3, .9]]])
        direct = np.zeros((7, 2, 2))
        for occupancies in product(range(3), repeat=3):
            value, ways = np.eye(2), 1
            for j in occupancies:
                value = value @ local[j]
                ways *= comb(2, j)
            q = sum(occupancies)
            direct[q] += ways/comb(6, q)*value
        actual, _ = screen.local_screen.placement(local, 6, epochs=3, windows=2)
        np.testing.assert_allclose(np.exp(actual), direct, rtol=2e-13, atol=2e-13)


if __name__ == '__main__':
    unittest.main()
