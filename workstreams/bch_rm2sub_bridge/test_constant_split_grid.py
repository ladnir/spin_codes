from fractions import Fraction as F
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, ctx

import constant_split_grid as grid
import search_split_certificates as search


class GridTests(unittest.TestCase):
    def test_ordinary_policy_rejects_constant_rows(self):
        self.assertEqual(len(grid.ordinary_bands()), 12)
        original = grid.sparse.BANDS
        for bad in (0, 256):
            with patch.object(grid.sparse, 'BANDS', [(bad,)]+list(original[1:])):
                with self.assertRaises(ValueError):
                    grid.ordinary_bands()

    def test_every_grid_offset_matches_single_fold(self):
        qmax = 7
        rng = np.random.default_rng(738)
        matrices = rng.random((qmax+1, 4, 4))
        powers = np.array([-100, -30, 0, 5, 40, 3, 20, 0], dtype=np.int64)
        left, right = np.array([.25, .75]), np.array([.8, .3])
        batches = list(grid.folds(matrices, powers, left, right))
        for q in range(1, qmax+1):
            for d, single, exponent in grid.single.constant_matrices(matrices[:q+1], powers[:q+1], left, right):
                np.testing.assert_array_equal(single, batches[d][1][q-d])
                self.assertEqual(exponent, int(batches[d][2][q-d]))

    def test_terminal_dominates_arb_all_entries_and_scales(self):
        ctx.prec = 512
        rng = np.random.default_rng(817)
        matrices = np.concatenate((rng.random((6, 4, 4)), np.zeros((1, 4, 4)),
            np.eye(4)[None, :, :], np.full((1, 4, 4), np.nextafter(0., 1.))))
        exponents = np.array([-10000, -200, 0, 3, 150, 10000, -20, 40, 0], dtype=np.int64)
        values, powers = grid.terminal_256(matrices, exponents)
        for a, e, v, p in zip(matrices, exponents, values, powers):
            exact = grid.density.engine.power(tuple(arb(float(x)) for x in a.flat), 256)
            expected = (sum(exact[:4], arb(0))*arb(2)**(256*int(e))).upper()
            actual = arb(float(v))*arb(2)**int(p)
            self.assertTrue(actual >= expected)

    def test_exact_dyadic_ceiling(self):
        for mantissa in (.5, 1., np.nextafter(1., 0.), np.nextafter(1., 2.), .12345):
            for exponent in (-200, -20, 0, 100):
                for count in (1, 3, 1 << 50, (1 << 400)+1):
                    exact = F(float(mantissa))*F(2)**exponent*count
                    power = grid.dyadic_ceiling(mantissa, exponent, count)
                    self.assertGreaterEqual(F(2)**power, exact)
                    if power > -80:
                        self.assertLess(F(2)**(power-1), exact)

    def test_missing_cases_rejected(self):
        with self.assertRaises(ValueError):
            grid.bounds([[-80]*3], 3, 3)
        with self.assertRaises(ValueError):
            grid.bounds([[-80, -80, None]], 2, 2)
        values = grid.bounds([[-80]*3, [-80]*4], 2, 3)
        self.assertEqual(values, [F(3, 1 << 80), F(4, 1 << 80)])

    def test_grid_single_evaluator_comparison(self):
        ctx.prec = 256
        spec = dict(rows=8, cutoff=2)
        # Small positive matrices keep this comparison independent of audit files.
        region = [tuple(arb((i+1)*(j+1))/256 for j in range(16)) for i in range(7)]
        ps = [grid.base.encode(F(1, 2))]*12
        values = grid.evaluate(region, ps, spec, 2, 6, -10)
        for q in range(2, 7):
            old = grid.single.evaluate(region[:q+1], ps, spec, q, -10)
            old_powers = [r['upper_power'] for r in reversed(old)]
            self.assertEqual(values[q-2], old_powers)

    def test_endpoint_pooled_fails_split_closes(self):
        # Actual selected map, actual BCH caps, no saved numerical values reused.
        ctx.prec = 256
        spec = grid.core.instance('t128_s19', 16)
        t, s, ac, kernel = grid.core.inputs.load('t128_s19')
        caps = grid.core.inputs.caps_module.caps()
        pooled_region = grid.density.engine.regions(t, s, ac, kernel, (-(arb(2)/10).exp()).exp(), 512, 512)
        pooled_ps = grid.density.choose(pooled_region, 512, caps)
        pooled = grid.density.evaluate(pooled_region, pooled_ps, 512, 512, spec, 2, caps)
        self.assertGreater(pooled[0]['upper_power'], 8000)
        best = [[None]*513]
        for tilt, q, h in search.DEFAULT_WITNESSES[4:]:
            region = grid.density.engine.regions(t, s, ac, kernel, (-(arb(tilt)/10).exp()).exp(), 512, 512)
            ps = grid.density.choose(region[h:q+1], q-h, caps)[:-1]
            grid.merge(best, grid.evaluate(region, ps, spec, 512, 512, tilt))
        self.assertEqual(best, [[-80]*513])
        self.assertLessEqual(grid.bounds(best, 512, 512)[0], F(2)**-70)


if __name__ == '__main__':
    unittest.main()
