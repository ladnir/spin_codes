"""Exact matrix and outward-arithmetic checks for dense marker coverage."""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import unittest
from unittest.mock import patch

from flint import arb, ctx, fmpq, fmpq_mat

import dense_complete as dense


def rational(value):
    value = F(value)
    return fmpq(value.numerator, value.denominator)


def exact_mixture(blocks, probability):
    return dense.marker_mixture(blocks, probability, matrix=fmpq_mat,
                                rational=rational, rounding=lambda x: x)


class DenseCompleteTests(unittest.TestCase):
    def setUp(self):
        self.previous = ctx.prec
        ctx.prec = 256

    def tearDown(self):
        ctx.prec = self.previous

    def test_noncommuting_coefficient_domination_retains_region_state(self):
        # Enumerate four potential slots, grouped into two ordered pairs.
        blocks = [fmpq_mat([[1, 1], [0, 1]]), fmpq_mat([[1, 0], [1, 1]]),
                  fmpq_mat([[2, 1], [0, 1]])]
        self.assertNotEqual(blocks[0]*blocks[1], blocks[1]*blocks[0])
        identity = fmpq_mat([[1, 0], [0, 1]])
        for p in (F(1, 5), F(1, 2), F(4, 5)):
            mixed = exact_mixture(blocks, p)
            for q in range(5):
                fixed = fmpq_mat(2, 2)
                for assignment in product((0, 1), repeat=4):
                    if sum(assignment) == q:
                        value = identity
                        for offset in (0, 2):
                            value = value*blocks[sum(assignment[offset:offset+2])]
                        fixed += value
                fixed /= comb(4, q)
                probability = F(comb(4, q))*p**q*(1-p)**(4-q)
                envelope = mixed**2/rational(probability)
                for regions in (1, 3):
                    left, right = fixed**regions, envelope**regions
                    self.assertTrue(all(left[i, j] <= right[i, j]
                                        for i in range(2) for j in range(2)))
        # The matrix product keeps state at seams; scalar region resets differ.
        powered = mixed**6
        self.assertNotEqual(sum(powered[0, j] for j in range(2)),
                            sum((mixed**2)[0, j] for j in range(2))**3)

    def test_log_formula_matches_direct_fractional_scalar(self):
        groups, regions, cutoff = 4, 3, 7
        beta, tilt, moment = F(5, 2), F(1, 20), F(3, 5)**6
        for q in range(1, 5):
            for p in (F(1, 5), F(1, 2), F(4, 5)):
                for alpha in (F(7, 10), F(1)):
                    actual = dense.coefficient_log_bound(q, p, dense.replay.aq(moment).log(),
                        arb(comb(groups, q)).log(), tilt, alpha, groups=groups,
                        regions=regions, cutoff=cutoff, beta=beta).exp()
                    mass = F(comb(groups, q))*p**q*(1-p)**(groups-q)
                    expected = (dense.replay.aq(F(comb(groups, q))*moment/mass**regions)
                        * dense.replay.aq(beta)**dense.replay.aq(alpha*q)
                        * dense.replay.aq(alpha*tilt*cutoff).exp())
                    self.assertTrue((actual/expected-1).abs_upper() < arb(2)**-220)
                    self.assertGreaterEqual(dense.replay.fraction(dense.replay.endpoint(actual)),
                                            dense.replay.fraction(dense.replay.endpoint(expected.lower())))

    def test_probability_one_endpoint_and_rejections(self):
        blocks = [fmpq_mat([[1, j], [j+1, 2]]) for j in range(3)]
        self.assertEqual(exact_mixture(blocks, F(0)), blocks[0])
        self.assertEqual(exact_mixture(blocks, F(1)), blocks[-1])
        actual = dense.coefficient_log_bound(4, F(1), arb(5).log(), arb(0),
            F(1, 20), F(7, 10), groups=4, regions=3, cutoff=7, beta=F(5, 2))
        expected = arb(5).log()+dense.replay.aq(F(49, 200))+dense.replay.aq(F(14, 5))*dense.replay.aq(F(5, 2)).log()
        self.assertTrue((actual-expected).abs_upper() < arb(2)**-220)
        for q, p in ((3, 1), (4, 0), (0, F(1, 2)), (5, F(1, 2))):
            with self.assertRaises(ValueError):
                dense.coefficient_log_bound(q, p, arb(0), arb(0), F(1, 20), F(1), groups=4)

    def test_validation_precedes_loading_physical_receipts(self):
        with patch.object(Path, 'exists', return_value=False), patch.object(dense.replay, 'load_local') as load:
            for options in (dict(minimum=2), dict(maximum=8193), dict(minimum=513, maximum=512),
                            dict(alphas=['0']), dict(markers=['1']), dict(markers=['0']),
                            dict(precision=128)):
                with self.assertRaises(ValueError):
                    dense.run(['unused.json'], 'unused-output.json', **options)
            load.assert_not_called()


if __name__ == '__main__':
    unittest.main()
