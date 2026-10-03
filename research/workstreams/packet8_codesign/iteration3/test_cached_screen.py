import unittest
from fractions import Fraction

import numpy as np

import cached_screen as cs


class CachedScreenTests(unittest.TestCase):
    def test_deterministic_full_occupancy(self):
        witness=cs.route_threshold(512)
        self.assertEqual(witness['h'],6144)
        self.assertTrue(witness['exact_integer_check_passed'])
        bound=cs.cg.cap_counts.rational_bad_route_bound(512,6144,'1',cap=3)
        self.assertEqual(bound,0)

    def test_fractional_point_parser(self):
        self.assertEqual(cs.parse_point('64:.2,3/10'),(64,(Fraction(1,5),Fraction(3,10))))
        for value in ('0:.2','513:.2','64:.2,1/5','64:0'):
            with self.assertRaises(ValueError):
                cs.parse_point(value)

    def test_batched_screen_against_log(self):
        rng=np.random.default_rng(1)
        families=rng.uniform(.01,.2,size=(3,9,2,2))
        got,backend=cs.evaluate(families,3,regions=2)
        expected=np.array([cs.cd.logarithmic_moment(local,3,regions=2) for local in families])
        np.testing.assert_allclose(got,expected,rtol=1e-12,atol=1e-12)


if __name__=='__main__':
    unittest.main()
