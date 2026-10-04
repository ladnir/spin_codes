"""Exact integer-union tests independent of the Arb evaluation."""
from fractions import Fraction as F
import unittest

import assemble_whole as assembly


def rational(pair):
    m, e = pair
    return F(m)*F(2)**e


class ExactUnionTests(unittest.TestCase):
    def test_order_matches_rationals(self):
        values = [(1,-60), (3,-62), (7,-63), (1,-10000), (1<<300,-320), (1,10000)]
        for a in values:
            for b in values:
                self.assertEqual(assembly.less(a,b), rational(a) < rational(b))

    def test_sum_retains_small_terms(self):
        values = [(1,-60), (3,-62), (7,-63), (1,-10000), (1<<300,-320), (1,10000)]
        actual = assembly.exact_sum(values)
        self.assertEqual(rational(actual), sum(map(rational,values),F(0)))
        encoded = assembly.encode(actual)
        self.assertEqual(rational((int(encoded['numerator_hex'],16),encoded['exponent'])), rational(actual))

    def test_strict_target_comparison(self):
        self.assertFalse(assembly.less((1,-50),(1,-50)))
        self.assertTrue(assembly.less(((1<<256)-1,-306),(1,-50)))
        above = assembly.exact_sum([(1,-50),(1,-10000)])
        self.assertFalse(assembly.less(above,(1,-50)))

    def test_reject_invalid_dyadics(self):
        for value in ([0,0],[-1,0],[1,1.0],[True,0],[1,10000001],[1]):
            with self.subTest(value=value), self.assertRaises(ValueError):
                assembly.dyadic(value)


if __name__ == '__main__':
    unittest.main()
