import unittest
from fractions import Fraction as Q
from math import log

from packet_outer_cost_whole import outer_options, exact_union


class CostWholeTests(unittest.TestCase):
    def test_parallel_outer_domination(self):
        ref, options = outer_options(512)
        self.assertEqual((ref.message_bits, ref.output_bits, ref.regions), (512, 1024, 256))
        self.assertEqual(len(options), 2)
        ratio = Q(options[1]['beta'])/Q(options[0]['beta'])
        self.assertEqual(ratio, Q((1 << 32)-1, (1 << 32)+1)**8)
        ref, options = outer_options(256)
        self.assertEqual(Q(options[0]['beta']), Q(options[1]['beta']))

    def test_exact_union_and_margin(self):
        value, margin, passed = exact_union([[1, -50], [1, -51]])
        m, e = value
        got = Q(m << max(e, 0), 1 << max(-e, 0))
        expected = Q(3, 1 << 51)
        self.assertGreaterEqual(got, expected)
        self.assertLess(got-expected, expected*Q(1, 1 << 240))
        self.assertTrue(passed)
        _, _, passed = exact_union([[1, -40]])
        self.assertFalse(passed)

    def test_invalid_union(self):
        for values in ([], [[0, -40]], [[1, 1.0]], [[True, -40]]):
            with self.assertRaises(ValueError):
                exact_union(values)


if __name__ == '__main__':
    unittest.main()
