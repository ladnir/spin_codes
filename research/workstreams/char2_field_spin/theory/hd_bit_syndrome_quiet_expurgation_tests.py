from fractions import Fraction as F
from math import comb
import unittest
from hd_bit_syndrome_quiet_expurgation import event_upper, quiet_histogram


class QuietExpurgationTests(unittest.TestCase):
    def test_exact_binary_affine_rank_histogram(self):
        expected = {(1, 0): 16, (2, 0): 120, (3, 0): 560,
                    (4, 0): 1680, (4, 1): 140,
                    (5, 0): 2688, (5, 1): 1680,
                    (6, 1): 7168, (6, 2): 840,
                    (7, 2): 11200, (7, 3): 240}
        self.assertEqual(dict(quiet_histogram()), expected)
        per, total = event_upper(expected)
        self.assertEqual(per[4], F(8192*140*(256**16-1), 255**33))
        self.assertEqual(per[1]+per[2]+per[3], 0)
        self.assertLess(total, F(1, 2**40))

    def test_zero_space_and_projective_factor(self):
        per, total = event_upper({(2, 0): 1, (3, 1): 2}, groups=3, q=4,
                                 columns=2, rs_dimension=1)
        self.assertEqual(per[2], 0)
        self.assertEqual(total, F(3*2*(4**2-1), 3**4))
        with self.assertRaises(ValueError):
            quiet_histogram(8)


if __name__ == '__main__':
    unittest.main()
