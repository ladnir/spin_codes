from fractions import Fraction as F
import unittest

from combine_ladder_dense import full_union


class UnionTests(unittest.TestCase):
    def test_exact_sum(self):
        self.assertEqual(full_union(8,{1:F(1,32),2:F(1,64)},[3,8],F(1,128)),F(7,128))

    def test_no_gap_or_wrong_endpoint(self):
        for values,interval in (({1:F(1,32)},[3,8]),({1:F(1,32),2:F(1,64)},[3,7]),
                                ({1:F(1,32),2:F(1,64),3:F(1,128)},[3,8])):
            with self.assertRaises(ValueError):
                full_union(8,values,interval,F(1,256))


if __name__ == '__main__':
    unittest.main()
