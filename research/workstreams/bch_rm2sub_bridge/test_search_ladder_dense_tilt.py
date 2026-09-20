from fractions import Fraction as F
import unittest

import search_ladder_dense_tilt as search


class CoverTests(unittest.TestCase):
    def test_mixed_splits_exact_cover(self):
        root = search.box(3,10,F(0),F(1),{})
        left,right = search.children(root,'q')
        a,b = search.children(left,'v')
        leaves = {'00':a,'01':b,'1':right}
        self.assertTrue(search.check_partition(leaves,{'':'q','0':'v'},3,10)['complete'])
        bad = dict(leaves)
        del bad['00']
        with self.assertRaises(ValueError):
            search.check_partition(bad,{'':'q','0':'v'},3,10)

    def test_wrong_geometry_rejected(self):
        root = search.box(3,10,F(0),F(1),{})
        a,b = search.children(root,'v')
        b['lo'] = 4
        with self.assertRaises(ValueError):
            search.check_partition({'0':a,'1':b},{'':'v'},3,10)


if __name__ == '__main__':
    unittest.main()
