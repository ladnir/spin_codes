from fractions import Fraction as F
import unittest

import search_ladder_dense_balanced as search


class Coordinates:
    def density(self,a,b):
        return F(1,4)+F(3,4)*a,F(1,4)+F(3,4)*b


class BalancedTests(unittest.TestCase):
    def test_both_dimensions_shrink(self):
        node = search.old.box(512,8192,F(0),F(1),{})
        axes = set()
        for _ in range(12):
            axis = search.split_axis(node,Coordinates())
            axes.add(axis)
            node = search.old.children(node,axis)[1]
        self.assertEqual(axes,{'q','v'})
        lo,hi,a,b = search.old.geometry(node)
        self.assertLess(hi-lo,256)
        self.assertLess(b-a,F(1,16))

    def test_stalled_v1_geometry_splits_occupancy(self):
        node = search.old.box(4353,8192,F(934895,1000000),F(934896,1000000),{})
        self.assertEqual(search.split_axis(node,Coordinates()),'q')


if __name__ == '__main__':
    unittest.main()
