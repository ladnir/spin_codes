import itertools
import unittest

import parameter_range_bridge as bridge
import type_box_coverage


class RangeBridgeTests(unittest.TestCase):
    def test_clipping_matches_integer_domain(self):
        for total in (3,6,9):
            root = dict(lower=[0,0,0],upper=[total-1,total,total],witness={'test':1})
            all_points = {p for p in itertools.product(range(total+1),repeat=3)
                          if sum(p)==total and p[0] <= total-1}
            for minimum in range(1,total+1):
                result = bridge.clip(root,total,minimum)
                expected = {p for p in all_points if total-p[0]>=minimum}
                actual = {p for p in all_points if all(a<=v<=b for a,v,b in
                          zip(result['lower'],p,result['upper']))}
                self.assertEqual(actual,expected)
                self.assertEqual(type_box_coverage.check([result],total,minimum,3),len(expected))
                self.assertEqual(root['upper'],[total-1,total,total])
                self.assertIsNot(result['witness'],root['witness'])

    def test_empty_and_singleton_intersections(self):
        singleton = dict(lower=[8,2,0],upper=[8,10,10],witness={})
        self.assertIsNone(bridge.clip(singleton,10,3))
        self.assertEqual(bridge.clip(singleton,10,2),dict(lower=[8,2,0],upper=[8,2,0],witness={}))
        inconsistent = dict(lower=[9,2,0],upper=[10,10,10],witness={})
        self.assertIsNone(bridge.clip(inconsistent,10,1))


if __name__ == '__main__':
    unittest.main()
