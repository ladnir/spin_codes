import math
import unittest

import gf8_fullstate_outer_screen as outer
import gf8_fullstate_screen as base


class FullStateOuterTests(unittest.TestCase):
    def test_same_length_middle(self):
        point = (-1., -.4, math.log(32), -.2)
        g = outer.profile(128, 96)
        for state in (7, 8):
            for density in (.03, .10, .28):
                w = round(density*g['n'])
                self.assertAlmostEqual(outer.middle_bound(point, w, 128, 96, state),
                                       base.middle_bound(point, w, state), places=8)

    def test_same_length_sparse(self):
        point = (math.log(640), math.log(16))
        for state in (7, 8):
            for active in (1, 4):
                self.assertAlmostEqual(outer.sparse_bound(active, point, 128, 96, state),
                                       base.sparse_bound(active, point, state), places=9)

    def test_geometry_and_limits(self):
        g = outer.profile(256, 192)
        self.assertEqual(g['rows'] % 64, 0)
        self.assertGreaterEqual(g['k'], 2**20)
        self.assertEqual(g['n'], g['rows']*256)
        self.assertEqual(4*g['k'], 3*g['n'])
        with self.assertRaises(ValueError):
            outer.profile(257, 192)


if __name__ == '__main__':
    unittest.main()
