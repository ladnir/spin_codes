from itertools import product
from math import log2
import unittest

from hd_pair_scalar_control import activity_tail_log, load_bank


class ControlTests(unittest.TestCase):
    def test_positive_activity_convolution(self):
        blocks=((1.,2.,3.),(1.,5.,7.),(1.,11.,13.))
        for minimum in range(7):
            expected=sum(blocks[0][a]*blocks[1][b]*blocks[2][c]
                         for a,b,c in product(range(3),repeat=3) if a+b+c>=minimum)
            self.assertAlmostEqual(activity_tail_log(blocks,minimum=minimum),log2(expected))
    def test_empty_activity_and_overflow_safe(self):
        self.assertEqual(activity_tail_log(((1.,0.,0.),)*8),float('-inf'))
        self.assertAlmostEqual(activity_tail_log(((1.,0.,2.**256),)*8,minimum=16),2048.)
    def test_pinned_regional_scope(self):
        saved,curves,prefix=load_bank()
        self.assertEqual(curves.shape,(4,19))
        self.assertEqual(saved['source_class'],'Fmax<=1')
        self.assertAlmostEqual(prefix,7932.104519645711)


if __name__=='__main__':
    unittest.main()
