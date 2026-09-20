"""Checks for support geometry and exact outer-certificate partition replay."""
import unittest
from unittest.mock import patch
from fractions import Fraction as F

import outer_majorant as outer


class RefinedTests(unittest.TestCase):
    def test_lower_envelope_and_inheritance(self):
        segments,new=outer.supports()
        self.assertEqual(len(new),5)
        self.assertEqual(len(segments),20)
        previous=outer.screen.load_segments(outer.screen.FROZEN/'golay_ba3_concave_majorant.json')
        for a,b in zip(segments,segments[1:]):
            self.assertEqual(a.omega_hi,b.omega_lo)
            self.assertGreater(a.slope,b.slope)
            x=a.omega_hi
            self.assertEqual(a.slope*x+a.intercept,b.slope*x+b.intercept)
        for segment in segments:
            for x in (segment.omega_lo,segment.omega_hi):
                self.assertLessEqual(segment.slope*x+segment.intercept,
                                     min(s.slope*x+s.intercept for s in previous))

    def test_partition_rejects_gaps_duplicates_overlaps(self):
        segment=outer.supports()[1][0]
        leaf=lambda path:dict(path=path,upper=-1.)
        with patch.object(outer.old.Box,'objective_upper',return_value=-1.):
            self.assertEqual(outer.replay_partition(segment,[leaf('0'),leaf('1')]),-1.)
            for paths in [['0'],['0','0','1'],['','0'],['00','01','1','10']]:
                with self.assertRaises(AssertionError):
                    outer.replay_partition(segment,[leaf(p) for p in paths])

    def test_partition_checks_inequality(self):
        segment=outer.supports()[1][0]
        with patch.object(outer.old.Box,'objective_upper',return_value=F(0)):
            with self.assertRaises(AssertionError):
                outer.replay_partition(segment,[dict(path='',upper=-1.)])


if __name__=='__main__':unittest.main()
