from fractions import Fraction as Q
from math import comb
import unittest
import numpy as np
from feedback_exact import fwht,annihilators,distribution,census


class ExactFeedbackTests(unittest.TestCase):
    def test_transform_and_overflow(self):
        x=np.array([1,-2,3,7,5,0,-4,11],dtype=np.int64)
        self.assertEqual(fwht(fwht(x)).tolist(),(8*x).tolist())
        with self.assertRaises(OverflowError):fwht([2**62,2**62])
        with self.assertRaises(ValueError):fwht([1,2,3])

    def test_every_small_feedback_target(self):
        for columns,bits in (([1,2,4,3,5,7,6,1],3),([1,2,4,8]*3,4),([0]*8,2)):
            windows=len(columns)//4;histograms=[[0]*(1<<bits) for _ in range(windows+1)]
            for x in range(1<<len(columns)):
                j=sum(bool((x>>(4*w))&15) for w in range(windows));target=0
                for b,c in enumerate(columns):
                    if x>>b&1:target^=c
                histograms[j][target]+=1
            annihilated=annihilators(columns,bits);summary=census(columns,bits,windows)
            for j,expected in enumerate(histograms):
                actual,denominator=distribution(annihilated,windows,j)
                self.assertEqual(actual.tolist(),expected)
                self.assertEqual(denominator,comb(windows,j)*15**j)
                self.assertEqual(summary[j]['zero'],expected[0])
                self.assertEqual(summary[j]['nonzero_peak'],max(expected[1:]))

    def test_large_count_guard(self):
        with self.assertRaises(OverflowError):distribution(np.array([32,32]),32,16)
        with self.assertRaises(ValueError):distribution(np.array([2,1]),2,3)
        with self.assertRaises(ValueError):annihilators([1,2,3],2)


if __name__=='__main__':unittest.main()
