import unittest
from fractions import Fraction as F

import field128_dense_certificate as certificate


class Field128DenseTests(unittest.TestCase):
    def test_integer_geometry_and_coverage(self):
        g = certificate.geometry()
        self.assertEqual(g,dict(rows=12096,k=1052352,n=1548288,
                                threshold=154828,dense_first=371590))
        self.assertLess(g['k']-2**20,87*64)
        self.assertLess(10*g['threshold'],g['n'])
        self.assertGreaterEqual(10*(g['threshold']+1),g['n'])
        intervals = certificate.weight_intervals(g['n'])
        self.assertEqual(len(intervals),152)
        self.assertEqual(intervals[0][0],371590)
        self.assertEqual(intervals[-1][1],1548288)
        for index,(lo,hi) in enumerate(intervals,48):
            self.assertGreaterEqual(200*lo,index*g['n'])
            self.assertLessEqual(200*hi,(index+1)*g['n'])
        self.assertTrue(all(left[1]+1==right[0] for left,right in zip(intervals,intervals[1:])))

    def test_bracket_monotonic_endpoints_and_extra_parity(self):
        for index in (48,49,100,168,199):
            tilt = certificate.TILTS[index-48]
            lower,upper = F(index,200),F(index+1,200)
            X = (1-lower)**32
            Y = F(1,256**10)*(1+(F(256,255)*tilt-1)*upper)**32
            target = certificate.interval_base(128,87,index,tilt)
            self.assertEqual(target,X**4+Y**4/256)
            self.assertLessEqual(target,(X+Y)**4)
            for alpha in (lower,(lower+upper)/2,upper):
                actual = ((1-alpha)**128+F(1,256**41)
                          *(1+(F(256,255)*tilt-1)*alpha)**128)
                self.assertLessEqual(actual,target)
            self.assertGreaterEqual(tilt,1)

    def test_invalid_witness_cannot_certify(self):
        self.assertFalse(certificate.certify_interval(48,F(1)))
        with self.assertRaises(ValueError):
            certificate.interval_base(128,87,48,F(1,2))

    def test_full_exact_replay(self):
        result = certificate.certify()
        self.assertTrue(result['dense_rigorous'])
        self.assertTrue(result['cutoff_cover_complete'])
        self.assertFalse(result['whole_code_certified'])
        self.assertEqual(result['covered_weights'],1548288-371590+1)
        self.assertEqual(result['dense_margin_bits_lower'],13690)


if __name__ == '__main__':
    unittest.main()
