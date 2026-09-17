"""Exact local cancellation checks and numerical coefficient cross-checks."""
from fractions import Fraction as F
import unittest

import numpy as np
import adaptive_step as study


class AdaptiveStepTests(unittest.TestCase):
    def test_retention_average(self):
        for h in (1,2,4,8):
            for gap in (0,1,2):
                direct = sum((F(1,2**(gap*h+h-a+b))
                              for a in range(h) for b in range(h)),F())/h**2
                self.assertEqual(direct,study.adjacent_retention(h,gap))

    def test_cancellation_against_enumerated_transvections(self):
        s,t = 3,4
        columns = [1,2,3,4]
        pairs = [(u,v) for u in range(1,2**s) for v in range(2**s)
                 if (u&v).bit_count()%2 == 0]
        p = [[F() for _ in range(2**s)] for _ in range(2**s)]
        for q in range(2**s):
            for u,v in pairs:
                p[q][q^(u if (q&v).bit_count()%2 else 0)] += F(1,len(pairs))
        for h,gap in ((1,0),(2,0),(2,1),(4,0)):
            powers = [[[F(i==j) for j in range(2**s)] for i in range(2**s)]]
            for _ in range((gap+2)*h):
                prev = powers[-1]
                powers.append([[sum((prev[i][k]*p[k][j] for k in range(2**s)),F())
                                for j in range(2**s)] for i in range(2**s)])
            actual = sum((powers[gap*h+h-a+b][x][y]
                          for a in range(h) for b in range(h)
                          for x in columns for y in columns),F())/(h*h*t*t)
            self.assertEqual(actual,study.cancel_probability(t,s,h,gap))

    def test_retained_cancellation_selects_short_boundary_pulses(self):
        for h in (1,2,4,8,16):
            mass = sum((F(1,2**(h-a+b)) for a in range(h) for b in range(h)),F())
            first = sum((F(h-a+b,2**(h-a+b)) for a in range(h) for b in range(h)),F())
            self.assertEqual(first/mass,3-F(2*h,2**h-1))

    def test_retained_pair_output_weight(self):
        a,columns = [1,2,3,4,5],[3,5,6,7,1]
        image = lambda q: sum(((q&c).bit_count()&1)<<i for i,c in enumerate(a))
        for p,q in enumerate(columns):
            for distance in range(1,8):
                inputs = [1<<p]+[0]*(distance-1)+[1<<p]
                state,total = 0,0
                for x in inputs:
                    total += (x^image(state)).bit_count()
                    if x: state ^= columns[p]
                self.assertEqual(state,0)
                self.assertEqual(total,distance*image(q).bit_count()+2-2*((image(q)>>p)&1))

    def test_short_maps_and_dimension_limits(self):
        for t,s in ((8,6),(16,10),(32,15)):
            record = study.inner(t,s)
            self.assertEqual(sum(record['spectrum'].values()),2**s-1)
            self.assertNotIn(t,record['spectrum'])
            self.assertEqual(len(set(record['feedback_columns'])),t)
        with self.assertRaises(ValueError): study.inner(32,20)

    def test_scaled_and_log_coefficient_products(self):
        record = study.inner(8,6)
        for sharp in (False,True):
            for refresh in (False,True):
                z,a = study.grid.wm.regions(*study.epoch_logs(record,np.array([-6.,-3.]),refresh,sharp),4)
                a -= np.log(4.)
                actual = study.grid.coefficients(z,a,16)
                expected = study.grid.wm.coefficients(z,a,16)
                np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-11)


if __name__ == '__main__':
    unittest.main()
