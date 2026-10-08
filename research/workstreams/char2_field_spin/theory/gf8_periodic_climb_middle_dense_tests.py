from fractions import Fraction as F
import math
import unittest

import numpy as np

from gf8_periodic_climb_profile import metadata
from gf8_periodic_climb_middle_certificate import (certify,WITNESSES_84_6667,
                                                  WITNESSES_74_56)
from gf8_periodic_climb_dense import certify as dense
from gf8_puncture_frontier import period_matrix
from gf8_rate_climb_middle_certificate import exact_transfer,exact_period_transfer,cell_upper,prefactor
from return_cover_certificate import interval


class PeriodicClimbComponentTests(unittest.TestCase):
    def test_exact_ordered_product(self):
        witness=WITNESSES_84_6667[0][1]
        for pattern in ((6,),(6,6,6,7),(5,6),(6,5)):
            exact=exact_period_transfer(witness,pattern)
            point=[math.log(float(value)) for value in
                   (witness.x,witness.y,witness.tau,witness.tilt)]
            matrix,scale=period_matrix(point,pattern)
            np.testing.assert_allclose(np.asarray(exact,dtype=float),matrix*math.exp(scale),rtol=2e-13)
        self.assertEqual(exact_period_transfer(witness,(6,)),exact_transfer(witness,6))
        self.assertNotEqual(exact_period_transfer(witness,(5,6)),exact_period_transfer(witness,(6,5)))
        with self.assertRaises(ValueError):
            exact_period_transfer(witness,())

    def test_periodic_geometry_and_integer_joins(self):
        for dimension,pattern,sparse_last,first,dense_first in (
                (84,(6,6,6,7),1070,48169,1364788),(74,(5,6),991,54559,1636762)):
            g=metadata(dimension,pattern)
            self.assertEqual(g['input_n']*sum(pattern)//(8*len(pattern)),g['output_n'])
            self.assertEqual(g['output_threshold'],(g['output_n']+9)//10-1)
            self.assertGreaterEqual((129-dimension)*(sparse_last+1),first)
            self.assertNotIn('outputs',g)
            self.assertEqual(g['packet_period'],len(pattern))
            self.assertLess(first,dense_first)

    def test_periodic_cell_dominates_singletons(self):
        g=metadata(74,(5,6))
        witness=WITNESSES_74_56[0][1]
        factor=prefactor(witness,g)
        cell=cell_upper(54559,54562,g,witness,factor)
        individual=sum((cell_upper(w,w,g,witness,factor) for w in range(54559,54563)),interval(F(0)))
        self.assertLessEqual(individual.hi,cell.hi)

    def test_complete_period4_replay(self):
        middle=certify()
        tail=dense()
        self.assertEqual(middle['weights'],[48169,1364787])
        self.assertEqual(tail['weights'],[1364788,1605632])
        self.assertGreaterEqual(middle['margin_bits_lower'],40)
        self.assertGreaterEqual(tail['margin_bits_lower'],40)
        self.assertTrue(middle['complete'] and middle['rigorous'])
        self.assertFalse(middle['whole_code'])

    def test_complete_period2_middle_replay(self):
        result=certify(74,(5,6))
        self.assertEqual(result['weights'],[54559,1636761])
        self.assertGreaterEqual(result['margin_bits_lower'],40)
        self.assertTrue(result['complete'] and result['rigorous'])
        self.assertFalse(result['whole_code'])


if __name__=='__main__':
    unittest.main()
