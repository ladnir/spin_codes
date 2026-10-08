from fractions import Fraction as F
from math import comb
import unittest

from flint import ctx
from gf8_punctured_dense import interval_margin as frozen_dense
from gf8_punctured_middle_certificate import (exact_transfer as frozen_transfer,
                                             WITNESSES as frozen_witnesses)
from gf8_rate_climb_dense import interval_margin,certify as certify_dense,point_frontier
from gf8_rate_climb_middle_certificate import (exact_transfer,cell_upper,prefactor,
                                              certify as certify_middle,WITNESSES_80_6,source_upper)
from gf8_rate_climb_profile import metadata
from return_cover_certificate import interval


class RateClimbComponentTests(unittest.TestCase):
    def test_frozen84_control_unchanged(self):
        self.assertEqual(exact_transfer(frozen_witnesses[0][1],7),frozen_transfer(frozen_witnesses[0][1]))
        previous=ctx.prec
        try:
            ctx.prec=192
            g=metadata(84,7)
            self.assertTrue((interval_margin(g,750000,752000,F(10))
                             -frozen_dense(g,750000,752000,F(10))).contains(0))
        finally:
            ctx.prec=previous

    def test_new_geometry_and_cell_domination(self):
        g=metadata(80,6)
        self.assertEqual((g['input_n'],g['output_n'],g['output_threshold']),
                         (1679360,1259520,125951))
        self.assertEqual(g['nominal_rate'],'5/6')
        w=WITNESSES_80_6[0][1]
        factor=prefactor(w,g)
        cell=cell_upper(50381,50384,g,w,factor)
        individual=sum((cell_upper(x,x,g,w,factor) for x in range(50381,50385)),interval(F(0)))
        self.assertLessEqual(individual.hi,cell.hi)

    def test_iterative_source_powers_exact_points(self):
        for lo,hi in ((F(3,100),F(1,20)),(F(2,5),F(9,20))):
            w=WITNESSES_80_6[0][1]
            envelope=source_upper(lo,hi,w)
            for alpha in (lo,(lo+hi)/2,hi):
                beta=w.tau*w.tilt
                exact=(1-alpha)**128+sum(F(comb(128,h))*beta**h*alpha**h*(1-alpha)**(128-h)
                                        for h in range(49,129))/w.tau**48
                self.assertLessEqual(exact,F(envelope.hi))

    def test_dense_frontier_and_directed_replay(self):
        control=point_frontier(80,6)
        self.assertGreater(control['approximate_first_density'],.79)
        self.assertLess(control['approximate_first_density'],.80)
        result=certify_dense()
        self.assertEqual(result['weights'],[1343488,1679360])
        self.assertEqual(result['completion_factor_log2'],3358720)
        self.assertGreaterEqual(result['margin_bits_lower'],40)
        self.assertTrue(result['rigorous'] and result['complete'])

    def test_middle_complete_replay(self):
        result=certify_middle()
        self.assertEqual(result['weights'],[50381,1343487])
        self.assertGreaterEqual(result['margin_bits_lower'],40)
        self.assertTrue(result['rigorous'] and result['complete'])
        self.assertFalse(result['whole_code'])

    def test_nondefault_requires_witnesses(self):
        with self.assertRaises(ValueError):
            certify_middle(88,7)
        with self.assertRaises(ValueError):
            certify_dense(first_density=.8)


if __name__=='__main__':
    unittest.main()
