from fractions import Fraction as F
import unittest

import numpy as np

from hd_shifted_alias_bulk import AliasWitness,uniform_regional
from hd_shifted_bulk_search import PositiveWitness
from hd_shifted_cover_backend import CachedBound
from hd_shifted_frontier_screen import SAVED_WITNESSES
from hd_shifted_unoriented_screen import regional,unsigned_histograms


class UniformAliasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        saved = SAVED_WITNESSES[(128,128,100,100)]
        cls.z,cls.alpha = saved['mark'],saved['alpha']
        row = saved['marker_log_odds'][0]
        cls.odds = (float(np.logaddexp2(*row[:2])),row[2])

    def test_uniform_numerator_dominates_every_checked_target(self):
        uniform = uniform_regional(self.z,self.alpha,self.odds)
        from hd_grouped_packet_screen import log_multinomial
        for r,a in ((256,100),(260,90),(250,90),(270,100)):
            histogram = unsigned_histograms(4096,r,a,a)[0]
            local,_ = regional(self.z,self.alpha,histogram,self.odds,grid=128)
            restored = local+log_multinomial(histogram)+np.dot(histogram[1:],self.odds)
            self.assertTrue(np.all(restored<=uniform+1e-8))

    def test_alias_witness_reduces_global_numerator_and_records_replay_mode(self):
        odds = (self.odds,self.odds)
        plain = PositiveWitness.create(self.z,self.alpha,odds)
        alias = AliasWitness.create(self.z,self.alpha,odds)
        self.assertLess(alias.log2_numerator,plain.log2_numerator-2000)
        saved = alias.serializable()
        self.assertEqual(saved['coefficient_method'],'positive-alias')
        self.assertTrue(saved['uniform_aliases'])
        self.assertEqual(saved['grid'],128)

    def test_directed_uniform_box_agrees_with_float_proposal(self):
        odds = (self.odds,self.odds)
        box = ((255,257),(99,101),(99,101))
        screen,_ = AliasWitness.create(self.z,self.alpha,odds).evaluate(box)
        row = tuple(F(2**v).limit_denominator(10**9) for v in self.odds)
        backend = CachedBound(self.z,self.alpha,(row,row),grid=128,uniform_aliases=True)
        receipt = backend.receipt(box)
        self.assertAlmostEqual(float(F(receipt['log2_upper_rational'])),screen,places=3)
        self.assertTrue(backend.witness()['uniform_aliases'])
        self.assertFalse(receipt['whole_code_certificate'])

    def test_uniform_alias_requires_grid(self):
        with self.assertRaises(ValueError):
            CachedBound(F(1,2),F(1,10),((F(1),F(1)),)*2,uniform_aliases=True)

    def test_whole_box_run_count_closes_without_profile_factor(self):
        from hd_shifted_run_count import stationary_markers
        saved = SAVED_WITNESSES[(256,256,200,200)]
        row = saved['marker_log_odds'][0]
        odds = (float(np.logaddexp2(*row[:2])),row[2])
        markers = tuple(F(2**v).limit_denominator(10**9) for v in odds)
        backend = CachedBound(saved['mark'],saved['alpha'],(markers,markers),grid=128,
                              run_markers=stationary_markers(4096,512,112))
        result = backend.receipt(((504,520),(196,204),(196,204)))
        self.assertLess(F(result['log2_upper_rational']),-50)
        self.assertEqual(result['integer_profile_triples'],1377)
        self.assertIn(result['support_count_method'],('whole-box-cyclic-runs','tilted-cyclic-run-sum'))
        self.assertFalse(result['whole_code_certificate'])


if __name__=='__main__':
    unittest.main()
