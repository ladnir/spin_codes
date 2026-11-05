from copy import deepcopy
import unittest

import numpy as np

from hd_shifted_4051_extend import recovered_bank, extend
from hd_shifted_bulk_search import LOG_DENSITY, PositiveWitness
from hd_shifted_cover_backend import profile_count
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION


class FrozenBoundaryTests(unittest.TestCase):
    def fixture(self):
        z, alpha = .9, .1
        odds = ((1., 2.), (3., 4.))
        (s0, f0), (s1, f1) = odds
        slopes = (alpha*LOG_DENSITY+1-128*(s0+s1), -2+256*s0-128*f0, 256*s1-128*f1)
        witness = PositiveWitness(z, alpha, odds, (4096, 1024, 1024), -1e8, slopes, 0.)
        boxes = [((4096, 4096), (0, 1023), (0, 4096)),
                 ((4096, 4096), (1024, 4096), (0, 4096))]
        return dict(construction=SWAP_CONSTRUCTION, requested_weight_range=[4096, 4096],
                    witnesses=[dict(witness.serializable(), construction=SWAP_CONSTRUCTION,
                                    pair_swaps=True, orientation_aggregation='independent-fair-swaps')],
                    accepted_boxes=[dict(ranges=box, witness_id=0,
                                         union_log2_upper_screen=witness.evaluate(box)[0]) for box in boxes])

    def test_recovery_cancels_geometry_without_building_an_operator(self):
        plan = self.fixture()
        _, intercepts, discrepancy = recovered_bank(plan)
        self.assertLess(discrepancy, 1e-7)
        self.assertAlmostEqual(intercepts[0], -1e8+1e-7, places=6)
        broken = deepcopy(plan)
        broken['accepted_boxes'][1]['union_log2_upper_screen'] += 1
        with self.assertRaises(ArithmeticError):
            recovered_bank(broken)

    def test_extended_partition_counts_full_domain_without_certificate(self):
        result = extend(self.fixture(), first=4051, last=4095)
        self.assertEqual(result['unresolved_boxes'], [])
        self.assertEqual(result['accepted_integer_triples'], profile_count(((4051, 4095), (0, 4096), (0, 4096))))
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['reconstruction']['operator_builds'], 0)


if __name__ == '__main__':
    unittest.main()
