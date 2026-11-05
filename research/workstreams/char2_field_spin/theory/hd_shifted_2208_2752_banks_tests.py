from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import unittest

from hd_shifted_2208_2752_banks import GAPS, SOURCE, COUNT_METHOD, run
from hd_shifted_outer_density import TOEPLITZ_OUTER_MEASURE
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION


class GapBankTests(unittest.TestCase):
    def test_assignment_is_disjoint_from_verified_2496_band(self):
        self.assertEqual(GAPS, {2208: (2154, 2274), 2752: (2684, 2827)})
        for center in (0, 2496, True):
            with self.assertRaises(ValueError):
                run(center, 'unused.json')

    def test_optional_frozen_proposals_and_receipts(self):
        if not SOURCE.exists():
            self.skipTest('Optional frozen source bank absent; tests never build a bank')
        source = json.loads(SOURCE.read_text(encoding='utf-8'))
        alpha = source['fixed_witness']['source_exact_witness']['alpha']
        for center, (first, last) in GAPS.items():
            stem = 'hd_shifted_potential_band_%d_%d_%d' % (center, first, last)
            proposal_path = SOURCE.with_name(stem+'.json')
            if not proposal_path.exists():
                continue
            with self.subTest(center=center):
                raw = proposal_path.read_bytes()
                proposal = json.loads(raw)
                self.assertEqual(proposal['requested_weight_range'], [first, last])
                self.assertEqual([row['active_groups'] for row in proposal['results']], list(range(first, last+1)))
                self.assertEqual(proposal['count_method'], COUNT_METHOD)
                self.assertFalse(proposal['directed'])
                self.assertEqual(proposal['z_choices'], 1)
                self.assertEqual(proposal['fixed_witness']['outer_measure'], TOEPLITZ_OUTER_MEASURE)
                exact = proposal['fixed_witness']['source_exact_witness']
                self.assertEqual(exact['alpha'], alpha)
                self.assertEqual(exact['density_scaled_from_weight'], 2496)
                self.assertLessEqual(len(exact['markers']), 100)
                self.assertNotIn('alpha_and_potential', exact)
                receipts = []
                for precision in (256, 384):
                    path = SOURCE.with_name(stem+'_directed%d.json' % precision)
                    if not path.exists():
                        continue
                    receipt = json.loads(path.read_text(encoding='utf-8'))
                    self.assertEqual(receipt['proposal_sha256'], sha256(raw).hexdigest())
                    self.assertEqual(receipt['weight_range'], [first, last])
                    self.assertEqual(receipt['weights_replayed'], last-first+1)
                    self.assertEqual([row['active_groups'] for row in receipt['results']],
                                     list(range(first, last+1)))
                    self.assertEqual(receipt['construction'], SWAP_CONSTRUCTION)
                    self.assertEqual(receipt['outer_measure'], TOEPLITZ_OUTER_MEASURE)
                    self.assertTrue(receipt['directed'])
                    self.assertTrue(receipt['covers_all_supports_at_each_weight'])
                    self.assertTrue(receipt['covers_all_even_odd_splits'])
                    self.assertFalse(receipt['whole_code_certificate'])
                    self.assertLess(F(receipt['log2_upper_rational']), -55)
                    receipts.append(receipt)
                if len(receipts) == 2:
                    self.assertEqual(receipts[0]['log2_upper_rational'], receipts[1]['log2_upper_rational'])
                    self.assertEqual([row['log2_upper_rational'] for row in receipts[0]['results']],
                                     [row['log2_upper_rational'] for row in receipts[1]['results']])


if __name__ == '__main__':
    unittest.main()
