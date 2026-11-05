from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import unittest

from hd_shifted_2496_bank import CENTER, CENTERS, Z_CHOICES, selected_interval
from hd_shifted_outer_density import TOEPLITZ_OUTER_MEASURE
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION


class Bank2496Tests(unittest.TestCase):
    def test_explicit_work_limits(self):
        self.assertEqual(CENTER, 2496)
        self.assertEqual(len(Z_CHOICES), 3)
        self.assertLessEqual(len(CENTERS), 100)
        self.assertEqual(CENTERS[0], 0)
        self.assertEqual(CENTERS[-1], CENTER//2)
        self.assertTrue(all(0 < z < 1 for z in Z_CHOICES))

    def test_selects_only_contiguous_passing_component(self):
        proposal = dict(results=[dict(active_groups=r, log2_upper_screen=v)
            for r, v in zip(range(2492, 2501), (-100, -40, -80, -70, -90, -65, -10, -80, -90))],
            fixed_witness={'source_exact_witness': {'z': '1/2'}}, directed=False,
            whole_code_certificate=False)
        selected = selected_interval(proposal)
        self.assertEqual(selected['requested_weight_range'], [2494, 2497])
        self.assertEqual(selected['weights_evaluated'], 4)
        self.assertFalse(selected['directed'])
        self.assertFalse(selected['whole_code_certificate'])
        self.assertGreater(selected['union_log2_upper_screen'], -65)
        self.assertLess(selected['union_log2_upper_screen'], -64)
        selected['fixed_witness']['source_exact_witness']['z'] = '1/3'
        self.assertEqual(proposal['fixed_witness']['source_exact_witness']['z'], '1/2')

    def test_rejects_missing_weight_and_returns_none_for_failed_center(self):
        with self.assertRaises(ValueError):
            selected_interval(dict(results=[dict(active_groups=r, log2_upper_screen=-100)
                                            for r in (2495, 2497)]))
        self.assertIsNone(selected_interval(dict(results=[dict(active_groups=2496, log2_upper_screen=-1)])))

    def test_frozen_proposal_scope(self):
        path = Path(__file__).with_name('hd_shifted_potential_band_2496_choice0_2154_2827_selected.json')
        if not path.exists():
            self.skipTest('Optional frozen research artifact absent')
        proposal = json.loads(path.read_text(encoding='utf-8'))
        self.assertEqual(proposal['requested_weight_range'], [2275, 2683])
        self.assertEqual([row['active_groups'] for row in proposal['results']], list(range(2275, 2684)))
        self.assertTrue(all(row['log2_upper_screen'] <= -60 for row in proposal['results']))
        self.assertEqual(proposal['construction'], SWAP_CONSTRUCTION)
        self.assertFalse(proposal['directed'])
        self.assertFalse(proposal['whole_code_certificate'])
        fixed = proposal['fixed_witness']
        self.assertEqual(fixed['outer_measure'], TOEPLITZ_OUTER_MEASURE)
        self.assertEqual(len(fixed['source_exact_witness']['markers']), 50)
        self.assertEqual(F(fixed['source_exact_witness']['z']), F(Z_CHOICES[0]))

    def test_optional_cross_precision_receipts(self):
        paths = [Path(__file__).with_name('hd_shifted_potential_band_2496_2275_2683_directed%d.json' % precision)
                 for precision in (256, 384)]
        if not all(path.exists() for path in paths):
            self.skipTest('Optional directed receipts absent; tests never launch a replay')
        receipts = [json.loads(path.read_text(encoding='utf-8')) for path in paths]
        proposal = paths[0].with_name('hd_shifted_potential_band_2496_choice0_2154_2827_selected.json')
        digest = sha256(proposal.read_bytes()).hexdigest()
        for receipt in receipts:
            self.assertEqual(receipt['proposal_sha256'], digest)
            self.assertEqual(receipt['weight_range'], [2275, 2683])
            self.assertEqual(receipt['weights_replayed'], 409)
            self.assertTrue(receipt['every_integer_weight_replayed'])
            self.assertTrue(receipt['covers_all_supports_at_each_weight'])
            self.assertTrue(receipt['covers_all_even_odd_splits'])
            self.assertTrue(receipt['directed'])
            self.assertFalse(receipt['whole_code_certificate'])
            self.assertEqual(receipt['outer_measure'], TOEPLITZ_OUTER_MEASURE)
            self.assertEqual(receipt['construction'], SWAP_CONSTRUCTION)
            self.assertEqual(receipt['field_order'], 256)
            self.assertLess(F(receipt['log2_upper_rational']), -60)
        self.assertEqual(receipts[0]['log2_upper_rational'], receipts[1]['log2_upper_rational'])
        self.assertEqual([row['log2_upper_rational'] for row in receipts[0]['results']],
                         [row['log2_upper_rational'] for row in receipts[1]['results']])


if __name__ == '__main__':
    unittest.main()
