"""Tiny coverage bookkeeping and recurrence tests, without numerical gates."""

from fractions import Fraction as F
from pathlib import Path
import unittest

import numpy as np
from flint import arb,ctx

import hd_field128_lift_certificate as base
import hd_field128_lift_tail_certificate as tail
import hd_field128_lift_alias_certificate as alias
import hd_field128_lift_alias_cover as cover
from hd_shifted_coverage import RegionalAliases


def fake_bank():
    bank = object.__new__(alias.AliasBank)
    bank.grid,bank.precision = 64,256
    bank.markers = ((F(1,4),F(1,8)),(F(1,8),F(1,2)))
    bank.atlases = []
    for offset in (1,4):
        values = np.empty((4,64,64))
        for i in range(4):
            values[i] = np.arange(4096).reshape((64,64))%7+i+offset
        bank.atlases.append(RegionalAliases(values,(arb(1),arb(2),arb(3),arb(4))))
    bank.numerators = [tuple(map(arb,(20,30,40,50)))]*2
    bank.constant,bank.slope = arb(-123),arb(2)
    bank.check_pins = lambda:None
    bank.metadata = lambda:dict(directed=True,whole_code_certificate=False)
    return bank


class AliasCoverTests(unittest.TestCase):
    def test_recurrence_matches_direct_normalization(self):
        bank = fake_bank()
        cache = cover.RegionalCache(bank)
        with ctx.workprec(256):
            for weight in (1,7,66,8191,8192):
                expected = bank.regional(weight)
                actual = cache.regional(weight,base.Budget(10))
                self.assertEqual(set(expected),set(actual))
                for full in expected:
                    for a,b in zip(actual[full],expected[full]):
                        self.assertTrue(abs(a-b)<=max(a,b)*arb(2)**-230)

    def test_band_replays_all_weights_and_sums_once(self):
        bank = fake_bank()
        receipt = cover.band_receipt(bank,1,3,(F(1,8),F(1)),max_seconds=20)
        self.assertEqual([row['active_groups'] for row in receipt['weight_receipts']],[1,2,3])
        self.assertTrue(receipt['every_integer_weight_evaluated'])
        self.assertEqual(receipt['pending_weight_ranges'],[])
        self.assertFalse(receipt['whole_code_certificate'])
        with ctx.workprec(256):
            for row in receipt['weight_receipts']:
                weight = row['active_groups']
                regions = bank.regional(weight)
                exact = min(base.rational_ceiling(bank.constant+bank.slope*weight+tail.log2arb(
                            tail.support_sum(base.PAIRS,weight,regions,r,method='runs')))
                            for r in (F(1,8),F(1)))
                self.assertLessEqual(abs(F(row['log2_upper_rational'])-exact),F(1,10**11))
        self.assertEqual(receipt['log2_upper_rational'],cover.union_ceiling(receipt['weight_receipts'],256))

    def test_timeout_never_claims_uncomputed_weight(self):
        receipt = cover.band_receipt(fake_bank(),1,3,(F(1),),max_seconds=1e-12)
        self.assertFalse(receipt['complete_requested_range'])
        self.assertEqual(receipt['pending_weight_ranges'],[[1,3]])
        self.assertIsNone(receipt['log2_upper_rational'])
        self.assertEqual(receipt['failure']['error_type'],'TimeoutError')

    def test_box_partition_checks_clipping_overlap_and_gaps(self):
        root = ((3,7),(0,4096),(0,4096))
        boxes = (((3,7),(0,1),(0,4096)),((3,7),(2,4096),(0,4096)))
        checked,total = cover.checked_partition(3,7,boxes)
        self.assertEqual(checked,boxes)
        self.assertEqual(total,tail.profile_count(root))
        for wrong in ((boxes[0],), (root,boxes[0])):
            with self.assertRaises(ValueError):
                cover.checked_partition(3,7,wrong)

    def test_witness_reader_ignores_old_ceilings(self):
        old = dict(log2_upper_rational='-999999',fixed_witness=dict(pair_swaps=True,
                   orientation_aggregation='independent-fair-swaps',grid=64,
                   source_exact_witness=dict(z='3/4',alpha='1/2',potential='1/8',markers=['1/4 1/8'])))
        actual = cover.witness_from_artifact(old,old=True)
        self.assertEqual(actual,dict(z=F(3,4)**16,alpha=F(1,32),ratio=F(1,8),markers=((F(1,4),F(1,8)),),grid=64))
        old['fixed_witness']['source_exact_witness']['z'] = .75
        with self.assertRaises(ValueError):
            cover.witness_from_artifact(old,old=True)

    def test_old_manifest_only_supplies_clipped_geometry(self):
        plan = cover.manifest_seeds(Path(__file__).with_name('hd_shifted_closed_cover_manifest.json'))
        self.assertEqual(plan['seeds'][0]['weight_range'],[65,256])
        self.assertEqual(plan['seeds'][-1]['weight_range'],[5120,5853])
        self.assertFalse(plan['old_numerical_ceilings_reused'])
        self.assertEqual(sum(hi-lo+1 for row in plan['seeds'] for lo,hi in [row['weight_range']]),5853-65+1)
        self.assertIn('hd_field128_lift_acb_alias_certificate.py',cover.source_pins())
        with self.assertRaises(ValueError):
            cover.checked_range(65,1024,256)


if __name__=='__main__':
    unittest.main()
