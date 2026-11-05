"""Construction identity, unsigned-only transfer, and explicit remaining coverage."""
from copy import deepcopy
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from hd_shifted_cover_backend import CONSTRUCTION
from hd_shifted_cover_ledger import Box, Ledger
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION
from hd_shifted_swap_ledger import (SwapLedger, inherit_unsigned_ledger,
                                    endpoint_receipt, seed_ledger, replay_swap_plan,
                                    replay_sparse_family, replay_sparse_batch,replay_swap_cover)


def old_receipt(box):
    return dict(ranges=box.ranges, log2_upper_rational='-100', directed=True,
                covers_all_even_odd_splits=True, construction=CONSTRUCTION,
                field_order=256, forbidden_weight=734003,
                support_count_method='matching-profile-envelope')


class SwapLedgerTests(unittest.TestCase):
    def test_normal_import_rejects_ensemble_mismatch(self):
        root = Box(((1,6),(0,3),(0,3)))
        source, target = Ledger(root,pairs=3), SwapLedger(root,pairs=3)
        source.accept(old_receipt(Box(((4,6),(0,3),(0,3)))))
        before = target.summary()
        with self.assertRaises(ValueError):
            target.accept(source.accepted[0])
        with self.assertRaises(ValueError):
            target.import_band(source,4,6)
        self.assertEqual(target.summary(),before)

    def test_explicit_inheritance_preserves_mass_and_unresolved_profiles(self):
        source = Ledger(Box(((1,6),(0,3),(0,3))),pairs=3)
        source.accept(old_receipt(Box(((4,6),(0,3),(0,3)))))
        target = inherit_unsigned_ledger(source)
        self.assertEqual(target.summary()['unresolved_profile_triples'],
                         source.summary()['unresolved_profile_triples'])
        self.assertEqual(target.summary()['accepted_union_log2_upper'],
                         source.summary()['accepted_union_log2_upper'])
        self.assertEqual(target.accepted[0]['inherited_from'],CONSTRUCTION)
        self.assertEqual(target.accepted[0]['construction'],SWAP_CONSTRUCTION)
        self.assertEqual(source.accepted[0]['construction'],CONSTRUCTION)
        self.assertFalse(target.summary()['whole_code_certificate'])

    def test_oriented_or_unidentified_receipts_never_inherit(self):
        target = SwapLedger(Box(((1,6),(0,3),(0,3))),pairs=3)
        receipt = old_receipt(Box(((4,6),(0,3),(0,3))))
        bad = []
        item = deepcopy(receipt); item.pop('support_count_method'); bad.append(item)
        item = deepcopy(receipt); item['covers_all_even_odd_splits'] = False; bad.append(item)
        item = deepcopy(receipt); item['construction'] = SWAP_CONSTRUCTION; bad.append(item)
        for item in bad:
            with self.assertRaises(ValueError):
                target.accept_inherited(item)
        self.assertEqual(target.accepted,[])

    def test_endpoint_has_one_unsigned_profile_but_two_supports(self):
        receipt = endpoint_receipt()
        self.assertEqual(receipt['integer_profile_triples'],1)
        self.assertEqual(receipt['exact_support_count'],2)
        self.assertLess(F(receipt['log2_upper_rational']),-95000)
        ledger = SwapLedger()
        ledger.accept(receipt)
        self.assertEqual(ledger.summary()['accepted_profile_triples'],1)
        self.assertFalse(ledger.summary()['whole_code_certificate'])

    def test_seed_keeps_missing_domain_visible(self):
        ledger = seed_ledger()
        summary = ledger.summary()
        self.assertEqual(summary['construction'],SWAP_CONSTRUCTION)
        self.assertEqual(summary['accepted_boxes'],2)
        self.assertEqual(summary['complete_top_weight_range'],(7649,8192))
        self.assertGreater(summary['unresolved_profile_triples'],0)
        self.assertFalse(summary['whole_code_certificate'])

    def test_swap_proposal_requires_explicit_identity(self):
        with self.assertRaises(ValueError):
            replay_swap_plan({})
        with self.assertRaises(ValueError):
            replay_swap_plan(dict(construction=SWAP_CONSTRUCTION,witnesses=[{}]))

    def test_sparse_replay_is_complete_only_on_its_stated_band(self):
        ledger,result = replay_sparse_family(1)
        summary = result['summary']
        self.assertEqual(ledger.construction,SWAP_CONSTRUCTION)
        self.assertTrue(summary['complete_local_cover'])
        self.assertFalse(summary['whole_code_certificate'])
        self.assertEqual(summary['total_profile_triples'],1)
        self.assertLess(F(summary['accepted_union_log2_upper']),-100)
        global_ledger = SwapLedger()
        global_ledger.import_band(ledger,1,1)
        self.assertFalse(global_ledger.summary()['complete_local_cover'])
        self.assertGreater(global_ledger.summary()['unresolved_profile_triples'],0)

    def test_combined_cover_rejects_overlapping_assignments_before_replay(self):
        proposal = ('unused',{},1,64)
        with self.assertRaises(ValueError):
            replay_swap_cover((proposal,),sparse_last=64)
        with self.assertRaises(ValueError):
            replay_swap_cover((proposal,),(proposal,))
        with self.assertRaises(ValueError):
            replay_swap_cover((('unused',{},500,524),),projected=True)
        with self.assertRaises(ValueError):
            replay_swap_cover((('unused',{},7600,7800),))
        for last in (-1,65,2.5):
            with self.assertRaises(ValueError):
                replay_swap_cover(sparse_last=last)

    def test_unsigned_projected_bounds_use_explicit_inheritance(self):
        target = SwapLedger(Box(((1,6),(0,3),(0,3))),pairs=3)
        receipt = old_receipt(Box(((4,6),(0,3),(0,3))))
        receipt['support_count_method'] = 'projected-normalization-with-cyclic-run-tilt'
        target.accept_inherited(receipt)
        self.assertEqual(target.accepted[0]['inheritance'],
                         'independent-swap-unsigned-coefficient-domination')

    def test_batched_replay_claims_only_the_stated_weight_range(self):
        ledger,result = replay_sparse_batch(1,2,F(599,600),F(1),4)
        self.assertTrue(result['summary']['complete_local_cover'])
        self.assertFalse(result['summary']['whole_code_certificate'])
        self.assertEqual(ledger.root.ranges,((1,2),(0,4096),(0,4096)))
        self.assertEqual(result['receipt']['witness']['coefficient_maximum'],4)
        self.assertEqual(len(result['receipt']['weight_receipts']),2)

    def test_explicit_batch_overlap_is_rejected_before_computation(self):
        batch = (65,72,F(333,350),F(1,10),128)
        with self.assertRaises(ValueError):
            replay_swap_cover(sparse_batches=(batch,batch))
        with self.assertRaises(ValueError):
            replay_swap_cover(sparse_last=64,sparse_batches=((64,*batch[1:]),))
        with self.assertRaises(ValueError):
            replay_swap_cover(swap_bands=(('unused',{},70,80),),sparse_batches=(batch,))
        with self.assertRaises(ValueError):
            replay_swap_cover(sparse_batches=((65,72,F(1),F(1,10),128),))

    def test_potential_bands_cannot_overlap_other_proof_domains(self):
        band = ('unused',dict(requested_weight_range=[1008,1040]))
        with self.assertRaises(ValueError):
            replay_swap_cover(potential_bands=(band,band))
        with self.assertRaises(ValueError):
            replay_swap_cover(potential_bands=(band,),projected=True)
        with self.assertRaises(ValueError):
            replay_swap_cover(potential_bands=(band,),swap_bands=(('unused',{},1000,1100),))

    def test_potential_band_is_recomputed_before_acceptance(self):
        proposal = dict(requested_weight_range=[512,513],log2_upper_rational='-10000')
        receipt = dict(ranges=((512,513),(0,4096),(0,4096)),log2_upper_rational='-120',
                       directed=True,covers_all_even_odd_splits=True,construction=SWAP_CONSTRUCTION,
                       field_order=256,forbidden_weight=734003)
        with patch('hd_shifted_potential_band.replay_saved_band',return_value=receipt) as replay:
            ledger,result = replay_swap_cover(potential_bands=(('test-only',proposal),))
            replay.assert_called_once_with(proposal,precision=256,progress=None,max_weights=33)
            self.assertIn(receipt,ledger.accepted)
            self.assertFalse(result['summary']['whole_code_certificate'])
        receipt['log2_upper_rational'] = '-20'
        with patch('hd_shifted_potential_band.replay_saved_band',return_value=receipt):
            with self.assertRaises(ArithmeticError):
                replay_swap_cover(potential_bands=(('test-only',proposal),))

    def test_potential_resource_budget_is_not_supplied_by_the_proposal(self):
        proposal = dict(requested_weight_range=[400,600],max_weights=201)
        with self.assertRaises(ValueError):
            replay_swap_cover(potential_bands=(('unused',proposal),))
        for budget in (0,8193,2.5):
            with self.assertRaises(ValueError):
                replay_swap_cover(potential_max_weights=budget)


if __name__ == '__main__':
    unittest.main()
