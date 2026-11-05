from types import SimpleNamespace
from pathlib import Path
import json
import unittest
from unittest.mock import patch

from hd_shifted_bulk_search import AdaptiveAlphaProposer, PositiveWitness
from hd_shifted_cover_backend import profile_count
from hd_shifted_swap_bulk import SwapProposer, restore_identity, subdivide
from hd_shifted_swap_screen import CONSTRUCTION, SwapPositiveWitness


class SwapBulkTests(unittest.TestCase):
    def test_default_provider_is_unchanged(self):
        self.assertIs(AdaptiveAlphaProposer.witness_type, PositiveWitness)
        self.assertIs(SwapProposer.witness_type, SwapPositiveWitness)

    def test_frontier_initial_bank_uses_swap_provider_and_budget(self):
        row = dict(profile=[4096, 1024, 1024],
                   positive=dict(mark=.04, alpha=.25, marker_log_odds=[[2., 3.]]*2))
        proposer = SwapProposer((row,), max_evaluations=3)
        with patch.object(SwapPositiveWitness, 'create', return_value='swap') as create:
            self.assertEqual(proposer.initial_bank(), ['swap']*3)
            self.assertEqual(create.call_count, 3)
        self.assertTrue(proposer.exhausted)

    def test_optimization_objective_also_uses_swap_provider(self):
        proposer = SwapProposer(max_evaluations=8, marker_maxfev=2, weight_maxfev=0)
        witness = SimpleNamespace(evaluate=lambda box: (-1000., (4096, 1024, 1024)))
        with patch.object(SwapPositiveWitness, 'create', return_value=witness) as create:
            actual = proposer.propose(((4096, 4097), (1024, 1025), (1024, 1025)), 4096, -96.)
        self.assertIs(actual, witness)
        self.assertEqual(create.call_count, 1)
        self.assertEqual(proposer.evaluations, 1)

    def test_compact_identity_is_restored_per_witness(self):
        plan = restore_identity(dict(witnesses=[dict(mark=.5, alpha=.1)],
                                     directed=False, whole_code_certificate=False))
        self.assertEqual(plan['construction'], CONSTRUCTION)
        self.assertEqual(plan['witnesses'][0]['construction'], CONSTRUCTION)
        self.assertTrue(plan['witnesses'][0]['pair_swaps'])
        self.assertFalse(plan['directed'])
        self.assertFalse(plan['whole_code_certificate'])

    def test_optional_saved_partition_is_complete_and_disjoint(self):
        path = Path(__file__).with_name('hd_shifted_swap_bulk_plan.json')
        if not path.exists():
            self.skipTest('Optional ignored fixture absent; regenerate explicitly with hd_shifted_swap_bulk.py')
        with path.open(encoding='utf-8') as source:
            plan = json.load(source)
        self.assertEqual(plan['construction'], CONSTRUCTION)
        self.assertFalse(plan['directed'])
        self.assertFalse(plan['whole_code_certificate'])
        boxes = [row['ranges'] for row in plan['accepted_boxes']+plan['unresolved_boxes']]
        self.assertEqual(sum(profile_count(tuple(map(tuple, box))) for box in boxes),
                         plan['feasible_domain_triples'])
        for index, box in enumerate(boxes):
            for other in boxes[index+1:]:
                self.assertTrue(any(hi < other_lo or other_hi < lo
                                    for (lo, hi), (other_lo, other_hi) in zip(box, other)))
        self.assertEqual(sum(row['integer_triples'] for row in plan['accepted_boxes']),
                         plan['accepted_integer_triples'])
        self.assertTrue(all(row['construction'] == CONSTRUCTION and row['pair_swaps']
                            for row in plan['witnesses']))
        self.assertLessEqual(plan['adaptive_proposer']['evaluations'],
                             plan['adaptive_proposer']['max_evaluations'])

    def test_bounded_subdivision_retains_a_complete_small_partition(self):
        saved = dict(construction=CONSTRUCTION, requested_weight_range=[1, 2],
                     witnesses=[dict(mark=.5, alpha=.1, marker_log_odds=[[0., 0.]]*2,
                                     proposal_centre=[1, 0, 0])], accepted_boxes=[],
                     unresolved_boxes=[dict(ranges=[[1, 2], [0, 1], [0, 1]])])
        witness = SimpleNamespace(serializable=lambda: saved['witnesses'][0])
        def fake_best(bank, box):
            score = -200. if all(lo == hi for lo, hi in box) else 200.
            return 0, score, tuple(lo for lo, _ in box)
        with patch.object(SwapPositiveWitness, 'create', return_value=witness), \
                patch('hd_shifted_swap_bulk.best_witness', side_effect=fake_best):
            full = subdivide(saved, first=1, last=2, max_splits=20, max_candidates=0)
            partial = subdivide(saved, first=1, last=2, max_splits=1, max_candidates=0)
        self.assertFalse(full['unresolved_boxes'])
        self.assertEqual(full['accepted_integer_triples'], full['feasible_domain_triples'])
        self.assertTrue(partial['unresolved_boxes'])
        for plan in (full, partial):
            boxes = [tuple(map(tuple, row['ranges']))
                     for row in plan['accepted_boxes']+plan['unresolved_boxes']]
            self.assertEqual(sum(map(profile_count, boxes)), plan['feasible_domain_triples'])
        self.assertEqual(partial['subdivision']['splits'], 1)
        self.assertEqual(partial['subdivision']['candidates'], 0)

    def test_optional_complete_band_fixture_counts_and_scope(self):
        path = Path(__file__).with_name('hd_shifted_swap_bulk_plan3.json')
        if not path.exists():
            self.skipTest('Optional generated full-band fixture absent; no search runs during unit tests')
        with path.open(encoding='utf-8') as source:
            plan = json.load(source)
        self.assertEqual(plan['requested_weight_range'], [4096, 4351])
        self.assertEqual(plan['accepted_integer_triples'], 1009047040)
        self.assertEqual(plan['accepted_integer_triples'], plan['feasible_domain_triples'])
        self.assertFalse(plan['unresolved_boxes'])
        self.assertFalse(plan['directed'])
        self.assertFalse(plan['whole_code_certificate'])
        self.assertEqual(sum(profile_count(tuple(map(tuple, row['ranges'])))
                             for row in plan['accepted_boxes']), 1009047040)
        self.assertLess(plan['union_log2_upper_screen'], -96.)

    def test_extra_bank_imports_witnesses_not_partition_and_deduplicates(self):
        row = dict(mark=.5, alpha=.1, marker_log_odds=[[0., 0.]]*2,
                   proposal_centre=[1, 0, 0])
        saved = restore_identity(dict(requested_weight_range=[1, 1], witnesses=[row],
            accepted_boxes=[], unresolved_boxes=[dict(ranges=[[1, 1], [0, 0], [0, 0]])]))
        extra = restore_identity(dict(witnesses=[row, dict(row, mark=.25)],
            requested_weight_range=[4096, 4351], accepted_boxes=[], unresolved_boxes=[]))
        def create(mark, alpha, odds, centre):
            return SimpleNamespace(serializable=lambda: dict(row, mark=mark))
        with patch.object(SwapPositiveWitness, 'create', side_effect=create) as provider, \
                patch('hd_shifted_swap_bulk.best_witness', return_value=(1, -200., (1, 0, 0))):
            plan = subdivide(saved, bank_plans=[extra], first=1, last=1,
                             max_candidates=0, max_splits=0)
        self.assertEqual(provider.call_count, 2)
        self.assertEqual(plan['subdivision']['initial_bank_size'], 2)
        self.assertEqual(plan['accepted_integer_triples'], 1)
        self.assertEqual(plan['requested_weight_range'], [1, 1])
        self.assertEqual(plan['witnesses'][0]['mark'], .25)
        extra['witnesses'][0]['pair_swaps'] = False
        with self.assertRaisesRegex(ValueError, 'independent fair swaps'):
            subdivide(saved, bank_plans=[extra], first=1, last=1)

    def test_optional_plan4_is_a_complete_disjoint_bounded_band(self):
        path = Path(__file__).with_name('hd_shifted_swap_bulk_plan4.json')
        if not path.exists():
            self.skipTest('Optional generated plan4 absent; no search runs during unit tests')
        with path.open(encoding='utf-8') as source:
            plan = json.load(source)
        self.assertEqual(plan['requested_weight_range'], [4352, 5119])
        self.assertEqual(plan['accepted_integer_triples'], 2305327616)
        self.assertEqual(plan['accepted_integer_triples'], plan['feasible_domain_triples'])
        self.assertFalse(plan['unresolved_boxes'])
        self.assertFalse(plan['directed'])
        self.assertFalse(plan['whole_code_certificate'])
        boxes = [tuple(map(tuple, row['ranges'])) for row in plan['accepted_boxes']]
        self.assertEqual(sum(map(profile_count, boxes)), plan['feasible_domain_triples'])
        for index, box in enumerate(boxes):
            self.assertTrue(4352 <= box[0][0] <= box[0][1] <= 5119)
            for other in boxes[index+1:]:
                self.assertTrue(any(hi < other_lo or other_hi < lo
                                    for (lo, hi), (other_lo, other_hi) in zip(box, other)))
        stats = plan['subdivision']
        self.assertTrue(stats['whole_box_centering'])
        self.assertLessEqual(stats['splits'], stats['max_splits'])
        self.assertLessEqual(stats['candidates'], stats['max_candidates'])
        self.assertEqual(stats['local_optimizations'], 0)
        self.assertLess(plan['union_log2_upper_screen'], -96.)
        self.assertTrue(all(row['construction'] == CONSTRUCTION and row['pair_swaps']
                            and row['orientation_aggregation'] == 'independent-fair-swaps'
                            for row in plan['witnesses']))

    def test_whole_box_centering_does_not_need_a_point_probe(self):
        row = dict(mark=.5, alpha=.1, marker_log_odds=[[0., 0.]]*2,
                   proposal_centre=[1, 0, 0])
        saved = restore_identity(dict(requested_weight_range=[1, 2], witnesses=[row],
            accepted_boxes=[], unresolved_boxes=[dict(ranges=[[1, 2], [0, 0], [0, 0]])]))
        witness = SimpleNamespace(z=.5, alpha=.1, odds=((0., 0.),)*2,
            centre=(1, 0, 0), serializable=lambda: row,
            evaluate=lambda box: (-200., (1, 0, 0)))
        with patch.object(SwapPositiveWitness, 'create', return_value=witness), \
                patch('hd_shifted_swap_bulk.best_witness', return_value=(0, 200., (1, 0, 0))) as best:
            plan = subdivide(saved, first=1, last=2, centre_boxes=True,
                             max_splits=0, max_candidates=1)
        self.assertEqual(best.call_count, 1)
        self.assertEqual(plan['subdivision']['candidates'], 1)
        self.assertFalse(plan['unresolved_boxes'])
        self.assertEqual(plan['accepted_integer_triples'], 2)

    def test_optional_plan5_retains_every_gap_with_bounded_work(self):
        path = Path(__file__).with_name('hd_shifted_swap_bulk_plan5.json')
        if not path.exists():
            self.skipTest('Optional generated plan5 absent; no search runs during unit tests')
        with path.open(encoding='utf-8') as source:
            plan = json.load(source)
        self.assertEqual(plan['requested_weight_range'], [3072, 4095])
        self.assertEqual(plan['construction'], CONSTRUCTION)
        self.assertFalse(plan['directed'])
        self.assertFalse(plan['whole_code_certificate'])
        accepted = sum(profile_count(tuple(map(tuple, row['ranges'])))
                       for row in plan['accepted_boxes'])
        unresolved = sum(profile_count(tuple(map(tuple, row['ranges'])))
                         for row in plan['unresolved_boxes'])
        self.assertEqual(accepted, 2409216130)
        self.assertEqual(unresolved, 903322494)
        self.assertEqual(accepted, plan['accepted_integer_triples'])
        self.assertEqual(accepted+unresolved, plan['feasible_domain_triples'])
        stats = plan['subdivision']
        self.assertTrue(stats['whole_box_centering'])
        self.assertLessEqual(stats['splits'], stats['max_splits'])
        self.assertLessEqual(stats['candidates'], stats['max_candidates'])
        self.assertEqual(stats['local_optimizations'], 0)


if __name__ == '__main__':
    unittest.main()
