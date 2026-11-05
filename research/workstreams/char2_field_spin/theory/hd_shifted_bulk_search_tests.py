from itertools import product
import json
from math import isfinite, log2
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from hd_shifted_bulk_search import (
    ALPHA_CHECKPOINTS, AdaptiveAlphaProposer, PositiveWitness, alpha_checkpoint_bounds,
    best_witness, box_geometry, initial_partition, split_box, volume,
)
from hd_shifted_unoriented_screen import global_box_bound
from hd_shifted_cover_backend import profile_count


class BulkSearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Plans are local replay artifacts excluded by the repository's JSON
        # policy. Their integration tests skip on a clean checkout; bounded
        # regeneration is an explicit command, never an implicit unit-test job.
        first = Path(__file__).with_name('hd_shifted_bulk_plan.json')
        second = Path(__file__).with_name('hd_shifted_bulk_plan2.json')
        third = Path(__file__).with_name('hd_shifted_bulk_plan3.json')
        fourth = Path(__file__).with_name('hd_shifted_bulk_plan4.json')
        cls.plan = json.loads(first.read_text()) if first.exists() else None
        cls.continuation = json.loads(second.read_text()) if second.exists() else None
        cls.lower_continuation = json.loads(third.read_text()) if third.exists() else None
        cls.adaptive_continuation = json.loads(fourth.read_text()) if fourth.exists() else None

    def test_alpha_correction_closes_fixed_boxes(self):
        values = alpha_checkpoint_bounds()
        self.assertEqual(len(values), 6)
        for row, (_, _, expected) in zip(values, ALPHA_CHECKPOINTS):
            self.assertLess(row['union_log2_upper_screen'], -40)
            self.assertAlmostEqual(row['union_log2_upper_screen'], expected, places=5)
            self.assertFalse(row['directed'])

    def test_adaptive_alpha_proposer_honors_hard_evaluation_budget(self):
        box = ((5120, 5127), (1152, 1167), (1408, 1423))
        witness = PositiveWitness(.5, .1, ((0., 0.),)*2, (0, 0, 0), 0., (0., 0., 0.), 0.)
        proposer = AdaptiveAlphaProposer(max_evaluations=5, marker_maxfev=20, weight_maxfev=10)
        with patch.object(PositiveWitness, 'create', return_value=witness) as create:
            with patch.object(PositiveWitness, 'evaluate', return_value=(100., (0., 0., 0.))):
                result = proposer.propose(box, 5247.5, -96.)
        self.assertIs(result, witness)
        self.assertEqual(proposer.evaluations, 5)
        self.assertEqual(create.call_count, 5)
        self.assertTrue(proposer.exhausted)

    def test_adaptive_alpha_stops_after_a_sufficient_whole_box_witness(self):
        box = ((5120, 5127), (1152, 1167), (1408, 1423))
        witness = PositiveWitness(.5, .1, ((0., 0.),)*2, (0, 0, 0), 0., (0., 0., 0.), 0.)
        proposer = AdaptiveAlphaProposer(max_evaluations=50)
        with patch.object(PositiveWitness, 'create', return_value=witness):
            with patch.object(PositiveWitness, 'evaluate', return_value=(-100., (0., 0., 0.))):
                self.assertIs(proposer.propose(box, 5247.5, -96.), witness)
        self.assertEqual(proposer.evaluations, 1)
        self.assertEqual(proposer.weight_optimizations, 0)
        self.assertTrue(set((.25, .5, .75, 1.)).issubset(proposer.alpha_family))
        self.assertNotEqual(*AdaptiveAlphaProposer.marker_seed(box, .5)[0].reshape(2, 2).tolist())

    def test_cached_witness_matches_direct_global_box(self):
        z, alpha = .053587729568108665, .11372673628845778
        odds = ((-.31797964871576334, .9361653072364948),)*2
        box = ((4032, 4160), (1536, 1664), (1536, 1664))
        witness = PositiveWitness.create(z, alpha, odds)
        expected = global_box_bound(box, z, alpha, odds)
        actual, worst = witness.evaluate(box)
        self.assertAlmostEqual(actual, expected['union_log2_upper_screen'], places=6)
        self.assertEqual(worst, tuple(expected['worst_corner']))
        # Evaluating boxes must not rebuild their regional matrices.
        with patch('hd_shifted_bulk_search.positive_pair_matrix', side_effect=AssertionError('rebuild')):
            self.assertEqual(witness.evaluate(box), (actual, worst))

    def test_best_witness_uses_the_same_witness_for_all_corners(self):
        box = ((4095, 4097), (1599, 1601), (1599, 1601))
        first = PositiveWitness(.5, .1, ((0., 0.),)*2, (0, 0, 0), 0., (2., -1., 3.), 0.)
        second = PositiveWitness(.5, .1, ((0., 0.),)*2, (0, 0, 0), 15., (-2., 1., -3.), 0.)
        expected = [w.evaluate(box)[0] for w in (first, second)]
        index, score, _ = best_witness((first, second), box)
        self.assertEqual(index, int(np.argmin(expected)))
        self.assertAlmostEqual(score, min(expected))

    def test_initial_partition_is_disjoint_and_contains_feasible_small_domain(self):
        pairs = 8
        interiors, boundaries = initial_partition(8, 14, slab_width=3, full_width=2, pairs=pairs)
        boxes = [box for box, _ in interiors]+boundaries
        seen = set()
        for box in boxes:
            for point in product(*(range(lo, hi+1) for lo, hi in box)):
                self.assertNotIn(point, seen)
                seen.add(point)
        for r in range(8, 15):
            for a in range(r-pairs, r//2+1):
                for b in range(r-pairs, r//2+1):
                    self.assertIn((r, a, b), seen)
        for box, _ in interiors:
            for r, a, b in product(*(tuple(dict.fromkeys(interval)) for interval in box)):
                self.assertGreaterEqual(min(pairs-r+a, r-2*a, a, pairs-r+b, r-2*b, b), 0)

    def test_integer_bisection_preserves_disjoint_partition(self):
        box = ((4096, 4099), (1000, 1004), (1200, 1201))
        left, right = split_box(box)
        self.assertEqual(volume(left)+volume(right), volume(box))
        expanded = lambda value: set(product(*(range(lo, hi+1) for lo, hi in value)))
        self.assertFalse(expanded(left) & expanded(right))
        self.assertEqual(expanded(left) | expanded(right), expanded(box))

    def test_invalid_corner_box_is_clipped(self):
        vertices, entropy = box_geometry(((4096, 4200), (0, 100), (0, 100)))
        self.assertEqual(entropy.shape, (2, len(vertices)))
        for r, a, b in vertices:
            self.assertGreaterEqual(min(4096-r+a, r-2*a, a, 4096-r+b, r-2*b, b), 0)
        self.assertIn((4096., 0., 0.), tuple(map(tuple, vertices)))

    def test_empty_histogram_domain_is_rejected(self):
        with self.assertRaises(ValueError):
            box_geometry(((4200, 4300), (0, 100), (0, 100)))

    def test_clipped_envelope_contains_integer_interior(self):
        box = ((4096, 4100), (0, 4), (0, 4))
        witness = PositiveWitness(.5, .1, ((0., 0.),)*2, (0, 0, 0), 0., (2., -1., 3.), 0.)
        bound, _ = witness.evaluate(box)
        one_profile_bound = bound-log2(profile_count(box))
        for point in product(*(range(lo, hi+1) for lo, hi in box)):
            singleton = tuple((n, n) for n in point)
            if profile_count(singleton):
                value, _ = witness.evaluate(singleton)
                self.assertLessEqual(value, one_profile_bound+1e-7)

    def test_saved_plan_is_a_disjoint_complete_partition(self):
        if self.plan is None:
            self.skipTest('Optional bulk plan absent; use the documented bounded-search command')
        result = self.plan
        accepted = result['accepted_boxes']
        leaves = accepted+result['unresolved_boxes']
        boxes = [tuple(map(tuple, row['ranges'])) for row in leaves]
        first, last = result['requested_weight_range']
        for i, box in enumerate(boxes):
            self.assertTrue(first <= box[0][0] <= box[0][1] <= last)
            self.assertTrue(all(0 <= lo <= hi <= 4096 for lo, hi in box[1:]))
            for previous in boxes[:i]:
                self.assertTrue(any(hi < otherlo or otherhi < lo
                                    for (lo, hi), (otherlo, otherhi) in zip(box, previous)))
        self.assertEqual(sum(profile_count(box) for box in boxes), result['feasible_domain_triples'])
        self.assertEqual(sum(row['integer_triples'] for row in accepted), result['accepted_integer_triples'])
        for row in accepted:
            self.assertEqual(profile_count(tuple(map(tuple, row['ranges']))), row['integer_triples'])
            self.assertTrue(0 <= row['witness_id'] < len(result['witnesses']))
            self.assertTrue(isfinite(row['union_log2_upper_screen']))
            self.assertLessEqual(row['union_log2_upper_screen'], result['threshold'])
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])

    def test_saved_dense_band_has_no_unresolved_feasible_profile(self):
        if self.plan is None:
            self.skipTest('Optional bulk plan absent; use the documented bounded-search command')
        result = self.plan
        for row in result['unresolved_boxes']:
            self.assertLess(row['ranges'][0][1], 7168)
        count = 0
        for row in result['accepted_boxes']:
            (lo, hi), a, b = row['ranges']
            if hi >= 7168:
                count += profile_count(((max(lo, 7168), hi), tuple(a), tuple(b)))
        self.assertEqual(count, 76513888)
        self.assertEqual(count, profile_count(((7168, 7648), (0, 4096), (0, 4096))))

    def test_saved_witness_values_replay(self):
        if self.plan is None:
            self.skipTest('Optional bulk plan absent; use the documented bounded-search command')
        result = self.plan
        for index in (0, len(result['accepted_boxes'])//2, len(result['accepted_boxes'])-1):
            row = result['accepted_boxes'][index]
            saved = result['witnesses'][row['witness_id']]
            witness = PositiveWitness.create(saved['mark'], saved['alpha'], saved['marker_log_odds'])
            actual, _ = witness.evaluate(tuple(map(tuple, row['ranges'])))
            self.assertAlmostEqual(actual, row['union_log2_upper_screen'], places=6)

    def test_continuation_closes_exact_remaining_dense_domain(self):
        if self.plan is None or self.continuation is None:
            self.skipTest('Optional continuation plan absent; use the documented bounded-search command')
        source, result = self.plan, self.continuation
        self.assertEqual(result['requested_weight_range'], [6144, 7167])
        self.assertEqual(result['proposals_used'], 17)
        self.assertFalse(result['unresolved_boxes'])
        accepted = result['accepted_boxes']
        new = [row for row in accepted if not row['inherited']]
        self.assertEqual(len(new), 19)
        self.assertEqual(sum(row['integer_triples'] for row in new), 95795009)
        self.assertEqual(sum(row['integer_triples'] for row in accepted), 627922944)
        self.assertEqual(result['accepted_integer_triples'], result['feasible_domain_triples'])
        for i, row in enumerate(accepted):
            box = tuple(map(tuple, row['ranges']))
            self.assertEqual(profile_count(box), row['integer_triples'])
            for previous in accepted[:i]:
                self.assertTrue(any(hi < otherlo or otherhi < lo
                                    for (lo, hi), (otherlo, otherhi) in zip(box, previous['ranges'])))
        for row in new:
            # Every newly accepted rectangle lies in a previous unresolved leaf.
            self.assertTrue(any(all(oldlo <= lo <= hi <= oldhi
                                    for (lo, hi), (oldlo, oldhi) in zip(row['ranges'], old['ranges']))
                                for old in source['unresolved_boxes']))
        for row in (new[0], new[-1]):
            saved = result['witnesses'][row['witness_id']]
            witness = PositiveWitness.create(saved['mark'], saved['alpha'], saved['marker_log_odds'])
            actual, _ = witness.evaluate(tuple(map(tuple, row['ranges'])))
            self.assertAlmostEqual(actual, row['union_log2_upper_screen'], places=6)

    def test_lower_continuation_preserves_partition_and_budget(self):
        if self.plan is None or self.lower_continuation is None:
            self.skipTest('Optional lower continuation absent; use the documented bounded-search command')
        result = self.lower_continuation
        self.assertEqual(result['proposals_used'], 400)
        self.assertEqual(result['requested_weight_range'], [5120, 6143])
        leaves = result['accepted_boxes']+result['unresolved_boxes']
        boxes = [tuple(map(tuple, row['ranges'])) for row in leaves]
        for i, box in enumerate(boxes):
            for previous in boxes[:i]:
                self.assertTrue(any(hi < otherlo or otherhi < lo
                                    for (lo, hi), (otherlo, otherhi) in zip(box, previous)))
        self.assertEqual(sum(profile_count(box) for box in boxes), 1702713344)
        self.assertEqual(result['accepted_integer_triples'], 1458395260)
        new = [row for row in result['accepted_boxes'] if not row['inherited']]
        self.assertEqual(len(new), 279)
        self.assertEqual(sum(row['integer_triples'] for row in new), 712394880)
        self.assertEqual(max(row['ranges'][0][1] for row in result['unresolved_boxes']), 6015)
        for row in new:
            self.assertTrue(any(all(oldlo <= lo <= hi <= oldhi
                                    for (lo, hi), (oldlo, oldhi) in zip(row['ranges'], old['ranges']))
                                for old in self.plan['unresolved_boxes']))

    def test_adaptive_continuation_closes_the_full_band(self):
        if self.lower_continuation is None or self.adaptive_continuation is None:
            self.skipTest('Optional adaptive plan absent; use the documented bounded-search command')
        result = self.adaptive_continuation
        self.assertFalse(result['unresolved_boxes'])
        self.assertEqual(result['accepted_integer_triples'], 1702713344)
        self.assertEqual(result['accepted_integer_triples'], result['feasible_domain_triples'])
        self.assertEqual(result['proposals_used'], 46)
        stats = result['adaptive_proposer']
        self.assertEqual(stats['evaluations'], 9153)
        self.assertLessEqual(stats['evaluations'], stats['max_evaluations'])
        boxes = [tuple(map(tuple, row['ranges'])) for row in result['accepted_boxes']]
        for i, box in enumerate(boxes):
            for previous in boxes[:i]:
                self.assertTrue(any(hi < otherlo or otherhi < lo
                                    for (lo, hi), (otherlo, otherhi) in zip(box, previous)))
        self.assertEqual(sum(profile_count(box) for box in boxes), result['feasible_domain_triples'])
        new = [row for row in result['accepted_boxes'] if not row['inherited']]
        self.assertEqual(len(new), 159)
        self.assertEqual(sum(row['integer_triples'] for row in new), 244318084)
        for row in new:
            self.assertTrue(any(all(oldlo <= lo <= hi <= oldhi
                                    for (lo, hi), (oldlo, oldhi) in zip(row['ranges'], old['ranges']))
                                for old in self.lower_continuation['unresolved_boxes']))
        parameters = {(result['witnesses'][row['witness_id']]['mark'],
                       result['witnesses'][row['witness_id']]['alpha']) for row in new}
        self.assertEqual(len(parameters), 2)


if __name__ == '__main__':
    unittest.main()
