import copy
from fractions import Fraction as Q
import unittest

import search_strategy as strategy


class SearchStrategyTests(unittest.TestCase):
    def test_preserves_accepted_and_input(self):
        original = dict(leaves={'0': dict(witness={'parameters': ['1', '2', '3']})},
            unresolved={'1': dict(cell=['1/2', '1'], witness={'variance_partition': ['old']}, proposal=999)},
            visited=9)
        before = copy.deepcopy(original)
        result, stats = strategy.presplit((0, 1), original, Q(1, 8), 3)
        self.assertEqual(original, before)
        self.assertEqual(result['leaves'], original['leaves'])
        self.assertEqual(result['visited'], 9)
        self.assertEqual(set(result['unresolved']), {'100', '101', '110', '111'})
        self.assertTrue(all('witness' not in row and 'proposal' not in row for row in result['unresolved'].values()))
        self.assertEqual(stats['skipped_internal_cells'], 3)
        self.assertEqual(stats['numerical_evaluations'], 0)

    def test_retains_unsplit_candidate_exactly(self):
        source = dict(leaves={'0': {}}, unresolved={'1': {'witness': {'tilt': '1'}}})
        result, stats = strategy.presplit((0, 1), source, Q(1, 2))
        self.assertEqual(result, source)
        self.assertIsNot(result['unresolved']['1'], source['unresolved']['1'])
        self.assertEqual(stats['skipped_internal_cells'], 0)

    def test_nonzero_rational_root(self):
        root = (Q(1683, 7936000), Q(1))
        original = dict(leaves={'1': {}}, unresolved={'0': {}})
        result, stats = strategy.presplit(root, original, Q(1, 1024), 10)
        self.assertEqual(stats['pending_after'], 512)
        self.assertEqual(stats['skipped_internal_cells'], 511)
        cells = strategy.partition(root, result['leaves'], result['unresolved'])
        self.assertTrue(all(cells[p][1]-cells[p][0] <= Q(1, 1024) for p in result['unresolved']))

    def test_bad_partitions_rejected(self):
        for leaves, unresolved in (({'0': {}}, {}), ({'0': {}}, {'00': {}, '1': {}}),
                ({'0': {}}, {'0': {}, '1': {}}), ({'x': {}}, {}),
                ({'0': {}}, {'1': {'cell': ['0', '1']}})):
            self.assertRaises(ValueError, strategy.partition, (0, 1), leaves, unresolved)

    def test_depth_limit_fails_without_mutation(self):
        source = dict(leaves={}, unresolved={'': {}}, visited=0)
        before = copy.deepcopy(source)
        self.assertRaises(ValueError, strategy.presplit, (0, 1), source, Q(1, 1024), 9)
        self.assertEqual(source, before)
        self.assertRaises(ValueError, strategy.presplit, (0, 1), source, Q(1, 2), True)

    def test_complete_cover_is_unchanged(self):
        source = dict(leaves={'': {'witness': {'tilt': '1'}}}, unresolved={})
        result, stats = strategy.presplit((0, 1), source)
        self.assertEqual(result, source)
        self.assertEqual(stats['pending_after'], 0)

    def test_retile_coalesces_only_unresolved_siblings(self):
        original = dict(leaves={'00': {'witness': {'important': [1, 2]}}},
            unresolved={p: {'witness': {'old': p}, 'proposal': 999, 'upper': [1, -999]}
                        for p in ('01', '10', '11')}, visited=55, extra={'retained': [1]})
        before = copy.deepcopy(original)
        result, stats = strategy.retile((0, 1), original, Q(1, 2), 4)
        self.assertEqual(original, before)
        self.assertEqual(result['leaves'], original['leaves'])
        self.assertEqual(result['visited'], 55)
        self.assertEqual(result['extra'], original['extra'])
        self.assertEqual(set(result['unresolved']), {'01', '1'})
        self.assertEqual(result['unresolved']['01'], original['unresolved']['01'])
        self.assertEqual(result['unresolved']['1'], {'cell': ['1/2', '1']})
        self.assertEqual(stats['pending_before'], 3)
        self.assertEqual(stats['pending_after_coalesce'], 2)
        self.assertEqual(stats['coalesced_internal_cells'], 1)
        self.assertEqual(stats['presplit_internal_cells'], 0)

    def test_merge_then_split_never_restores_old_witnesses(self):
        original = dict(leaves={}, unresolved={format(i, '03b'): {'witness': {'i': i}, 'proposal': -100}
                                               for i in range(8)}, visited=8)
        result, stats = strategy.retile((0, 1), original, Q(1, 8), 3)
        self.assertEqual(set(result['unresolved']), set(original['unresolved']))
        self.assertTrue(all('witness' not in row and 'proposal' not in row
                            for row in result['unresolved'].values()))
        self.assertEqual(stats['coalesced_internal_cells'], 7)
        self.assertEqual(stats['pending_after_coalesce'], 1)
        self.assertEqual(stats['presplit_internal_cells'], 7)
        self.assertEqual(stats['numerical_evaluations'], 0)

    def test_retile_all_small_accepted_patterns(self):
        root = (Q(1, 17), Q(9, 10))
        for mask in range(256):
            leaves, pending = {}, {}
            for i in range(8):
                target = leaves if mask & (1 << i) else pending
                target[format(i, '03b')] = {'witness': {'old': i}}
            original = dict(leaves=leaves, unresolved=pending, visited=12)
            before = copy.deepcopy(original)
            result, stats = strategy.retile(root, original, Q(1, 4), 5)
            cells = strategy.partition(root, result['leaves'], result['unresolved'])
            self.assertEqual(original, before)
            self.assertEqual(result['leaves'], leaves)
            self.assertEqual(result['visited'], 12)
            self.assertTrue(all(cells[p][1]-cells[p][0] <= Q(1, 4) for p in result['unresolved']))
            self.assertEqual(stats['numerical_evaluations'], 0)

    def test_retile_bad_partition_or_depth_preserves_input(self):
        original = dict(leaves={}, unresolved={'0': {}, '1': {}}, visited=1)
        before = copy.deepcopy(original)
        self.assertRaises(ValueError, strategy.retile, (0, 1), original, Q(1, 1024), 9)
        self.assertEqual(original, before)
        self.assertRaises(ValueError, strategy.retile, (0, 1), original, 0, 9)
        self.assertRaises(ValueError, strategy.retile, (0, 1), original, Q(1, 2), True)
        invalid = dict(leaves={'0': {}}, unresolved={'00': {}, '1': {}})
        self.assertRaises(ValueError, strategy.retile, (0, 1), invalid, Q(1, 2), 9)

    def test_retile_complete_cover_is_unchanged(self):
        source = dict(leaves={'': {'witness': {'tilt': '1'}}}, unresolved={}, visited=14)
        result, stats = strategy.retile((0, 1), source)
        self.assertEqual(result, source)
        self.assertEqual(stats['coalesced_internal_cells'], 0)
        self.assertEqual(stats['pending_after'], 0)


if __name__ == '__main__':
    unittest.main()
