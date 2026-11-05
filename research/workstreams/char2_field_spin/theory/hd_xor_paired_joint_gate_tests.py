"""Exact source-geometry integration tests; synthetic caps are not certificates."""
from fractions import Fraction as F
from math import log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_xor_paired_joint_gate as gate


class JointGateTests(unittest.TestCase):
    def test_schedule_partition(self):
        _, complete = gate.checked_schedule(gate.SCHEDULE)
        self.assertTrue(complete)
        self.assertFalse(gate.checked_schedule(gate.SCHEDULE[:1])[1])
        for schedule in (gate.SCHEDULE+gate.SCHEDULE[:1],
                         [(0, 17, '1/2', '-1')],
                         [(1, 17, '1', '-1')],
                         [(1, 17, '1/2', '1')]):
            with self.assertRaises(ValueError):
                gate.checked_schedule(schedule)

    def test_saved_witnesses_keep_source_identity(self):
        rows = [dict(source_matching_counts=list(source), mu_range=[1, 17],
                     mark=str(z), log2_potential='-2', failure_log2_screen=-999)
                for source, z in (((0, 0), F(999, 1000)), ((8, 8), F(998, 1000)))]
        schedules, complete = gate.literal_schedules({(0, 0), (8, 8)}, witness_rows=rows)
        self.assertTrue(complete)
        self.assertEqual(schedules[0, 0][0][2], F(999, 1000))
        self.assertEqual(schedules[8, 8][0][2], F(998, 1000))
        with self.assertRaises(ValueError):
            gate.literal_schedules({(0, 0), (8, 8)}, witness_rows=rows[:1])
        with self.assertRaises(ValueError):
            gate.literal_schedules({(0, 0)}, witness_rows=rows+rows[:1])

    def test_order_and_subset_lower_bound(self):
        raw = {(m, f): np.array([[f, -f], [m, -m]], dtype=float)
               for m in range(18) for f in range(m//2+1)}
        tables = gate.ordered_tables(raw, 0, 8)
        self.assertEqual(tables.shape, (8, 18, 2, 2))
        for position in range(8):
            full = 8 if position % 2 else 0
            for m in range(18):
                low, high = max(0, m-17+full), min(full, m//2)
                np.testing.assert_array_equal(tables[position, m],
                    np.array([[high, -low], [m, -m]], dtype=float))
        # At full source occupancy the exact matching count cannot be dropped.
        self.assertEqual(tables[1, 17, 0, 1], -8)

    def test_arb_order_and_upper_endpoints(self):
        from flint import arb
        raw = {(m, f): tuple(arb(v) for v in (f, 20-f, m, 20-m))
               for m in range(18) for f in range(m//2+1)}
        tables = gate.ordered_tables(raw, 0, 8, directed=True)
        for position in range(8):
            full = 8 if position % 2 else 0
            for m in range(18):
                low, high = max(0, m-17+full), min(full, m//2)
                self.assertEqual(tables[position][m],
                                 tuple(arb(v) for v in (high, 20-low, m, 20-m)))

    def test_counts_and_scope_with_synthetic_caps(self):
        sources = {(0, 0): 31, (8, 8): 7}
        pairs = {i: i+1 for i in range(1, 18)}
        schedule = [(1, 17, '999/1000', '-2')]
        with patch.object(gate, 'matching_joint_counts', return_value=sources), \
             patch.object(gate, 'geometry_stratum_counts', return_value=pairs), \
             patch.object(gate, 'pair_regions_screen', return_value=None), \
             patch.object(gate, 'ordered_tables', return_value=None), \
             patch.object(gate, 'MultiplicityCapsFloat', return_value=lambda mu: {}), \
             patch.object(gate, 'joint_log_count', return_value=0.):
            result = gate.run(schedule=schedule)
            prefix = -log2(gate.Q-1)-gate.LIMIT*log2(.999)+2
            self.assertAlmostEqual(result['failure_log2_screen'],
                prefix+log2(sum(sources.values()))+log2(sum(pairs.values())), places=9)
            self.assertTrue(result['complete_geometry_partition'])
            self.assertFalse(result['whole_code_certificate'])
            partial = gate.run(schedule=schedule, source_pairs=[(8, 8)])
            self.assertFalse(partial['complete_geometry_partition'])
            self.assertAlmostEqual(partial['failure_log2_screen'],
                prefix+log2(7)+log2(sum(pairs.values())), places=9)


if __name__ == '__main__':
    unittest.main()
