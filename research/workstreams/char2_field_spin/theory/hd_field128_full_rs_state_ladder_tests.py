"""Mocked resource/provenance tests; no numerical search or distance gate."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import hd_field128_full_rs_state_ladder as ladder


def fake_weight(r, z, alpha, *, state, q, distance, max_seconds):
    return dict(construction=ladder.state_provider.construction(state), state_dimension=state,
                active_groups=r, z=z, alpha=alpha, field_order=q, target_distance=str(distance),
                log2_upper_screen=float(-state-r), worst_profile=(-state-r, 0, 0),
                resource_limit_reached=False, directed=False, whole_code_certificate=False)


def fake_proposal(provider, r, z, alpha, distance, evaluations, seconds):
    row = provider.weight_screen(r, z, alpha, distance=distance, max_seconds=seconds)
    row.update(proposal_evaluations=1, max_evaluations=evaluations)
    return row


class LadderTests(unittest.TestCase):
    def test_plan_seeds_and_scope_guards(self):
        plan = ladder.checked_plan(ladder.DEFAULT_CANDIDATES, ladder.DEFAULT_WEIGHTS)
        self.assertEqual(len(plan), 14)
        self.assertEqual(plan[0], (24, F(2, 5), 1))
        self.assertEqual(plan[-1], (32, F(9, 20), 64))
        z, alpha = ladder.initial_witness(1)
        self.assertEqual(z, float(F(1199, 1200)**16))
        self.assertEqual(alpha, 1/16)
        self.assertEqual(ladder.initial_witness(1, 256), (float(F(1199, 1200)), 1.))
        for candidates, weights in ((((24, F(1, 2)),), (1,)), (((20, F(2, 5)),), (1,)),
                                    (ladder.DEFAULT_CANDIDATES, (1, 1)), (ladder.DEFAULT_CANDIDATES, (65,))):
            with self.assertRaises(ValueError):
                ladder.checked_plan(candidates, weights)

    def test_one_snapshot_per_weight_and_state_adapter(self):
        snapshots = []
        with patch.object(ladder.search, 'proposal', side_effect=fake_proposal), \
             patch.object(ladder.state_provider, 'weight_screen', side_effect=fake_weight):
            result = ladder.run_ladder(weights=(1, 4), evaluations=3, max_seconds=10, checkpoint=snapshots.append)
        self.assertEqual([row['completed_selected_weights'] for row in snapshots], [0, 1, 2, 3, 4])
        self.assertEqual(result['pending_selected_weights'], [])
        self.assertTrue(result['every_requested_weight_evaluated'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['directed'])
        self.assertEqual([row['requested_state'] for row in result['results']], [24, 24, 32, 32])
        for row in result['results']:
            self.assertEqual(row['state_dimension'], row['requested_state'])
            self.assertEqual(row['target_distance'], row['requested_distance'])
            self.assertEqual(row['max_evaluations'], 3)
        for name in ('hd_field128_candidate_proposal.py', 'hd_field128_full_rs_state_screen.py',
                     'hd_field128_full_rs_state_ladder.py'):
            self.assertIn(name, result['source_sha256'])

    def test_failed_weight_does_not_erase_earlier_or_later_results(self):
        def sometimes(provider, r, *args):
            if r == 2:
                raise TimeoutError('Mocked weight had no complete evaluation')
            return fake_proposal(provider, r, *args)
        snapshots = []
        with patch.object(ladder.search, 'proposal', side_effect=sometimes), \
             patch.object(ladder.state_provider, 'weight_screen', side_effect=fake_weight):
            result = ladder.run_ladder(candidates=((32, F(12, 25)),), weights=(1, 2, 4),
                                      max_seconds=10, checkpoint=snapshots.append)
        self.assertEqual([row['active_groups'] for row in result['results']], [1, 4])
        self.assertEqual([row['active_groups'] for row in result['pending_selected_weights']], [2])
        self.assertEqual(result['attempted_job_indices'], [0, 1, 2])
        self.assertEqual(result['failures'][0]['job_index'], 1)
        self.assertEqual(len(snapshots), 4)
        self.assertFalse(result['every_requested_weight_evaluated'])

    def test_global_budget_stops_before_any_search(self):
        snapshots = []
        with patch.object(ladder.search, 'proposal') as mocked:
            result = ladder.run_ladder(max_seconds=1e-12, checkpoint=snapshots.append)
        mocked.assert_not_called()
        self.assertEqual(result['completed_selected_weights'], 0)
        self.assertEqual(len(result['pending_selected_weights']), 14)
        self.assertEqual(result['failures'][0]['error_type'], 'TimeoutError')
        self.assertEqual(len(snapshots), 1)
        with self.assertRaises(ValueError):
            ladder.run_ladder(initial_z=.5)
        with self.assertRaises(ValueError):
            ladder.run_ladder(max_weight_seconds=3601)

    def test_checkpoint_bytes_are_exclusive_and_hash_pinned(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)/'checkpoint.json'
            receipt = ladder.write_new(path, dict(directed=False, scope='mock only'))
            self.assertEqual(receipt['sha256'], sha256(path.read_bytes()).hexdigest())
            self.assertEqual(json.loads(path.read_bytes())['scope'], 'mock only')
            with self.assertRaises(FileExistsError):
                ladder.write_new(path, {})


if __name__ == '__main__':
    unittest.main()
