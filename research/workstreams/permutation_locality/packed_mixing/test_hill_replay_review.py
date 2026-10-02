"""Independent small replay-boundary tests; no numerical proof is run."""
import copy
from fractions import Fraction as Q
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import arb, ctx
import hill_cover
import whole_hill_replay as whole
from test_hill_cover import Model, record
from test_whole_hill_replay import fixture, fresh


class IndependentReplayBoundaryTests(unittest.TestCase):
    def setUp(self):
        old_precision = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', old_precision)
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        dense_path, sparse_path = (Path(directory.name)/name for name in ('dense.json', 'sparse.json'))
        dense, sparse = fixture(updates=4)
        dense_path.write_text(json.dumps(dense))
        sparse_path.write_text(json.dumps(sparse))
        self.prepared = whole.prepare(dense_path, sparse_path)

    def test_strict_whole_threshold_depends_on_every_fresh_endpoint(self):
        sparse, dense = copy.deepcopy(fresh(self.prepared))
        # Dense=2^-42; sparse=3*2^-42. Each row separately passes the child
        # target41, but the full sum equals2^-40 and must be rejected.
        for row in sparse['results']:
            row['upper'] = [3, -47]
        sparse['aggregate_upper'] = [3, -42]
        for row in dense['fresh_replay']['checked']:
            row['upper'] = [1, -43]
        dense['fresh_replay']['aggregate_upper'] = [1, -42]
        with self.assertRaises(ValueError):
            whole._assemble_fresh(self.prepared, sparse, dense, 384, 40)
        # A tiny change in the final occupancy must be included exactly,
        # even though all displayed low-precision margins would round alike.
        sparse['results'][-1]['upper'] = [3*(1 << 53)-1, -100]
        result = whole._assemble_fresh(self.prepared, sparse, dense, 384, 40)
        self.assertEqual(whole.dyadic(result['upper']), Q(2)**-40-Q(2)**-100)
        self.assertEqual(whole.dyadic(result['sparse_upper']), Q(3)*Q(2)**-42-Q(2)**-100)
        self.assertEqual(result['construction']['inner']['updates'], 4)

    def test_final_upward_rounding_cannot_preserve_an_unrepresentable_strict_margin(self):
        sparse, dense = copy.deepcopy(fresh(self.prepared))
        sparse['precision'] = dense['precision'] = dense['fresh_replay']['precision'] = 128
        for row in sparse['results']:
            row['upper'] = [3, -47]
        sparse['results'][-1]['upper'] = [3*(1 << 453)-1, -500]
        sparse['aggregate_upper'] = [3, -42]
        for row in dense['fresh_replay']['checked']:
            row['upper'] = [1, -43]
        dense['fresh_replay']['aggregate_upper'] = [1, -42]
        exact = Q(2)**-40-Q(2)**-500
        self.assertLess(exact, Q(2)**-40)
        self.assertEqual(whole.dyadic(whole.compact_endpoint(exact, 160)), Q(2)**-40)
        with self.assertRaises(ValueError):
            whole._assemble_fresh(self.prepared, sparse, dense, 128, 40)

    def test_dense_replay_rejects_actual_update_change_before_completing(self):
        class MutatedInner(Model):
            def outward(self, cell, witness):
                self.data['updates'] = 2
                return arb(2)**-100

        snapshots = []
        with patch.object(hill_cover, 'fresh_model', return_value=MutatedInner(4)), \
                self.assertRaises(ArithmeticError):
            hill_cover.replay(record(4, complete=True), 384, 41, snapshots.append)
        self.assertTrue(snapshots)
        self.assertTrue(all(not snapshot['complete_dense'] for snapshot in snapshots))
        self.assertTrue(all(not snapshot['passed'] for snapshot in snapshots))


if __name__ == '__main__':
    unittest.main()
