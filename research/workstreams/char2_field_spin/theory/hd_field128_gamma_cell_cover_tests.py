"""Tiny geometry and fake-evaluator control tests; no real alias banks/gates."""

from contextlib import redirect_stdout
from copy import deepcopy
from fractions import Fraction as F
from io import StringIO
from itertools import product
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import hd_field128_gamma_cell_cover as c


def points(cell):
    return {point for point in product(*(range(lo, hi+1) for lo, hi in cell.ranges)) if cell.contains(point)}


class FakeEvaluator:
    def __init__(self, policy=None):
        self.policy, self.calls, self.pin_checks = policy, [], 0

    def check_pins(self):
        self.pin_checks += 1

    def receipt(self, cell, *, marker_pairs, max_seconds, failure_bits):
        self.calls.append(cell)
        value = self.policy(cell) if self.policy else -80
        row = dict(cell=cell.descriptor(), status='complete', complete_cell=True, directed=True,
                   zero=False, log2_upper_rational=str(value), failure_bits=failure_bits,
                   meets_cell_budget=value <= -failure_bits)
        if value is None:
            raise AssertionError('Use an explicit pending-policy subclass')
        return row


class CellCoverTests(unittest.TestCase):
    def run_cover(self, evaluator, weight=3, pairs=4, **kwargs):
        # Pin machinery is separately tested; avoid repeated closure scans in tiny mocks.
        with patch.object(c, 'source_pins', return_value={}), patch.object(c, 'check_source_pins'):
            return c.cover(evaluator, weight, pairs=pairs, **kwargs)

    def test_canonical_roots_cover_every_tiny_support_exactly_once(self):
        for pairs in range(2, 6):
            for weight in range(1, 2*pairs):
                roots = c.selected_weight_roots(weight, pairs)
                profiles = [{cell.profile(point) for point in points(cell)} for cell in roots]
                for first, group in enumerate(profiles):
                    self.assertFalse(any(group & other for other in profiles[first+1:]))
                for mask in range(1, (1 << (2*pairs))-1):
                    active = [(mask >> i) & 1 for i in range(2*pairs)]
                    if sum(active) != weight:
                        continue
                    a = sum(active[2*i]*active[2*i+1] for i in range(pairs))
                    b = sum(active[(2*i+1) % (2*pairs)]*active[(2*i+2) % (2*pairs)] for i in range(pairs))
                    self.assertEqual(sum((weight, a, b) in group for group in profiles), 1)

    def test_complete_union_and_zero_bound(self):
        fake = FakeEvaluator()
        result = self.run_cover(fake, weight=4)
        self.assertTrue(result['complete_selected_weight'])
        self.assertFalse(result['pending_frontier'])
        self.assertEqual(len(result['accepted_leaves']), 3)
        self.assertTrue(F(-79) < F(result['log2_upper_rational']) < F(-78))
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['directed'])  # Injected interfaces cannot produce production certificates.
        class Zero(FakeEvaluator):
            def receipt(self, cell, **kwargs):
                row = super().receipt(cell, **kwargs)
                row.update(zero=True, log2_upper_rational=None)
                return row
        zero = self.run_cover(Zero())
        self.assertTrue(zero['zero'])
        self.assertIsNone(zero['log2_upper_rational'])
        self.assertTrue(zero['meets_requested_selected_weight_target'])

    def test_splits_complete_partition_and_immutable_delta_replay(self):
        fake = FakeEvaluator(lambda cell: -90 if len(points(cell)) <= 1 else 0)
        checkpoints = []
        result = self.run_cover(fake, checkpoint=lambda row: checkpoints.append(deepcopy(row)), max_nodes=32)
        self.assertTrue(result['complete_selected_weight'])
        self.assertEqual(result['evaluated_nodes'], 7)
        self.assertEqual(len(result['accepted_leaves']), 4)
        replay = {}
        for sequence, row in enumerate(checkpoints):
            self.assertEqual(row['sequence'], sequence)
            for node in row['changed_nodes']:
                replay[node['id']] = node
            accepted, pending = c.audit_tree(list(replay.values()), 3, 4, 64)
            self.assertEqual(set(pending), set(row['pending_frontier']))
            self.assertEqual(set(accepted), set(row['accepted_leaves']))
            self.assertFalse(row['complete_selected_weight'])
            self.assertIsNone(row['log2_upper_rational'])
        self.assertEqual(replay, {node['id']: node for node in result['nodes']})
        self.assertEqual(checkpoints[0]['changed_nodes'][0]['status'], 'pending')
        # Damaging a child or accepting a missing receipt is rejected independently.
        damaged = deepcopy(result['nodes'])
        damaged[0]['children'][1] = damaged[0]['children'][0]
        with self.assertRaises(ValueError):
            c.audit_tree(damaged, 3, 4, 64)

    def test_node_limit_pending_receipt_and_unsplittable_failure(self):
        limited = self.run_cover(FakeEvaluator(lambda cell: 0), max_nodes=1)
        self.assertEqual(limited['stop_reason'], 'max_nodes')
        self.assertEqual(limited['evaluated_nodes'], 1)
        self.assertEqual(len(limited['pending_frontier']), 2)
        self.assertFalse(limited['complete_selected_weight'])
        self.assertIsNone(limited['log2_upper_rational'])
        class Pending(FakeEvaluator):
            def receipt(self, cell, **kwargs):
                return dict(cell=cell.descriptor(), complete_cell=False,
                            log2_upper_rational='-10000', meets_cell_budget=True)
        pending = self.run_cover(Pending())
        self.assertFalse(pending['accepted_leaves'])
        self.assertEqual(pending['nodes'][0]['pending_reason'], 'incomplete-cell-receipt')
        self.assertIsNone(pending['log2_upper_rational'])
        terminal = self.run_cover(FakeEvaluator(lambda cell: 0), weight=1)
        self.assertEqual(terminal['nodes'][0]['pending_reason'], 'unsplittable-over-budget')
        self.assertFalse(terminal['complete_selected_weight'])

    def test_time_budget_and_pin_failure(self):
        clock = [0.]
        class Slow(FakeEvaluator):
            def receipt(self, cell, **kwargs):
                clock[0] += 2.
                return super().receipt(cell, **kwargs)
        with patch.object(c, 'monotonic', side_effect=lambda: clock[0]):
            result = self.run_cover(Slow(), weight=4, max_seconds=1.)
        self.assertEqual(result['evaluated_nodes'], 1)
        self.assertEqual(result['stop_reason'], 'max_seconds')
        self.assertEqual(len(result['pending_frontier']), 2)
        self.assertIsNone(result['log2_upper_rational'])
        class Changed(FakeEvaluator):
            def check_pins(self):
                super().check_pins()
                if self.pin_checks == 2:
                    raise RuntimeError('Changed test source')
        with self.assertRaises(RuntimeError):
            self.run_cover(Changed())
        with self.assertRaises(RuntimeError):
            self.run_cover(FakeEvaluator(), input_check=lambda: (_ for _ in ()).throw(RuntimeError('Changed input')))

    def test_witness_literals_guards_and_exclusive_checkpoints(self):
        artifact = {'witness': {'z': '1/2', 'alpha': '1/3', 'markers': [['1/4', '1/5']]},
                    'log2_upper_rational': '-1000000000000000'}
        self.assertEqual(c.witness_literals(artifact), (F(1, 2), F(1, 3), ((F(1, 4), F(1, 5)),)))
        bad = deepcopy(artifact)
        bad['witness']['z'] = .5
        with self.assertRaises(ValueError):
            c.witness_literals(bad)
        for kwargs in (dict(max_nodes=0), dict(max_seconds=0), dict(cell_seconds=float('nan')),
                       dict(failure_bits=True), dict(precision=64)):
            with self.assertRaises(ValueError):
                self.run_cover(FakeEvaluator(), **kwargs)
        with self.assertRaises(ValueError):
            c.selected_weight_roots(8, 4)
        with TemporaryDirectory() as temporary:
            path = Path(temporary)/'checkpoint.json'
            c.write_new(path, {'a': 1})
            with self.assertRaises(FileExistsError):
                c.write_new(path, {'a': 2})
            self.assertEqual(json.loads(path.read_text()), {'a': 1})

    def test_mismatched_or_false_complete_receipt_is_rejected(self):
        class WrongCell(FakeEvaluator):
            def receipt(self, cell, **kwargs):
                row = super().receipt(cell, **kwargs)
                row['cell']['pairs'] += 1
                return row
        class WrongBudget(FakeEvaluator):
            def receipt(self, cell, **kwargs):
                row = super().receipt(cell, **kwargs)
                row['log2_upper_rational'] = '0'
                return row
        for fake in (WrongCell(), WrongBudget()):
            with self.assertRaises(ValueError):
                self.run_cover(fake)

    def test_public_cli_builds_fresh_bank_from_literals_and_partial_exit(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root/'witness.json'
            source.write_text(json.dumps({'witness': {'z': '1/2', 'alpha': '1/3',
                                                      'markers': [['1/4', '1/5']]}}))
            output = root/'output'
            bank = SimpleNamespace(pins={}, metadata=lambda: {'mock': True})
            evaluator = FakeEvaluator(lambda cell: 0)
            argv = ['cell-cover', '--witness-source', str(source), '--r', '3', '--state', '30',
                    '--distance', '2/5', '--max-nodes', '1', '--output-dir', str(output)]
            with patch('sys.argv', argv), patch.object(c.provider, 'AliasBank', return_value=bank) as factory, \
                    patch.object(c.numeric, 'Evaluator', return_value=evaluator), \
                    patch.object(c, 'source_pins', return_value={}), patch.object(c, 'check_source_pins'), \
                    redirect_stdout(StringIO()):
                # Exercise public setup/build/finalization without any numerical bank.
                with patch.object(c, 'cover', side_effect=lambda *args, **kwargs: {'status': 'partial',
                        'evaluated_nodes': 1, 'complete_selected_weight': False, 'log2_upper_rational': None,
                        'stop_reason': 'max_nodes', 'whole_code_certificate': False,
                        'pending_frontier': ['n000001'], 'accepted_leaves': []}):
                    with self.assertRaises(SystemExit) as exit_status:
                        c.main()
                self.assertEqual(exit_status.exception.code, 2)
            self.assertEqual(factory.call_args.args, (F(1, 2), F(1, 3), ((F(1, 4), F(1, 5)),)))
            summary = json.loads((output/'summary.json').read_text())
            self.assertFalse(summary['complete_selected_weight'])
            self.assertIsNone(summary['log2_upper_rational'])
            self.assertIn('sha256', summary['witness_source'])


if __name__ == '__main__':
    unittest.main()
