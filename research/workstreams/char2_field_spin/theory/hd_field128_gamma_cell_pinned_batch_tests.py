"""Tiny mocked-bank pin-lifecycle tests; no coefficient or distance gates."""

from fractions import Fraction as F
from pathlib import Path
import unittest
from unittest.mock import Mock, patch

from flint import arb, ctx

import hd_field128_gamma_cell_pinned_batch as batch


def evaluator():
    regional = batch.regional
    bank = object.__new__(regional.provider.AliasBank)
    bank.precision, bank.constant, bank.slope = 256, arb(-1000), arb(2)
    bank.construction, bank.state, bank.q = 'synthetic-test-only', 30, 1 << 128
    bank.distance, bank.forbidden = F(2, 5), 838860
    bank.z, bank.alpha, bank.grid, bank.bits = F(1, 2), F(1, 8), 128, 80
    bank.markers = ((F(1, 2), F(1, 3)),)
    bank.pins, bank.check_pins = {'mock': 'not-a-gate'}, Mock()
    with ctx.workprec(256):
        rows = {f: (regional.base.rational(F(1, 4)),)*4 for f in range(3)}
    with patch.object(regional.restored, 'RegionalCache') as cache:
        cache.return_value.regional.return_value = rows
        result = batch.PinnedBatchEvaluator(bank, 4)
    bank.check_pins.reset_mock()
    return result, bank


class PinnedBatchTests(unittest.TestCase):
    def test_frozen_backbone_with_exactly_two_bank_checks_and_partial_checkpoints(self):
        value, bank = evaluator()
        saved = []
        with patch.object(batch.frozen_cover, 'cover', wraps=batch.frozen_cover.cover) as driver:
            result = value.run_cover(max_nodes=10, max_seconds=30., cell_seconds=5., checkpoint=saved.append)
        driver.assert_called_once()
        self.assertIs(driver.call_args.args[0], value)
        self.assertEqual(driver.call_args.args[1], 4)
        self.assertEqual(bank.check_pins.call_count, 2)
        self.assertFalse(value._batch_active)
        self.assertTrue(result['complete_selected_weight'])
        self.assertTrue(result['batch_source_verified_at_exit'])
        self.assertTrue(result['directed'])
        self.assertTrue(saved)
        for checkpoint in saved:
            self.assertFalse(checkpoint['directed'])
            self.assertFalse(checkpoint['complete_selected_weight'])
            self.assertFalse(checkpoint['batch_source_verified_at_exit'])
            self.assertIsNone(checkpoint['log2_upper_rational'])
            self.assertTrue(checkpoint['cell_receipts_conditional_on_final_source_check'])
            for node in checkpoint['changed_nodes']:
                if 'receipt' in node:
                    self.assertEqual(node['receipt']['source_verification'], 'pending-batch-exit')
        self.assertIn('hd_field128_gamma_cell_pinned_batch.py', result['batch_source_sha256'])
        self.assertIn('hd_field128_gamma_cell_cover.py', result['batch_source_sha256'])
        self.assertIn('hd_field128_gamma_cell_regional_table.py', result['batch_source_sha256'])

    def test_standalone_receipts_and_checks_still_verify_sources(self):
        value, bank = evaluator()
        cell = batch.frozen_cover.selected_weight_roots(4)[0]
        result = value.receipt(cell)
        self.assertEqual(result['source_verification'], 'verified-standalone')
        self.assertEqual(bank.check_pins.call_count, 1)
        value.check_pins()
        self.assertEqual(bank.check_pins.call_count, 2)

    def test_source_mutation_during_batch_rejects_final_return(self):
        value, bank = evaluator()
        changed, saved = [False], []
        read = Path.read_bytes
        def read_with_changed_wrapper(path):
            data = read(path)
            return data+b'\n# synthetic source mutation\n' if changed[0] and path.name == Path(batch.__file__).name else data
        def checkpoint(row):
            saved.append(row)
            changed[0] = True
        with patch.object(Path, 'read_bytes', read_with_changed_wrapper):
            with self.assertRaisesRegex(RuntimeError, 'pinned'):
                value.run_cover(max_nodes=10, max_seconds=30., checkpoint=checkpoint)
        self.assertFalse(value._batch_active)
        self.assertTrue(saved)
        self.assertTrue(all(row['complete_selected_weight'] is False for row in saved))
        value.check_pins()

    def test_exceptions_and_reentrancy_reset_flag_and_check_on_failure(self):
        value, bank = evaluator()
        with patch.object(batch.frozen_cover, 'cover', side_effect=ArithmeticError('mock driver failure')):
            with self.assertRaisesRegex(ArithmeticError, 'mock driver failure'):
                value.run_cover()
        self.assertFalse(value._batch_active)
        self.assertEqual(bank.check_pins.call_count, 2)
        bank.check_pins.reset_mock()
        def reenter(row):
            value.run_cover()
        with self.assertRaisesRegex(RuntimeError, 'already active'):
            value.run_cover(checkpoint=reenter)
        self.assertFalse(value._batch_active)
        self.assertEqual(bank.check_pins.call_count, 2)
        value.check_pins()

    def test_scope_and_checkpoint_completion_guards(self):
        value, bank = evaluator()
        for kwargs in (dict(weight=3), dict(weight=True), dict(pairs=4), dict(checkpoint=1)):
            with self.assertRaises(ValueError):
                value.run_cover(**kwargs)
        def bad_driver(evaluator, weight, *, checkpoint, **kwargs):
            checkpoint(dict(directed=True, complete_selected_weight=True,
                            whole_code_certificate=False, log2_upper_rational='-100'))
        with patch.object(batch.frozen_cover, 'cover', side_effect=bad_driver):
            with self.assertRaisesRegex(RuntimeError, 'checkpoint must remain'):
                value.run_cover()
        self.assertFalse(value._batch_active)
        self.assertEqual(bank.check_pins.call_count, 2)


if __name__ == '__main__':
    unittest.main()
