"""Pure guards and mocked forwarding: no banks, coefficients or output files."""
from contextlib import contextmanager, ExitStack
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import hd_field128_gamma_portable_band_runner_v27 as runner


def complete(bound='-60'):
    return dict(status='complete', complete_requested_band=True, directed=True,
        whole_code_certificate=False, requested_failure_bits=54, pending_frontier=[],
        meets_requested_band_target=True, zero=False, log2_upper_rational=bound,
        requested_weight_range=[960, 999], evaluated_nodes=4, stop_reason=None)


class Tests(unittest.TestCase):
    def test_scope_and_threshold_guards_before_any_work(self):
        for threshold in (54, 58, 64, 72, 128):
            runner.validate(960, 999, 4096, threshold)
        for threshold in (True, 53, 129, 64., '64'):
            with patch.object(runner, 'hard_deadline', side_effect=AssertionError('unexpected work')):
                with self.assertRaisesRegex(ValueError, 'cell_failure_bits'):
                    runner.run('seed.json', '0'*64, 960, 999, 'unused', cell_failure_bits=threshold)
        for first, last, nodes in ((0, 1, 4096), (1, 65, 4096), (8192, 8192, 4096), (960, 999, 4097)):
            with self.assertRaises(ValueError):
                runner.validate(first, last, nodes, 64)

    def test_record_changes_only_wrapper_identity_and_cell_acceptance(self):
        pins = {runner.NAME:'1'*64, runner.previous.NAME:'2'*64, runner.controller.NAME:'3'*64}
        witness = dict(path='seed.json', sha256='4'*64)
        old = runner.previous.driver_record(pins, witness, 960, 999, 4096, 960)
        new = runner.driver_record(pins, witness, 960, 999, 4096, 960, 64)
        different = {key for key in old.keys() & new.keys() if old[key] != new[key]}
        self.assertEqual(different, {'schema', 'driver', 'cell_failure_bits', 'wrapper_scope'})
        self.assertEqual(new['driver'], dict(path=runner.NAME, sha256=pins[runner.NAME]))
        self.assertEqual((new['cell_failure_bits'], new['aggregate_failure_bits']), (64, 54))
        self.assertEqual(new['partition_controller_sha256'], pins[runner.controller.NAME])
        self.assertEqual(new['witness_source'], witness)
        self.assertEqual(new['requested_weight_range'], [960, 999])
        self.assertEqual(new['max_output_bytes'], 32 << 20)
        self.assertEqual(new['max_table_bytes'], 64 << 20)
        self.assertEqual(new['max_address_space_bytes'], 4 << 30)
        self.assertEqual(new['disk_reserve_bytes'], 20 << 30)
        self.assertEqual(new['hard_wall_seconds'], 600.)
        self.assertFalse(new['fresh_bank_rebuilt'])
        self.assertFalse(new['whole_code_certificate'])

    def test_admission_requires_actual_54_bit_union_and_complete_tree(self):
        self.assertTrue(runner.admitted(complete('-54')))
        self.assertTrue(runner.admitted(dict(complete(), zero=True, log2_upper_rational=None)))
        for fields in (dict(status='partial'), dict(complete_requested_band=False), dict(directed=False),
                dict(whole_code_certificate=True), dict(requested_failure_bits=53),
                dict(pending_frontier=[1]), dict(meets_requested_band_target=False),
                dict(log2_upper_rational='-53999/1000'), dict(log2_upper_rational='nan'),
                dict(zero=True), dict(zero=None), dict(log2_upper_rational=None)):
            self.assertFalse(runner.admitted(dict(complete(), **fields)), fields)

    def test_source_closure_and_frozen_wrapper_mutation(self):
        pins = runner.source_pins()
        self.assertEqual(pins[runner.previous.NAME], runner.FROZEN_PREVIOUS)
        self.assertEqual(pins[runner.controller.NAME], runner.previous.FROZEN_CONTROLLER)
        self.assertEqual(pins[runner.NAME], sha256(Path(runner.__file__).read_bytes()).hexdigest())
        runner.controller.check_source_pins(pins)
        read, target = Path.read_bytes, Path(runner.previous.__file__)
        def changed(path):
            raw = read(path)
            return raw+b'\n# changed\n' if path == target else raw
        with patch.object(Path, 'read_bytes', changed):
            with self.assertRaisesRegex(RuntimeError, 'Frozen v26'):
                runner.source_pins()

    def test_linux_hard_600_second_exit_and_reset(self):
        fake_signal = SimpleNamespace(SIGALRM=14, ITIMER_REAL=0, getsignal=lambda unused:'old',
            signal=Mock(), setitimer=Mock())
        with patch.object(runner.sys, 'platform', 'linux'), patch.object(runner, 'signal', fake_signal), \
                patch.object(runner.os, '_exit') as exit_process:
            with runner.hard_deadline():
                fake_signal.signal.call_args_list[0].args[1](None, None)
                exit_process.assert_called_once_with(124)
                self.assertEqual(fake_signal.setitimer.call_args_list[0].args, (0, 600.))
            self.assertEqual(fake_signal.setitimer.call_args_list[-1].args, (0, 0))
            self.assertEqual(fake_signal.signal.call_args_list[-1].args, (14, 'old'))

    def mock_run(self, result=None, fail_cover=False):
        writes, active = [], []
        pins = {runner.NAME:'1'*64, runner.previous.NAME:'2'*64, runner.controller.NAME:'3'*64}
        witness = dict(path='seed.json', sha256='4'*64)
        seed = dict(active_groups=960)
        bank = SimpleNamespace(metadata=lambda:dict(witness={'mock':True}), check_pins=Mock())
        evaluator = object()
        @contextmanager
        def deadline():
            active.append(True)
            try:
                yield
            finally:
                active.pop()
        def build(*args, **kwargs):
            self.assertTrue(active)
            self.assertEqual(args, (F(1, 4), F(1, 8), [(F(2), F(3))]))
            self.assertEqual((kwargs['state'], kwargs['grid'], kwargs['precision'], kwargs['bits']), (30, 128, 256, 80))
            self.assertEqual(kwargs['max_table_bytes'], 64 << 20)
            return bank
        def cover(actual, **kwargs):
            self.assertIs(actual, evaluator)
            self.assertEqual((kwargs['cell_failure_bits'], kwargs['aggregate_bits'], kwargs['max_nodes']), (64, 54, 4096))
            kwargs['context_checkpoint'](dict(partial_context=True))
            kwargs['checkpoint'](dict(sequence=0, complete_requested_band=False))
            if fail_cover:
                raise TimeoutError('mock deadline')
            return deepcopy(result or complete())
        output = SimpleNamespace(write_new=lambda path, value:writes.append((path.name, deepcopy(value))))
        with ExitStack() as stack:
            mocks = {}
            for owner, name, replacement in ((runner, 'hard_deadline', deadline),
                    (runner, 'source_pins', lambda:dict(pins)),
                    (runner.controller, 'load_seed', lambda *args:(seed, (F(1,4), F(1,8), [(F(2),F(3))]), witness)),
                    (runner.controller, 'OutputBudget', lambda unused:output),
                    (runner.controller.cert, 'AliasBank', build),
                    (runner.controller.bands, 'Evaluator', Mock(return_value=evaluator)),
                    (runner.controller, 'cover', cover)):
                stack.enter_context(patch.object(owner, name, replacement))
            for owner, name in ((runner.previous, 'check_pins'), (runner.controller, 'reserve_guard'),
                    (runner.controller, 'install_memory_limit'), (Path, 'mkdir')):
                mocks[name] = stack.enter_context(patch.object(owner, name))
            if fail_cover:
                with self.assertRaises(TimeoutError):
                    runner.run('seed.json', '4'*64, 960, 999, 'never-created', cell_failure_bits=64)
            else:
                saved = runner.run('seed.json', '4'*64, 960, 999, 'never-created', cell_failure_bits=64)
                self.assertEqual(saved['runner_admitted'], runner.admitted(saved))
            mocks['reserve_guard'].assert_called_once_with(runner.HERE, 32 << 20)
            mocks['install_memory_limit'].assert_called_once()
            self.assertGreaterEqual(mocks['check_pins'].call_count, 3)
            self.assertGreaterEqual(bank.check_pins.call_count, 2)
        self.assertEqual(active, [])
        return writes

    def test_mock_forwarding_limits_and_truthful_metadata(self):
        writes = dict(self.mock_run())
        self.assertEqual(set(writes), {'driver.json', 'bank.json', 'context.json', 'node_000000.json', 'summary.json'})
        self.assertEqual(writes['driver.json']['cell_failure_bits'], 64)
        self.assertTrue(writes['bank.json']['fresh_bank_rebuilt'])
        self.assertTrue(writes['summary.json']['runner_admitted'])
        self.assertFalse(writes['summary.json']['whole_code_certificate'])

    def test_partial_failure_and_nonpassing_complete_are_preserved_honestly(self):
        for result in (dict(complete(), status='partial', complete_requested_band=False,
                pending_frontier=[1], meets_requested_band_target=False, log2_upper_rational=None),
                dict(complete('-53'), meets_requested_band_target=False)):
            writes = dict(self.mock_run(result))
            self.assertFalse(writes['summary.json']['runner_admitted'])
        writes = dict(self.mock_run(fail_cover=True))
        self.assertIn('node_000000.json', writes)
        self.assertNotIn('summary.json', writes)

    def test_cli_exit_requires_aggregate_not_just_completeness(self):
        argv = ['runner', '--seed', 'seed.json', '--expected-sha256', '4'*64, '--first', '960',
                '--last', '999', '--output-dir', 'never-created', '--cell-failure-bits', '64']
        for value, code in ((complete(), 0), (dict(complete('-53'), meets_requested_band_target=False), 2)):
            with patch.object(runner.sys, 'argv', argv), patch.object(runner, 'run', return_value=value) as run, \
                    patch('builtins.print'):
                self.assertEqual(runner.main(), code)
                self.assertEqual(run.call_args.kwargs['cell_failure_bits'], 64)
        with patch.object(runner.sys, 'argv', argv[:-2]), patch.object(runner, 'run', return_value=complete()) as run, \
                patch('builtins.print'):
            self.assertEqual(runner.main(), 0)
            self.assertEqual(run.call_args.kwargs['cell_failure_bits'], 72)


if __name__ == '__main__':
    unittest.main()
