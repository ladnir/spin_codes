"""Small mocked-process replay tests; no expensive proof computation."""
import copy
from contextlib import redirect_stdout
from fractions import Fraction as Q
import io
import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from flint import arb, ctx
import hill_cover as hill
import whole_hill_replay as whole
from test_hill_cover import record
from test_whole_hill_replay import fixture, fresh


class Model:
    def __init__(self, scope):
        self.root = tuple(map(Q, scope['root']))
        self.data = dict(updates=scope['updates'], windows=32)
        self.threshold = scope['threshold']
        self.q_min = scope['minimum_groups']
        self.tilt = Q(scope['base_tilt'])
        self.variance_bins = scope['variance_bins']
        self.variance_shuffle = self.regional_count = True
        self.calls = []

    def outward(self, cell, witness):
        self.calls.append((cell, copy.deepcopy(witness)))
        return arb(2)**-60


def tasks(source):
    return [(path, row['cell'], row['witness']) for path, row in sorted(source['cover']['leaves'].items())]


class HillParallelReplayTests(unittest.TestCase):
    def setUp(self):
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)
        quiet = redirect_stdout(io.StringIO())
        quiet.__enter__()
        self.addCleanup(quiet.__exit__, None, None, None)

    def rows(self, source, jobs=None):
        return [hill._evaluate_replay(Model(source['scope']), source['scope'], 384, task)
                for task in tasks(source) if jobs is None or task[0] in jobs]

    def test_parent_spawn_fresh_rows_exact_sum_and_unchanged_inputs(self):
        source = record(4, complete=True)
        original = copy.deepcopy(source)
        executor = Mock()
        def mapped(function, assigned, chunksize):
            self.assertIs(function, hill._evaluate_replay_worker)
            self.assertEqual(chunksize, 1)
            self.assertEqual(assigned, tasks(source))
            return iter(self.rows(source))
        executor.map.side_effect = mapped
        checkpoints = []
        path, environment = sys.path[:], dict(os.environ)
        with patch.object(hill, 'ProcessPoolExecutor', return_value=executor) as pool, \
                patch.object(hill, 'get_context', return_value='SPAWN') as context, \
                patch.object(hill, 'fresh_model') as parent_build:
            result = hill.replay(source, 384, 41, checkpoints.append, workers=3)
        parent_build.assert_not_called()
        context.assert_called_once_with('spawn')
        self.assertEqual(pool.call_args.kwargs['max_workers'], 3)
        self.assertIs(pool.call_args.kwargs['initializer'], hill._initialize_replay_worker)
        self.assertEqual(pool.call_args.kwargs['initargs'], (source['scope'], 384, None))
        self.assertTrue(result['complete_dense'] and result['passed'])
        self.assertEqual(sum(hill.cell_search.dense.dyadic(r['upper']) for r in result['checked']), Q(1, 1 << 59))
        self.assertTrue(all(not r['complete_dense'] for r in checkpoints[:-1]))
        self.assertEqual(source, original)
        self.assertEqual(sys.path, path)
        self.assertEqual(dict(os.environ), environment)
        executor.shutdown.assert_called_once_with(wait=True, cancel_futures=True)
        executor.terminate_workers.assert_not_called()

    def test_worker_result_coverage_and_all_scope_fields_fail_closed(self):
        source = record(4, complete=True)
        good = self.rows(source)
        corruptions = [dict(path='1'), dict(cell=['0', '1']), dict(precision=256),
            dict(updates=2), dict(updates=True), dict(scope_sha256='wrong'),
            dict(witness_sha256='wrong'), dict(upper=[0, 1]), dict(upper=[True, -60])]
        cases = [good[:1], good+[good[0]], [good[0], good[0]], good[::-1]]
        for corruption in corruptions:
            changed = copy.deepcopy(good)
            changed[0].update(corruption)
            cases.append(changed)
        for rows in cases:
            with self.subTest(rows=rows), self.assertRaises(ArithmeticError):
                list(hill._checked_replay_rows(source['scope'], 384, tasks(source), rows))

    def test_worker_model_and_precision_drift_rejected(self):
        source = record(4, complete=True)
        for field, value in [('threshold', 0), ('q_min', 34), ('variance_bins', 64),
                             ('root', (Q(0), Q(1, 2)))]:
            model = Model(source['scope'])
            setattr(model, field, value)
            with self.subTest(field=field), self.assertRaises(ArithmeticError):
                hill._evaluate_replay(model, source['scope'], 384, tasks(source)[0])
        for mutation in ('precision', 'updates', 'witness', 'nonfinite'):
            model = Model(source['scope'])
            def outward(cell, witness):
                if mutation == 'precision': ctx.prec = 128
                elif mutation == 'updates': model.data['updates'] = 2
                elif mutation == 'witness': witness['changed'] = True
                return arb('nan') if mutation == 'nonfinite' else arb(2)**-60
            model.outward = outward
            with self.subTest(mutation=mutation), self.assertRaises(ArithmeticError):
                hill._evaluate_replay(model, source['scope'], 384, tasks(source)[0])

    def test_initializer_authenticates_independently_once_and_requires_plain_model(self):
        source = record(4, complete=True)
        model = Model(source['scope'])
        def build(scope, precision):
            ctx.prec = precision
            return model
        with patch.object(hill, 'fresh_model', side_effect=build) as fresh_model, \
                patch.dict(sys.modules, scalar_cover=SimpleNamespace(Model=Model)):
            hill._initialize_replay_worker(source['scope'], 384, None)
            rows = [hill._evaluate_replay_worker(task) for task in tasks(source)]
        fresh_model.assert_called_once_with(source['scope'], 384)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(model.calls), 2)
        with patch.object(hill, 'fresh_model', side_effect=build), \
                patch.dict(sys.modules, scalar_cover=SimpleNamespace(Model=object)), \
                self.assertRaises(ArithmeticError):
            hill._initialize_replay_worker(source['scope'], 384, None)

    def test_worker_exception_terminates_pool_and_never_marks_complete(self):
        source = record(4, complete=True)
        def failed(*args, **kwargs):
            yield self.rows(source)[0]
            raise RuntimeError('worker failed')
        executor = Mock()
        executor.map.side_effect = failed
        checkpoints = []
        with patch.object(hill, 'ProcessPoolExecutor', return_value=executor), \
                self.assertRaisesRegex(RuntimeError, 'worker failed'):
            hill.replay(source, 384, 41, checkpoints.append, workers=2)
        self.assertTrue(checkpoints)
        self.assertTrue(all(not row['complete_dense'] and not row['passed'] for row in checkpoints))
        executor.terminate_workers.assert_called_once_with()
        executor.shutdown.assert_called_once_with(wait=True, cancel_futures=True)

    def test_invalid_partition_worker_count_and_existing_logs_never_launch(self):
        with patch.object(hill, 'ProcessPoolExecutor') as pool:
            for workers in (0, 5, True, 2.0):
                with self.assertRaises(ValueError):
                    hill.replay(record(4, complete=True), workers=workers)
            with self.assertRaises(ValueError):
                hill.replay(record(4), workers=3)
            malformed = record(4, complete=True)
            malformed['cover']['leaves']['1']['cell'] = ['0', '1']
            with self.assertRaises(ValueError):
                hill.replay(malformed, workers=3)
            with TemporaryDirectory() as directory, self.assertRaises(ValueError):
                hill.replay(record(4, complete=True), workers=3, logdir=directory)
            pool.assert_not_called()

    def test_whole_forwards_only_dense_workers_and_retains_serial_child_order(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            d, s = fixture()
            dp, sp = root/'d.json', root/'s.json'
            dp.write_text(json.dumps(d)); sp.write_text(json.dumps(s))
            prepared = whole.prepare(dp, sp)
            sparse, dense = fresh(prepared)
            dense['fresh_replay']['workers'] = 3
            commands = []
            def child(command, **kwargs):
                commands.append(command)
                output = Path(command[command.index('--output')+1])
                output.write_text(json.dumps(sparse if len(commands) == 1 else dense))
            with patch.object(whole.subprocess, 'run', side_effect=child):
                result = whole.run(dp, sp, root/'out', dense_workers=3)
            self.assertNotIn('--replay-workers', commands[0])
            self.assertEqual(commands[1][-2:], ['--replay-workers', '3'])
            self.assertEqual(result['dense_replay_workers'], 3)
            self.assertTrue(result['verified'])
            with patch.object(whole.subprocess, 'run') as launch:
                for workers in (0, 5, True, 3.0):
                    with self.assertRaises(ValueError):
                        whole.run(dp, sp, root/'bad', dense_workers=workers)
                launch.assert_not_called()


if __name__ == '__main__':
    unittest.main()
