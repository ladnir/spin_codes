import copy
from fractions import Fraction as Q
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from flint import arb, ctx
import replay_pool as pool


def candidate():
    return dict(tilt='3/16', parameters=['1/100', '0', '0'], variance_dual=['0', '0'])


def jobs():
    return [('1', (Q(1, 2), Q(1)), candidate()), ('0', (Q(0), Q(1, 2)), candidate())]


def prepared():
    return pool.prepare([('100', '1/2'), (Q(50), Q(3, 4))], 199229, 33, 384, jobs())


class ToyModel:
    root = (Q(0), Q(1))

    def outward(self, cell, witness):
        self.seen = (cell, witness, ctx.prec)
        return arb(2)**-50


class ReplayPoolTests(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', self.precision)

    def fresh_rows(self, scope, tasks):
        return [pool.evaluate(ToyModel(), scope, task) for task in tasks]

    def test_exact_inputs_copied_and_deterministically_ordered(self):
        original = jobs()
        snapshot = copy.deepcopy(original)
        scope, tasks = pool.prepare([dict(mass='200/2', activity='2/4')], 199229, 33, 384, original)
        self.assertEqual([task[0] for task in tasks], ['0', '1'])
        self.assertEqual(scope['mixture'], [dict(mass='100', activity='1/2')])
        tasks[0][2]['parameters'][0] = '1/50'
        self.assertEqual(original, snapshot)
        scope2, _ = pool.prepare([(100, Q(1, 2))], 199229, 33, 384, original)
        self.assertEqual(scope['scope_sha256'], scope2['scope_sha256'])

    def test_bad_scope_and_mixtures_fail_before_worker_launch(self):
        for threshold, qmin, precision in ((True, 33, 384), (-1, 33, 384),
                (1 << 20, 33, 384), (1, 0, 384), (1, 2049, 384),
                (1, True, 384), (1, 33, 127), (1, 33, 384.0)):
            with self.assertRaises(ValueError):
                pool.prepare([(1, '1/2')], threshold, qmin, precision, jobs())
        for mixture in ([], [(0, '1/2')], [(1, 0)], [(1, 2)], [(1, True)],
                [(1.0, '1/2')], [(1, '1/2'), (2, '2/4')], [(1, 'nan')],
                [dict(mass='1', activity='1/2', extra=True)]):
            with self.assertRaises(ValueError):
                pool.prepare(mixture, 199229, 33, 384, jobs())

    def test_reject_malformed_duplicate_or_overlapping_tasks(self):
        variants = [[], [jobs()[0], jobs()[0]], [('x', (0, 1), candidate())],
            [('', (0, 0), candidate())], [('', (0.0, 1), candidate())],
            [('', (0, 1), dict(candidate(), parameters=['0', '0', '0']))],
            [('', (0, 1), dict(candidate(), parameters=['1', '0', '-1']))],
            [('', (0, 1), dict(candidate(), ignored_float=1.2))],
            [('', (0, 1), candidate()), ('0', (0, Q(1, 2)), candidate())],
            [('0', (0, Q(3, 4)), candidate()), ('1', (Q(1, 2), 1), candidate())]]
        for tasks in variants:
            with self.subTest(tasks=tasks), self.assertRaises(ValueError):
                pool.prepare([(1, '1/2')], 199229, 33, 384, tasks)

    def test_fresh_worker_calls_ordinary_outward_and_checks_geometry(self):
        scope, tasks = prepared()
        model = ToyModel()
        result = pool.evaluate(model, scope, tasks[0])
        self.assertEqual(model.seen, (tasks[0][1], tasks[0][2], 384))
        self.assertEqual(result['upper'], [1, -50])
        self.assertEqual(result['cell'], ['0', '1/2'])
        with self.assertRaises(ArithmeticError):
            pool.evaluate(model, scope, ('0', (Q(0), Q(1)), candidate()))

    def test_precision_change_nonfinite_and_nonpositive_fail_closed(self):
        scope, tasks = prepared()
        for value in (arb(0), arb(-1), arb('nan'), arb('inf')):
            model = SimpleNamespace(root=(Q(0), Q(1)), outward=lambda *args: value)
            with self.assertRaises(ArithmeticError):
                pool.evaluate(model, scope, tasks[0])
        def changed(*args):
            ctx.prec = 192
            return arb(1)
        with self.assertRaises(ArithmeticError):
            pool.evaluate(SimpleNamespace(root=(Q(0), Q(1)), outward=changed), scope, tasks[0])

    def test_exact_scope_coverage_returns_rows_without_new_claim(self):
        scope, tasks = prepared()
        result = pool.validate_results(scope, tasks, self.fresh_rows(scope, tasks))
        self.assertEqual(len(result), 2)
        self.assertTrue(all(row['upper'] == [1, -50] for row in result))
        self.assertEqual([row['path'] for row in result], ['0', '1'])

    def test_parent_rejects_any_scope_or_witness_mismatch(self):
        scope, tasks = prepared()
        corruptions = [dict(path='1'), dict(cell=['0', '1']), dict(precision=256),
            dict(threshold=199230), dict(minimum_groups=34), dict(scope_sha256='x'),
            dict(witness_sha256='x'), dict(output_tilt='1/99'), dict(upper=[0, 0]),
            dict(upper=[True, -50]), dict(upper=[1, -50.0])]
        for changes in corruptions:
            rows = self.fresh_rows(scope, tasks)
            rows[0].update(changes)
            with self.subTest(changes=changes), self.assertRaises((ArithmeticError, ValueError)):
                pool.validate_results(scope, tasks, rows)
        one_scope, one_job = pool.prepare([(1, '1/2')], 1, 1, 384,
            [('', (0, 1), candidate())])
        rows = self.fresh_rows(one_scope, one_job)
        rows[0]['threshold'] = True
        with self.assertRaises(ArithmeticError):
            pool.validate_results(one_scope, one_job, rows)

    def test_parent_rejects_missing_extra_duplicate_and_reordered_results(self):
        scope, tasks = prepared()
        rows = self.fresh_rows(scope, tasks)
        for bad in (rows[:1], rows+[rows[0]], [rows[0], rows[0]], rows[::-1]):
            with self.assertRaises(ArithmeticError):
                pool.validate_results(scope, tasks, bad)

    def test_worker_failure_propagates_without_partial_success(self):
        scope, tasks = prepared()
        def failing_rows():
            yield pool.evaluate(ToyModel(), scope, tasks[0])
            raise RuntimeError('worker arithmetic failed')
        with self.assertRaisesRegex(RuntimeError, 'worker arithmetic failed'):
            pool.validate_results(scope, tasks, failing_rows())

    def test_model_construction_uses_fresh_actual_r2_and_fixed_scope(self):
        scope, _ = prepared()
        original_path = sys.path[:]
        self.addCleanup(lambda: sys.path.__setitem__(slice(None), original_path))
        data = {'fresh': True}
        def actual(updates):
            self.assertEqual(updates, 2)
            ctx.prec = 192
            return data
        birth = SimpleNamespace(actual=Mock(side_effect=actual))
        mixture = SimpleNamespace(as_components=Mock(side_effect=lambda rows: ('fresh', rows)))
        class PlainModel:
            def __init__(self, *args, **kwargs):
                self.args, self.kwargs, self.precision = args, kwargs, ctx.prec
        with patch.dict(sys.modules, birth_classes=birth, shared_mixture=mixture,
                scalar_cover=SimpleNamespace(Model=PlainModel)):
            model = pool._build_model(scope)
        self.assertIs(type(model), PlainModel)
        self.assertEqual(model.precision, 384)
        self.assertIs(model.args[1], data)
        self.assertEqual(model.args[2:], (199229, 33, Q(3, 16)))
        self.assertEqual(model.kwargs, dict(inner=birth, variance_shuffle=True,
            variance_bins=16, regional_count=True))
        birth.actual.assert_called_once_with(2)
        mixture.as_components.assert_called_once()

    def test_one_model_per_initializer_not_per_cell(self):
        scope, tasks = prepared()
        with patch.object(pool, '_build_model', return_value=ToyModel()) as build:
            pool._initialize(scope, None)
            rows = [pool._evaluate_worker(task) for task in tasks]
        build.assert_called_once_with(scope)
        self.assertEqual([row['path'] for row in rows], ['0', '1'])

    def test_spawn_bounded_workers_determinism_and_no_actual_processes(self):
        scope, tasks = prepared()
        executor = Mock()
        executor.map.side_effect = lambda fn, tasks, chunksize: iter(self.fresh_rows(scope, tasks))
        constructor = Mock()
        constructor.return_value.__enter__ = Mock(return_value=executor)
        constructor.return_value.__exit__ = Mock(return_value=False)
        original_path = sys.path[:]
        with patch.object(pool, 'ProcessPoolExecutor', constructor), patch.object(pool, 'get_context') as context:
            result = pool.replay_dense(scope['mixture'], 199229, 33, 384, jobs(), workers=3)
        context.assert_called_once_with('spawn')
        self.assertEqual(constructor.call_args.kwargs['max_workers'], 3)
        self.assertIs(constructor.call_args.kwargs['initializer'], pool._initialize)
        self.assertEqual(executor.map.call_args.kwargs, dict(chunksize=1))
        self.assertEqual([row['upper'] for row in result], [[1, -50], [1, -50]])
        self.assertEqual(sys.path, original_path)

    def test_worker_limits_and_existing_logs_rejected_before_launch(self):
        with patch.object(pool, 'ProcessPoolExecutor') as launch:
            for workers in (0, 5, True, 2.0):
                with self.assertRaises(ValueError):
                    pool.replay_dense([(1, '1/2')], 1, 33, 384, jobs(), workers=workers)
            with TemporaryDirectory() as directory:
                with self.assertRaises(ValueError):
                    pool.replay_dense([(1, '1/2')], 1, 33, 384, jobs(), logdir=Path(directory))
            launch.assert_not_called()


if __name__ == '__main__':
    unittest.main()
