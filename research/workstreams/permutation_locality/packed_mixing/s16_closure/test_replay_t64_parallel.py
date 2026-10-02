"""Tiny/mock worker tests; no production model preparation or cover replay."""
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor as RealProcessPoolExecutor
from fractions import Fraction as Q
from multiprocessing import get_context
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from flint import arb, ctx
import replay_t64_parallel as replay
import test_whole_t64


class ToyModel:
    def __init__(self, root):
        self.root = tuple(map(Q, root))

    def outward(self, cell, witness):
        self.seen = cell, witness, ctx.prec
        return arb(2)**-60


def tiny_spawn_row(scope, job):
    # Real process/import/pickle smoke test, without a production initializer.
    return replay.evaluate(ToyModel(scope['root']), scope, 192, job)


class ParallelReplayTests(unittest.TestCase):
    # Reuse only the synthetic fixture, not its entire inherited test suite.
    setUp = test_whole_t64.WholeT64Tests.setUp
    save = test_whole_t64.WholeT64Tests.save
    prepared = test_whole_t64.WholeT64Tests.prepared
    fresh = test_whole_t64.WholeT64Tests.fresh
    mock_manifest = test_whole_t64.WholeT64Tests.mock_manifest

    def worker(self, scope, proposal, pid=123):
        return dict(pid=pid, precision=384, scope_sha256=replay.dense_t64.prior.fingerprint(scope),
            source=proposal, initialization='fresh dense_t64.fresh_model in spawned process')

    def setup_jobs(self):
        return replay.prepare(self.dp, self.proposal_path, 384, 2)

    def fake_executor(self, *, edit=None, omit=False, extra=False):
        case = self
        class FakeExecutor:
            def __init__(self, **kwargs):
                case.assertEqual(kwargs['max_workers'], 2)
                case.assertEqual(kwargs['mp_context'].get_start_method(), 'spawn')
                case.assertIs(kwargs['initializer'], replay._initialize)
                self.proposal, self.scope, self.precision = kwargs['initargs']

            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

            def map(self, function, jobs, chunksize):
                case.assertIs(function, replay._evaluate_worker)
                case.assertEqual(chunksize, 1)
                model = ToyModel(self.scope['root'])
                messages = [dict(row=replay.evaluate(model, self.scope, self.precision, job),
                    worker=case.worker(self.scope, self.proposal, 123+i%2)) for i, job in enumerate(jobs)]
                if edit:
                    edit(messages)
                if omit:
                    messages.pop()
                if extra:
                    messages.append(deepcopy(messages[0]))
                return iter(messages)
        return FakeExecutor

    def test_initializer_authenticates_fresh_full_scope_source_and_precision(self):
        old, _, proposal, _, _ = self.setup_jobs()
        scope = old['scope']
        model = ToyModel(scope['root'])
        ctx.prec = 384
        with patch.object(replay.dense_t64, 'fresh_model', return_value=(model, deepcopy(scope), proposal)) as fresh:
            replay._initialize(proposal, scope, 384)
        fresh.assert_called_once_with(proposal['path'], precision=384, variance_bins=16,
            distance=Q(1,10), birth_density='capped')
        self.assertIs(replay._model, model)
        self.assertEqual(replay._worker['source'], proposal)
        for change in (dict(physical_t=128), dict(root=['0', '1']), dict(maps={})):
            with patch.object(replay.dense_t64, 'fresh_model', return_value=(model, dict(scope, **change), proposal)), \
                 self.assertRaises(ValueError):
                replay._initialize(proposal, scope, 384)
        with patch.object(replay.dense_t64, 'fresh_model', return_value=(model, scope, dict(proposal, sha256='x'))), \
             self.assertRaises(ValueError):
            replay._initialize(proposal, scope, 384)

    def test_worker_outward_frontdoor_precision_and_geometry(self):
        old, _, _, _, jobs = self.setup_jobs()
        model, scope = ToyModel(old['scope']['root']), old['scope']
        row = replay.evaluate(model, scope, 384, jobs[0])
        self.assertEqual(set(row), {'path', 'cell', 'upper'})
        self.assertEqual(row['upper'], [1, -60])
        self.assertEqual(model.seen, (jobs[0][1], jobs[0][2], 384))
        with self.assertRaises(ValueError):
            replay.evaluate(model, scope, 384, (jobs[0][0], (Q(0), Q(1)), {}))
        for invalid in (arb(0), arb(-1), arb('nan'), arb('inf')):
            bad = SimpleNamespace(root=model.root, outward=lambda *args: invalid)
            with self.assertRaises(ArithmeticError):
                replay.evaluate(bad, scope, 384, jobs[0])
        def changed(*args):
            ctx.prec = 192
            return arb(1)
        with self.assertRaises(ValueError):
            replay.evaluate(SimpleNamespace(root=model.root, outward=changed), scope, 384, jobs[0])

    def test_result_rejects_wrong_cell_precision_source_scope_and_duplicate(self):
        old, _, proposal, cells, jobs = self.setup_jobs()
        scope = old['scope']
        row = replay.evaluate(ToyModel(scope['root']), scope, 384, jobs[0])
        message = dict(row=row, worker=self.worker(scope, proposal))
        replay.check_result(message, cells, set(), scope, proposal, 384)
        edits = [lambda m: m['row'].update(path='11'), lambda m: m['row'].update(cell=['0', '1']),
            lambda m: m['row'].update(upper=[0,0]), lambda m: m['row'].update(saved_score=0),
            lambda m: m['worker'].update(precision=256), lambda m: m['worker'].update(scope_sha256='x'),
            lambda m: m['worker'].update(source={}), lambda m: m['worker'].update(pid=True),
            lambda m: m['worker'].update(initialization='inherited fork')]
        for edit in edits:
            altered = deepcopy(message)
            edit(altered)
            with self.subTest(edit=edit), self.assertRaises(ValueError):
                replay.check_result(altered, cells, set(), scope, proposal, 384)
        with self.assertRaises(ValueError):
            replay.check_result(message, cells, {row['path']}, scope, proposal, 384)

    def test_parent_normal_dense_schema_exact_sum_and_spawn_metadata(self):
        output = self.directory/'parallel.json'
        with patch.object(replay, 'ProcessPoolExecutor', self.fake_executor()), \
             patch.object(replay.whole, 'proof_source_manifest', side_effect=self.mock_manifest):
            result = replay.run(self.dp, self.proposal_path, output, workers=2)
        self.assertEqual(result['schema'], 'packed-gl32-t64-s16-dense-replay-1')
        self.assertTrue(result['complete_dense'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['workers'], 2)
        self.assertEqual(result['worker_start_method'], 'spawn')
        self.assertEqual(len(result['worker_metadata']), 2)
        self.assertEqual(replay.whole.dyadic(result['aggregate_upper']), Q(1, 1 << 59))
        self.assertEqual(result['source_manifest_before'], result['source_manifest_after'])
        sparse, _ = self.fresh(self.prepared())
        whole = replay.whole._assemble_fresh(self.prepared(), sparse, result, 384, 40)
        self.assertTrue(whole['verified'])

    def test_missing_duplicate_extra_or_mutated_source_never_complete(self):
        for name, executor in [('missing', self.fake_executor(omit=True)),
                ('extra', self.fake_executor(extra=True)),
                ('duplicate', self.fake_executor(edit=lambda rows: rows.__setitem__(1, deepcopy(rows[0]))))]:
            path = self.directory/(name+'.json')
            with patch.object(replay, 'ProcessPoolExecutor', executor), \
                 patch.object(replay.whole, 'proof_source_manifest', side_effect=self.mock_manifest), \
                 self.assertRaises(ValueError):
                replay.run(self.dp, self.proposal_path, path, workers=2)
            record, _ = replay.whole.read_source(path)
            self.assertFalse(record['complete_dense'])
        def mutate(rows):
            self.map_path.write_text('changed')
        with patch.object(replay, 'ProcessPoolExecutor', self.fake_executor(edit=mutate)), \
             patch.object(replay.whole, 'proof_source_manifest', side_effect=self.mock_manifest), \
             self.assertRaises(ValueError):
            replay.run(self.dp, self.proposal_path, self.directory/'mutation.json', workers=2)

    def test_bad_workers_partial_search_and_stale_output_prevent_launch(self):
        for workers in (0, 5, True, 2.0):
            with patch.object(replay, 'ProcessPoolExecutor') as executor, self.assertRaises(ValueError):
                replay.run(self.dp, self.proposal_path, self.directory/'none.json', workers=workers)
            executor.assert_not_called()
        with patch.object(replay, 'ProcessPoolExecutor') as executor, self.assertRaises(ValueError):
            replay.run(self.dp, self.proposal_path, self.dp, workers=2)
        executor.assert_not_called()
        self.dense['cover']['leaves'].pop('1')
        self.save()
        with patch.object(replay, 'ProcessPoolExecutor') as executor, self.assertRaises(ValueError):
            replay.run(self.dp, self.proposal_path, self.directory/'none.json', workers=2)
        executor.assert_not_called()

    def test_cli_has_only_supported_parallel_flags(self):
        arguments = ['parallel', '--input', str(self.dp), '--source', str(self.proposal_path),
            '--output', str(self.directory/'parallel.json'), '--workers', '3', '--precision', '384']
        with patch.object(replay.sys, 'argv', arguments), patch.object(replay, 'run') as run:
            replay.main()
        run.assert_called_once_with(self.dp, self.proposal_path, self.directory/'parallel.json', 384, 3)

    def test_real_spawn_tiny_outward_smoke(self):
        with RealProcessPoolExecutor(max_workers=1, mp_context=get_context('spawn')) as executor:
            row = executor.submit(tiny_spawn_row, dict(root=['0', '1']),
                ('0', (Q(0), Q(1,2)), {})).result(timeout=60)
        self.assertEqual(row, dict(path='0', cell=['0', '1/2'], upper=[1, -60]))


if __name__ == '__main__':
    unittest.main()
