"""Synthetic scope/aggregation tests; no numerical proof jobs are launched."""
import copy
from contextlib import redirect_stdout
from fractions import Fraction as Q
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import whole_replay as whole


def fixture():
    premises = dict(schema='test-only-authenticated-by-mocked-children', block_width=8)
    dense = dict(schema=whole.dense_cover.SCHEMA, ensemble=whole.dense_cover.ENSEMBLE,
        K=1 << 20, N=1 << 21, updates=2, block_width=8, minimum_groups=33,
        maximum_groups=2048, distance='19/200', threshold=int(Q(19, 200)*(1 << 21)),
        comparison='direct-expected-shell-majorant', last_lp=104, refined=True,
        base_tilt='3/16', variance_bins=16, regional_count=True,
        root=['1/1000', '1'], outer_premises=premises,
        expected_cdf_sha256='a'*64, comparison_caps_sha256='b'*64,
        cover=dict(leaves={'0': dict(witness={}), '1': dict(witness={})}, unresolved={}))
    # Saved numeric endpoints are intentionally nonsensical: only fresh
    # child endpoints may contribute to the whole-code sum.
    sparse = dict(whole.SPARSE_SCOPE, schema='packed-gl32-sparse-1',
        distance='12/125', threshold=int(Q(12, 125)*(1 << 21)),
        count_premises=premises, count_sha256='c'*64,
        results=[dict(occupancy=q, passed=True, upper='NOT A PROOF INPUT', details={})
                 for q in range(1, 33)])
    return dense, sparse


def fresh(prepared, precision=384):
    sparse = dict(schema='packed-gl32-sparse-replay-1', scope=whole.SPARSE_SCOPE,
        distance=str(prepared['distance']), threshold=prepared['threshold'], precision=precision,
        sources=prepared['sparse_sources'], count_premises=prepared['premises'], count_sha256='c'*64,
        requested=list(range(1, 33)), uncovered_sparse=[], complete_requested=True,
        complete_sparse_prefix=True, results=[dict(occupancy=q, upper=[1, -40]) for q in range(1, 33)],
        aggregate_upper=[1, -35])
    dense = copy.deepcopy(prepared['dense'])
    dense.update(precision=precision, source=prepared['dense_source'])
    cells = whole.partition_cells(dense['root'], dense['cover']['leaves'])
    dense['fresh_replay'] = dict(complete_dense=True, unresolved=0,
        checked=[dict(path=p, cell=list(map(str, cell)), upper=[1, -40]) for p, cell in cells.items()],
        upper=[1, -39])
    return sparse, dense


class WholeReplayTests(unittest.TestCase):
    def setUp(self):
        # Mocked child runs must not print a synthetic certificate claim
        # into a transcript beside real proof-worker output.
        self.quiet = redirect_stdout(io.StringIO())
        self.quiet.__enter__()
        self.addCleanup(self.quiet.__exit__, None, None, None)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.dense, self.sparse = fixture()
        self.dp, self.sp = self.directory/'dense-source.json', self.directory/'sparse-source.json'
        self.save_sources()

    def save_sources(self):
        self.dp.write_text(json.dumps(self.dense), encoding='utf-8')
        self.sp.write_text(json.dumps(self.sparse), encoding='utf-8')

    def prepared(self):
        return whole.prepare(self.dp, [self.sp])

    def test_exact_sum_and_old_numbers_ignored(self):
        prepared = self.prepared()
        sparse, dense = fresh(prepared)
        result = whole._assemble_fresh(prepared, sparse, dense, 384, 20)
        self.assertEqual(whole.dyadic(result['upper']), Q(34, 1 << 40))
        self.assertTrue(result['verified'])
        self.assertEqual(result['threshold'], self.dense['threshold'])
        self.assertLess(Q(self.dense['distance']), Q(self.sparse['distance']))

    def test_input_dense_partition_rejections(self):
        for leaves, unresolved in (({'0': dict(witness={})}, {}),
                                   ({'': dict(witness={}), '0': dict(witness={})}, {}),
                                   ({'0': dict(witness={})}, {'1': {}}), ({}, {})):
            with self.subTest(leaves=leaves, unresolved=unresolved):
                trial = copy.deepcopy(self.dense)
                trial['cover'] = dict(leaves=leaves, unresolved=unresolved)
                with self.assertRaises(ValueError):
                    whole.complete_dense(trial)
        self.dense['minimum_groups'] = 32
        self.save_sources()
        with self.assertRaises(ValueError):
            self.prepared()

    def test_missing_and_duplicate_sparse_scopes(self):
        self.sparse['results'].pop()
        self.save_sources()
        with self.assertRaises(ValueError):
            self.prepared()
        self.dense, self.sparse = fixture()
        self.save_sources()
        other = self.directory/'overlap.json'
        other.write_text(json.dumps(self.sparse), encoding='utf-8')
        with self.assertRaises(ValueError):
            whole.prepare(self.dp, [self.sp, other])
        self.sparse['results'].append(copy.deepcopy(self.sparse['results'][0]))
        self.save_sources()
        with self.assertRaises(ValueError):
            self.prepared()

    def test_mismatched_source_construction_cutoff_and_counts(self):
        for change in (dict(block_columns=16), dict(threshold=5),
                       dict(count_premises={'different': True}),
                       dict(distance='9/100', threshold=int(Q(9, 100)*(1 << 21)))):
            with self.subTest(change=change):
                _, self.sparse = fixture()
                self.sparse.update(change)
                self.save_sources()
                with self.assertRaises(ValueError):
                    self.prepared()

    def test_fresh_scope_and_provenance_rejections(self):
        prepared = self.prepared()
        edits = [
            lambda s, d: s.update(complete_sparse_prefix=False),
            lambda s, d: s['results'].pop(),
            lambda s, d: s['results'].__setitem__(0, dict(occupancy=2, upper=[1, -40])),
            lambda s, d: s.update(requested=list(range(2, 33))),
            lambda s, d: s.update(sources=[]),
            lambda s, d: d.update(source=dict(path='elsewhere', sha256='0'*64)),
            lambda s, d: s.update(count_premises={'different': True}),
            lambda s, d: s.update(count_sha256='d'*64),
            lambda s, d: d.update(expected_cdf_sha256='d'*64),
            lambda s, d: s.update(threshold=s['threshold']+1),
            lambda s, d: d['fresh_replay'].update(complete_dense=False),
            lambda s, d: d['fresh_replay'].update(unresolved=1),
            lambda s, d: d['fresh_replay']['checked'].pop(),
            lambda s, d: d['fresh_replay']['checked'].append(copy.deepcopy(d['fresh_replay']['checked'][0])),
            lambda s, d: d['fresh_replay']['checked'][0].update(path='00'),
            lambda s, d: d['fresh_replay']['checked'][0].update(cell=['0', '1']),
        ]
        for edit in edits:
            sparse, dense = fresh(prepared)
            edit(sparse, dense)
            with self.subTest(edit=edit), self.assertRaises(ValueError):
                whole._assemble_fresh(prepared, sparse, dense, 384, 20)

    def test_aggregate_consistency_and_strict_twenty_bit_goal(self):
        prepared = self.prepared()
        for which in ('sparse', 'dense'):
            sparse, dense = fresh(prepared)
            if which == 'sparse':
                sparse['aggregate_upper'] = [1, -40]
            else:
                dense['fresh_replay']['upper'] = [1, -40]
            with self.assertRaises(ValueError):
                whole._assemble_fresh(prepared, sparse, dense, 384, 20)
        sparse, dense = fresh(prepared)
        # Each sparse term passes20bits, but their sum does not.
        for row in sparse['results']:
            row['upper'] = [1, -24]
        sparse['aggregate_upper'] = [1, -19]
        with self.assertRaises(ValueError):
            whole._assemble_fresh(prepared, sparse, dense, 384, 20)
        # Equality at2^-20 must also fail the strict margin.
        for row in sparse['results']:
            row['upper'] = [1, -26]
        sparse['aggregate_upper'] = [1, -21]
        for row in dense['fresh_replay']['checked']:
            row['upper'] = [1, -22]
        dense['fresh_replay']['upper'] = [1, -21]
        with self.assertRaises(ValueError):
            whole._assemble_fresh(prepared, sparse, dense, 384, 20)
        with self.assertRaises(ValueError):
            whole._assemble_fresh(prepared, *fresh(prepared), 384, 19)

    def test_compact_endpoint_rounds_up(self):
        value = Q((1 << 700)+3, 1 << 900)
        endpoint = whole.compact_endpoint(value, 64)
        self.assertGreaterEqual(whole.dyadic(endpoint), value)
        self.assertLessEqual(endpoint[0].bit_length(), 65)
        for endpoint in ([0, -2], [-1, -2], [True, -2], [1, 0.5]):
            with self.assertRaises(ValueError):
                whole.dyadic(endpoint)

    def test_serial_children_produce_only_fresh_final_receipt(self):
        prepared = self.prepared()
        sparse, dense = fresh(prepared)
        names = []
        def worker(command, *, cwd, check):
            name = Path(command[1]).name
            names.append(name)
            self.assertTrue(check)
            self.assertEqual(cwd, whole.REPO)
            output = Path(command[command.index('--output')+1])
            self.assertFalse(output.exists())
            output.write_text(json.dumps(sparse if name == 'sparse.py' else dense), encoding='utf-8')
        output = self.directory/'fresh'
        with patch.object(whole.subprocess, 'run', side_effect=worker):
            result = whole.run(self.dp, [self.sp], output)
        self.assertEqual(names, ['sparse.py', 'dense_cover.py'])
        self.assertTrue((output/'whole.json').exists())
        self.assertTrue(result['verified'])
        with patch.object(whole.subprocess, 'run') as child:
            with self.assertRaises(ValueError):
                whole.run(self.dp, [self.sp], output)
            child.assert_not_called()

    def test_failed_child_never_produces_a_whole_claim(self):
        output = self.directory/'failed'
        with patch.object(whole.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'test')) as child:
            with self.assertRaises(subprocess.CalledProcessError):
                whole.run(self.dp, [self.sp], output)
            self.assertEqual(child.call_count, 1)
        self.assertFalse((output/'whole.json').exists())

    def test_parallel_dense_option_preserves_serial_stage_order(self):
        prepared = self.prepared()
        sparse, dense = fresh(prepared)
        commands = []
        def worker(command, **kwargs):
            commands.append(command)
            name = Path(command[1]).name
            output = Path(command[command.index('--output')+1])
            output.write_text(json.dumps(sparse if name == 'sparse.py' else dense), encoding='utf-8')
        with patch.object(whole.subprocess, 'run', side_effect=worker):
            result = whole.run(self.dp, [self.sp], self.directory/'parallel', dense_workers=3)
        self.assertEqual([Path(c[1]).name for c in commands], ['sparse.py', 'dense_cover.py'])
        self.assertNotIn('--replay-workers', commands[0])
        self.assertEqual(commands[1][-2:], ['--replay-workers', '3'])
        self.assertEqual(result['dense_replay_workers'], 3)
        for workers in (0, 5, True, 2.0):
            with patch.object(whole.subprocess, 'run') as child:
                with self.assertRaises(ValueError):
                    whole.run(self.dp, [self.sp], self.directory/'unused', dense_workers=workers)
                child.assert_not_called()

    def test_source_mutation_during_children_rejects_claim(self):
        prepared = self.prepared()
        sparse, dense = fresh(prepared)
        def worker(command, **kwargs):
            name = Path(command[1]).name
            output = Path(command[command.index('--output')+1])
            output.write_text(json.dumps(sparse if name == 'sparse.py' else dense), encoding='utf-8')
            if name == 'dense_cover.py':
                self.dp.write_text(self.dp.read_text()+'\n', encoding='utf-8')
        output = self.directory/'changed'
        with patch.object(whole.subprocess, 'run', side_effect=worker):
            with self.assertRaises(ValueError):
                whole.run(self.dp, [self.sp], output)
        self.assertFalse((output/'whole.json').exists())


if __name__ == '__main__':
    unittest.main()
