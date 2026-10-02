"""Synthetic orchestration tests only; child proof computations are mocked."""
import copy
from contextlib import redirect_stdout
from fractions import Fraction as Q
import io
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import whole_hill_replay as whole


def fixture(updates=4):
    premises = dict(schema=whole.hill_cover.EKR_SCHEMA,
        canonical_cdf=[0]*32+[(1 << 512)-1], block_width=8, block_rows=4, groups_per_outer_word=32)
    dense_hash, sparse_hash, floor = whole.count_fingerprints(premises)
    scope = dict(schema=whole.hill_cover.CONTEXT_SCHEMA,
        ensemble=f'canonical-gl32-width8-shared4-r{updates}', K=1 << 20, N=1 << 21,
        updates=updates, block_width=8, distance='1/10', threshold=209715,
        minimum_groups=33, maximum_groups=2048, root=['0', '1'],
        comparison='direct-expected-shell-majorant', last_lp=104, refined=True,
        base_tilt='3/16', variance_bins=16, regional_count=True, outer_premises=premises,
        expected_cdf_sha256=dense_hash, comparison_caps_sha256='b'*64,
        mixture=[dict(mass='1', activity='1/2')])
    dense = dict(schema=whole.hill_cover.SCHEMA, scope=scope, precision=256, cell_target_bits=52,
        whole_code_certificate=False, cover=dict(unresolved={}, visited=2,
            leaves={p: dict(cell=cell, witness={}, upper='OLD BOUND IS NOT INPUT')
                    for p, cell in [('0', ['0', '1/2']), ('1', ['1/2', '1'])]}))
    sparse = dict(schema='packed-gl32-sparse-hill-variant-1', K=1 << 20, N=1 << 21,
        group_count=2048, outer='BCH256128', block_rows=4, block_columns=8,
        mixing='independent uniform GL32 per group/block; no additional GF16 stage',
        routing='independent shared column shuffle per group and independent regional shuffles',
        inner=dict(t=128, s=19, updates=updates), support_min=floor,
        distance='1/10', threshold=209715, count_refinement='ekr',
        count_premises=premises, count_sha256=sparse_hash,
        joint_return_through=2, lazy_density_through=4, tilts=['.001'],
        requested=list(range(1, 33)), complete_requested=True, complete_sparse_prefix=True,
        uncovered_sparse=[], results=[])
    for q in range(1, 33):
        details = (dict(method='all-support exact placement', support_choices=['.001']*257)
                   if q == 1 else dict(method='full support-box CDF cover', leaves=[dict(tilt='.001')]))
        details['support_min'] = floor
        sparse['results'].append(dict(occupancy=q, passed=True, upper='OLD BOUND IS NOT INPUT', details=details))
    return dense, sparse


def fresh(prepared):
    sparse = dict(schema='packed-gl32-sparse-hill-variant-replay-1', precision=384, target_bits=41,
        distance='1/10', threshold=209715, scope=prepared['construction'],
        sources=[prepared['sparse_source']], count_premises=prepared['count_premises'],
        count_sha256=prepared['sparse_count_sha256'], count_refinement='ekr',
        joint_return_through=2, lazy_density_through=4,
        requested=list(range(1, 33)), uncovered_sparse=[], complete_requested=True,
        complete_sparse_prefix=True, aggregate_upper=[1, -55],
        results=[dict(occupancy=q, passed=True, upper=[1, -60]) for q in range(1, 33)])
    dense = copy.deepcopy(prepared['dense'])
    dense.update(precision=384, source=prepared['dense_source'])
    dense['fresh_replay'] = dict(precision=384, target_bits=41, scope=prepared['scope'],
        complete_dense=True, passed=True, whole_code_certificate=False, aggregate_upper=[1, -59],
        checked=[dict(path=p, cell=row['cell'], upper=[1, -60])
                 for p, row in dense['cover']['leaves'].items()])
    return sparse, dense


class WholeHillReplayTests(unittest.TestCase):
    def setUp(self):
        quiet = redirect_stdout(io.StringIO())
        quiet.__enter__()
        self.addCleanup(quiet.__exit__, None, None, None)
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.dp, self.sp = self.directory/'dense.json', self.directory/'sparse.json'
        self.dense, self.sparse = fixture()
        self.save()

    def save(self):
        self.dp.write_text(json.dumps(self.dense))
        self.sp.write_text(json.dumps(self.sparse))

    def prepared(self):
        return whole.prepare(self.dp, self.sp)

    def test_native_hash_formats_share_one_exact_cdf(self):
        prepared = self.prepared()
        self.assertNotEqual(prepared['dense_count_sha256'], prepared['sparse_count_sha256'])
        self.assertEqual(prepared['scope']['updates'], 4)
        sparse, dense = fresh(prepared)
        result = whole._assemble_fresh(prepared, sparse, dense, 384, 40)
        self.assertEqual(whole.dyadic(result['upper']), Q(34, 1 << 60))
        self.assertTrue(result['verified'])

    def test_scope_cutoff_counts_and_completeness_rejected(self):
        edits = [lambda d, s: s['inner'].update(updates=3),
            lambda d, s: s.update(distance='11/100', threshold=int(Q(11, 100)*(1 << 21))),
            lambda d, s: s.update(count_refinement='incidence'),
            lambda d, s: s.update(count_sha256='d'*64),
            lambda d, s: s['count_premises'].update(extra=True),
            lambda d, s: s['results'].pop(),
            lambda d, s: s['results'].__setitem__(0, copy.deepcopy(s['results'][1])),
            lambda d, s: s.update(complete_sparse_prefix=False),
            lambda d, s: d['scope'].update(minimum_groups=34),
            lambda d, s: d['cover']['leaves'].pop('1')]
        original_dense, original_sparse = copy.deepcopy(self.dense), copy.deepcopy(self.sparse)
        for edit in edits:
            self.dense, self.sparse = copy.deepcopy(original_dense), copy.deepcopy(original_sparse)
            edit(self.dense, self.sparse)
            self.save()
            with self.subTest(edit=edit), self.assertRaises((ValueError, KeyError)):
                self.prepared()

    def test_fresh_source_rows_and_geometry_fail_closed(self):
        prepared = self.prepared()
        edits = [lambda s, d: s.update(sources=[]),
            lambda s, d: d.update(source={}),
            lambda s, d: s['scope']['inner'].update(updates=3),
            lambda s, d: s.update(precision=256),
            lambda s, d: s.update(threshold=1),
            lambda s, d: s.update(count_sha256='f'*64),
            lambda s, d: s['results'].pop(),
            lambda s, d: s['results'].__setitem__(0, dict(occupancy=2, passed=True, upper=[1, -60])),
            lambda s, d: d['fresh_replay'].update(complete_dense=False),
            lambda s, d: d['fresh_replay']['checked'].pop(),
            lambda s, d: d['fresh_replay']['checked'][0].update(cell=['0', '1']),
            lambda s, d: d['fresh_replay']['checked'][0].update(path='1')]
        for edit in edits:
            sparse, dense = copy.deepcopy(fresh(prepared))
            edit(sparse, dense)
            with self.subTest(edit=edit), self.assertRaises(ValueError):
                whole._assemble_fresh(prepared, sparse, dense, 384, 40)

    def test_aggregate_consistency_strict_goal_and_compaction(self):
        prepared = self.prepared()
        for which in ('sparse', 'dense'):
            sparse, dense = fresh(prepared)
            if which == 'sparse':
                sparse['aggregate_upper'] = [1, -60]
            else:
                dense['fresh_replay']['aggregate_upper'] = [1, -60]
            with self.assertRaises(ValueError):
                whole._assemble_fresh(prepared, sparse, dense, 384, 40)
        sparse, dense = fresh(prepared)
        for row in sparse['results']:
            row['upper'] = [1, -44]
        sparse['aggregate_upper'] = [1, -39]
        with self.assertRaises(ValueError):
            whole._assemble_fresh(prepared, sparse, dense, 384, 40)
        value = Q((1 << 700)+3, 1 << 900)
        compact = whole.compact_endpoint(value, 64)
        self.assertGreaterEqual(whole.dyadic(compact), value)
        self.assertLessEqual(compact[0].bit_length(), 65)

    def test_serial_children_and_existing_output_rejection(self):
        prepared = self.prepared()
        sparse, dense = fresh(prepared)
        commands = []
        def child(command, *, cwd, env, check):
            commands.append(command)
            self.assertEqual(cwd, whole.REPO)
            self.assertEqual(env['OMP_NUM_THREADS'], '1')
            self.assertTrue(check)
            output = Path(command[command.index('--output')+1])
            self.assertFalse(output.exists())
            output.write_text(json.dumps(sparse if len(commands) == 1 else dense))
        output = self.directory/'new'
        with patch.object(whole.subprocess, 'run', side_effect=child):
            result = whole.run(self.dp, self.sp, output)
        self.assertEqual([Path(command[1]).name for command in commands],
            ['sparse_hill_variant.py', 'hill_cover.py'])
        self.assertTrue(result['verified'])
        self.assertTrue((output/'whole.json').is_file())
        self.assertEqual(result['proof_sources'], whole.proof_source_manifest())
        self.assertTrue(any(row['path'].endswith('BchCircuit.h')
                            for row in result['proof_sources']['files']))
        self.assertTrue(any(row['path'].endswith('SelectedMaps.h')
                            for row in result['proof_sources']['files']))
        with patch.object(whole.subprocess, 'run') as child, self.assertRaises(ValueError):
            whole.run(self.dp, self.sp, output)
        child.assert_not_called()

    def test_failed_missing_or_mutated_child_never_writes_whole(self):
        output = self.directory/'failed'
        with patch.object(whole.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'test')), \
                self.assertRaises(subprocess.CalledProcessError):
            whole.run(self.dp, self.sp, output)
        self.assertFalse((output/'whole.json').exists())
        with patch.object(whole.subprocess, 'run'), self.assertRaises(ValueError):
            whole.run(self.dp, self.sp, self.directory/'missing')
        prepared = self.prepared()
        sparse, dense = fresh(prepared)
        def child(command, **kwargs):
            output = Path(command[command.index('--output')+1])
            output.write_text(json.dumps(sparse if Path(command[1]).name == 'sparse_hill_variant.py' else dense))
            self.dp.write_text(self.dp.read_text()+'\n')
        output = self.directory/'mutated'
        with patch.object(whole.subprocess, 'run', side_effect=child), self.assertRaises(ValueError):
            whole.run(self.dp, self.sp, output)
        self.assertFalse((output/'whole.json').exists())

    def test_proof_source_mutations_fail_before_whole_claim(self):
        prepared = self.prepared()
        sparse, dense = fresh(prepared)
        for mutation in ('python', 'header', 'added', 'removed'):
            with self.subTest(mutation=mutation):
                root = self.directory/mutation
                here = root/'research/workstreams/permutation_locality/packed_mixing'
                here.mkdir(parents=True)
                source = here/'proof.py'
                source.write_text('original = True\n')
                generated = root/'spin/src/kernels/generated'
                generated.mkdir(parents=True)
                for name in ('BchCircuit.h', 'SelectedMaps.h'):
                    (generated/name).write_text('// original\n')
                calls = []
                def child(command, **kwargs):
                    calls.append(command)
                    output = Path(command[command.index('--output')+1])
                    output.write_text(json.dumps(sparse if len(calls) == 1 else dense))
                    if len(calls) == 2:
                        if mutation == 'python':
                            source.write_text('original = False\n')
                        elif mutation == 'header':
                            (generated/'SelectedMaps.h').write_text('// changed\n')
                        elif mutation == 'added':
                            (here/'new_proof.py').write_text('# newly added\n')
                        else:
                            source.unlink()
                output = root/'out'
                with patch.object(whole, 'HERE', here), patch.object(whole, 'REPO', root), \
                        patch.object(whole.subprocess, 'run', side_effect=child), \
                        self.assertRaises(ValueError):
                    whole.run(self.dp, self.sp, output)
                self.assertEqual(len(calls), 2)
                self.assertFalse((output/'whole.json').exists())


if __name__ == '__main__':
    unittest.main()
