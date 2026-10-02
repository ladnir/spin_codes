"""Synthetic whole-replay orchestration; no numerical covers are run."""
from contextlib import redirect_stdout
from copy import deepcopy
from fractions import Fraction as Q
import io
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import whole_t64 as whole


class WholeT64Tests(unittest.TestCase):
    def setUp(self):
        quiet = redirect_stdout(io.StringIO())
        quiet.__enter__()
        self.addCleanup(quiet.__exit__, None, None, None)
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.dp, self.sp = self.directory/'dense-search.json', self.directory/'sparse-search.json'
        self.proposal_path, self.map_path = self.directory/'proposal.json', self.directory/'map.json'
        self.map_path.write_text('{"test_map": true}')
        self.map_record = dict(schema='synthetic-map-for-orchestration-only', t=64, s=16,
            source=whole.source_metadata(self.map_path), map_sha256='a'*64)
        self.map_patch = patch.object(whole, 'fresh_map_record', return_value=deepcopy(self.map_record))
        self.map_patch.start()
        self.addCleanup(self.map_patch.stop)
        premises = dict(schema=whole.sparse_t64.SCOPE['outer_count_refinement'],
            canonical_cdf=[0]*32+[(1 << 512)-1], block_width=8, block_rows=4, groups_per_outer_word=32)
        cdf_hash, shell_hash, sparse_hash, floor = whole.count_fingerprints(premises)
        minimum = Q(33, 2048)*Q(3, 19)
        midpoint = (minimum+1)/2
        scope = dict(whole.DENSE_SCOPE, distance='1/10', threshold=209715,
            root=[str(minimum), '1'], base_tilt='3/16', variance_bins=16,
            maps=deepcopy(self.map_record), maps_sha256=whole.dense_t64.prior.fingerprint(self.map_record),
            outer_premises=premises, expected_cdf_sha256=cdf_hash, comparison_caps_sha256=shell_hash,
            mixture=[dict(mass='1', activity='1/2')], mixture_verification=dict(synthetic=True))
        self.proposal_path.write_text(json.dumps(dict(scope=scope)))
        self.dense = dict(schema=whole.dense_t64.SCHEMA, scope=scope,
            source=whole.source_metadata(self.proposal_path), precision=256,
            whole_code_certificate=False, complete_search_partition=True,
            cover=dict(leaves={p: dict(cell=cell, witness={}, upper='OLD ENDPOINT NOT INPUT')
                for p, cell in [('0', [str(minimum), str(midpoint)]), ('1', [str(midpoint), '1'])]},
                unresolved={}, visited=3))
        self.sparse = dict(schema='t64-s16-uniform-gl-sparse-search-1',
            K=1 << 20, N=1 << 21, group_count=2048, outer='BCH256128', block_rows=4, block_columns=8,
            inner=deepcopy(whole.sparse_t64.INNER), kernel_variant=whole.sparse_t64.KERNEL_VARIANT,
            map_record=deepcopy(self.map_record), distance='1/10', threshold=209715,
            count_premises=deepcopy(premises), count_sha256=sparse_hash, support_min=floor,
            requested=list(range(1, 33)), tilts=['.001'], complete_requested=True,
            complete_sparse_prefix=True, whole_code_certificate=False, results=[], **whole.sparse_t64.SCOPE)
        for q in range(1, 33):
            details = (dict(method='all-support exact placement', support_choices=['.001']*257)
                       if q == 1 else dict(method='full support-box CDF cover', leaves=[dict(tilt='.001')]))
            details['support_min'] = floor
            self.sparse['results'].append(dict(occupancy=q, passed=True, upper='OLD ENDPOINT NOT INPUT', details=details))
        self.save()

    def save(self):
        self.dp.write_text(json.dumps(self.dense))
        self.sp.write_text(json.dumps(self.sparse))

    def prepared(self):
        return whole.prepare(self.dp, self.sp)

    def fresh(self, prepared):
        sparse = dict(schema='t64-s16-uniform-gl-sparse-replay-1',
            K=1 << 20, N=1 << 21, inner=deepcopy(whole.sparse_t64.INNER),
            kernel_variant=whole.sparse_t64.KERNEL_VARIANT, map_record=deepcopy(prepared['scope']['maps']),
            precision=384, target_bits=41, distance='1/10', threshold=209715,
            count_premises=deepcopy(prepared['count_premises']),
            count_sha256=prepared['sparse_count_sha256'], support_min=prepared['sparse']['support_min'],
            sources=[prepared['sparse_source']], requested=list(range(1, 33)),
            complete_requested=True, complete_sparse_prefix=True, whole_code_certificate=False,
            aggregate_upper=[1, -55], results=[dict(occupancy=q, passed=True, upper=[1, -60]) for q in range(1, 33)],
            **whole.sparse_t64.SCOPE)
        dense = dict(schema='packed-gl32-t64-s16-dense-replay-1', scope=deepcopy(prepared['scope']),
            source=prepared['proposal_source'], input_source=prepared['dense_source'], precision=384,
            complete_dense=True, whole_code_certificate=False, aggregate_upper=[1, -59],
            cells=[dict(path=path, cell=list(map(str, cell)), upper=[1, -60]) for path, cell in prepared['cells'].items()])
        return sparse, dense

    def mock_manifest(self, inputs):
        return dict(files=[whole.source_metadata(row['path']) for row in inputs])

    def child(self, prepared, edit=None):
        sparse, dense = self.fresh(prepared)
        if edit:
            edit(sparse, dense)
        def run(command, **kwargs):
            self.assertTrue(kwargs['check'])
            self.assertEqual(kwargs['env']['OPENBLAS_NUM_THREADS'], '1')
            self.assertEqual(kwargs['env']['OMP_NUM_THREADS'], '1')
            self.assertIn('-B', command)
            target = Path(command[command.index('--output')+1])
            self.assertFalse(target.exists())
            if 'sparse_t64.py' in command[2]:
                result = sparse
                self.assertIn('--replay', command)
            else:
                result = dense
                self.assertEqual(command[3], 'replay')
            target.write_text(json.dumps(result))
            return subprocess.CompletedProcess(command, 0)
        return run

    def test_exact_row_sum_ignores_all_search_endpoints(self):
        prepared = self.prepared()
        sparse, dense = self.fresh(prepared)
        result = whole._assemble_fresh(prepared, sparse, dense, 384, 40)
        self.assertEqual(whole.dyadic(result['upper']), Q(34, 1 << 60))
        self.assertTrue(result['whole_code_certificate'])
        self.assertFalse(result['implementation_claim'])
        self.assertNotEqual(prepared['dense_count_sha256'], prepared['sparse_count_sha256'])

    def test_search_scope_map_counts_sources_and_claim_rejected(self):
        edits = [lambda d, s: d.update(schema='old-s19-search'),
            lambda d, s: s.update(schema='s16-t128-search'),
            lambda d, s: d['scope'].update(physical_t=128),
            lambda d, s: s['inner'].update(macro_steps=32768),
            lambda d, s: s.update(distance='11/100', threshold=230686),
            lambda d, s: s['map_record'].update(map_sha256='b'*64),
            lambda d, s: d['scope'].update(maps_sha256='c'*64),
            lambda d, s: s['count_premises'].update(extra=True),
            lambda d, s: s.update(count_sha256='d'*64),
            lambda d, s: d['scope'].update(comparison_caps_sha256='e'*64),
            lambda d, s: s.update(support_min=s['support_min']+1),
            lambda d, s: d['source'].update(sha256='f'*64),
            lambda d, s: d['scope'].update(minimum_groups=34),
            lambda d, s: d['scope'].update(birth_density='classes')]
        original = deepcopy((self.dense, self.sparse))
        for edit in edits:
            self.dense, self.sparse = deepcopy(original)
            edit(self.dense, self.sparse)
            self.save()
            with self.subTest(edit=edit), self.assertRaises((ValueError, KeyError)):
                self.prepared()

    def test_partial_duplicate_and_false_complete_coverage_rejected(self):
        edits = [lambda d, s: s['results'].pop(),
            lambda d, s: s['results'].__setitem__(0, deepcopy(s['results'][1])),
            lambda d, s: s['results'][0].update(passed=False),
            lambda d, s: s['requested'].pop(),
            lambda d, s: s['results'][1]['details'].update(leaves=[]),
            lambda d, s: d['cover']['leaves'].pop('1'),
            lambda d, s: d['cover']['unresolved'].update({'1': d['cover']['leaves'].pop('1')}),
            lambda d, s: d['cover']['leaves'].update({'': dict(witness={})}),
            lambda d, s: d['cover']['leaves']['0'].update(cell=['0', '1'])]
        original = deepcopy((self.dense, self.sparse))
        for edit in edits:
            self.dense, self.sparse = deepcopy(original)
            edit(self.dense, self.sparse)
            self.save()
            with self.subTest(edit=edit), self.assertRaises((ValueError, KeyError)):
                self.prepared()
        # Even a fully covered *narrower* root must be rejected.
        self.dense, self.sparse = deepcopy(original)
        self.dense['scope']['root'] = ['1/2', '1']
        self.dense['cover']['leaves'] = {'': dict(witness={})}
        self.save()
        with self.assertRaises(ValueError):
            self.prepared()

    def test_fresh_scope_geometry_provenance_and_rows_rejected(self):
        prepared = self.prepared()
        edits = [lambda s, d: s.update(schema='old-replay'),
            lambda s, d: s.update(sources=[]), lambda s, d: d.update(input_source={}),
            lambda s, d: d.update(source={}), lambda s, d: s.update(precision=256),
            lambda s, d: s.update(target_bits=40), lambda s, d: s.update(threshold=1),
            lambda s, d: s['inner'].update(s=19),
            lambda s, d: s['map_record'].update(map_sha256='c'*64),
            lambda s, d: s.update(count_sha256='c'*64),
            lambda s, d: s['results'].pop(),
            lambda s, d: s['results'].__setitem__(0, deepcopy(s['results'][1])),
            lambda s, d: s['results'][0].update(passed=False),
            lambda s, d: d['cells'].pop(),
            lambda s, d: d['cells'][0].update(path='1'),
            lambda s, d: d['cells'][0].update(cell=['0', '1']),
            lambda s, d: d.update(complete_dense=False),
            lambda s, d: d['scope'].update(flush=True)]
        for edit in edits:
            sparse, dense = self.fresh(prepared)
            edit(sparse, dense)
            with self.subTest(edit=edit), self.assertRaises(ValueError):
                whole._assemble_fresh(prepared, sparse, dense, 384, 40)

    def test_exact_strict_target_and_aggregate_consistency(self):
        prepared = self.prepared()
        sparse, dense = self.fresh(prepared)
        sparse['aggregate_upper'] = [1, -60]
        with self.assertRaises(ValueError):
            whole._assemble_fresh(prepared, sparse, dense, 384, 40)
        sparse, dense = self.fresh(prepared)
        # Exactly 2^-40 must fail even when every individual sparse row passes.
        desired = Q(2)**-40-Q(32, 1 << 60)
        for row in dense['cells']:
            row['upper'] = whole.exact_endpoint(desired/2)
        dense['aggregate_upper'] = whole.exact_endpoint(desired)
        with self.assertRaises(ValueError):
            whole._assemble_fresh(prepared, sparse, dense, 384, 40)
        for bad in ([0, 0], [-1, -10], [True, -10], ['1', '-10'], [1], '0.1'):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                whole.dyadic(bad)
        self.assertEqual(whole.dyadic(whole.exact_endpoint(Q(3, 1 << 700))), Q(3, 1 << 700))

    def test_upward_compaction_is_json_safe_for_tiny_dense_cells(self):
        tiny = Q(2)**-100000
        for value in (tiny, Q(2)**-50+tiny, Q(2)**-50-tiny, Q(3, 1 << 200)):
            endpoint = whole.compact_endpoint(value, 416)
            self.assertGreaterEqual(whole.dyadic(endpoint), value)
            self.assertLessEqual(endpoint[0].bit_length(), 416)
            self.assertLess(len(json.dumps(endpoint)), 200)
        prepared = self.prepared()
        sparse, dense = self.fresh(prepared)
        dense['cells'][0]['upper'] = [1, -100000]
        # Its existing 2^-59 aggregate still dominates the fresh dense sum.
        result = whole._assemble_fresh(prepared, sparse, dense, 384, 40)
        exact = Q(33, 1 << 60)+tiny
        self.assertGreaterEqual(whole.dyadic(result['upper']), exact)
        self.assertLess(whole.dyadic(result['upper']), Q(2)**-40)
        self.assertLess(len(json.dumps(result)), 30000)

    def test_upward_compaction_rechecks_strict_target(self):
        prepared = self.prepared()
        sparse, dense = self.fresh(prepared)
        sparse_sum = Q(32, 1 << 60)
        dense_sum = Q(2)**-40-sparse_sum-Q(2)**-100000
        for row in dense['cells']:
            row['upper'] = whole.exact_endpoint(dense_sum/2)
        dense['aggregate_upper'] = whole.compact_endpoint(dense_sum, 416)
        self.assertLess(sparse_sum+dense_sum, Q(2)**-40)
        self.assertEqual(whole.dyadic(whole.compact_endpoint(sparse_sum+dense_sum, 416)), Q(2)**-40)
        with self.assertRaisesRegex(ValueError, 'upward-rounded'):
            whole._assemble_fresh(prepared, sparse, dense, 384, 40)

    def test_fresh_two_process_orchestration_and_no_saved_acceptance_mode(self):
        prepared = self.prepared()
        output = self.directory/'fresh'
        with patch.object(whole, 'proof_source_manifest', side_effect=self.mock_manifest), \
             patch.object(whole.subprocess, 'run', side_effect=self.child(prepared)) as runner:
            result = whole.run(self.dp, self.sp, output)
            self.assertEqual(runner.call_count, 2)
        self.assertEqual(whole.dyadic(result['upper']), Q(34, 1 << 60))
        self.assertEqual(result['source_manifest_before'], result['source_manifest_after'])
        self.assertEqual(set(result['fresh_receipts']), {'dense', 'sparse'})
        self.assertTrue((output/'whole.json').is_file())
        self.assertTrue((output/'sparse.log').is_file())
        with patch.object(whole.subprocess, 'run') as runner, self.assertRaises(ValueError):
            whole.run(self.dp, self.sp, output)
        runner.assert_not_called()

    def test_parallel_dense_dispatch_and_worker_metadata(self):
        prepared = self.prepared()
        sparse, dense = self.fresh(prepared)
        dense.update(workers=3, worker_start_method='spawn', worker_metadata=[dict(pid=123,
            precision=384, source=prepared['proposal_source'],
            scope_sha256=whole.dense_t64.prior.fingerprint(prepared['scope']),
            initialization='fresh dense_t64.fresh_model in spawned process')])
        commands = []
        def child(command, **kwargs):
            commands.append(command)
            output = Path(command[command.index('--output')+1])
            if len(commands) == 1:
                self.assertTrue(command[2].endswith('sparse_t64.py'))
                record = sparse
            else:
                self.assertTrue(command[2].endswith('replay_t64_parallel.py'))
                self.assertEqual(command[command.index('--workers')+1], '3')
                self.assertNotIn('--variance-bins', command)
                self.assertNotIn('--birth-density', command)
                self.assertNotIn('--distance', command)
                record = dense
            output.write_text(json.dumps(record))
        with patch.object(whole, 'proof_source_manifest', side_effect=self.mock_manifest), \
             patch.object(whole.subprocess, 'run', side_effect=child):
            result = whole.run(self.dp, self.sp, self.directory/'parallel', dense_workers=3)
        self.assertEqual(len(commands), 2)
        self.assertEqual(result['dense_replay_workers'], 3)
        self.assertEqual(result['dense_worker_metadata'], dense['worker_metadata'])
        dense['worker_start_method'] = 'fork'
        commands.clear()
        with patch.object(whole, 'proof_source_manifest', side_effect=self.mock_manifest), \
             patch.object(whole.subprocess, 'run', side_effect=child), self.assertRaises(ValueError):
            whole.run(self.dp, self.sp, self.directory/'badparallel', dense_workers=3)
        for workers in (0, 5, True, 2.0):
            with patch.object(whole.subprocess, 'run') as runner, self.assertRaises(ValueError):
                whole.run(self.dp, self.sp, self.directory/'invalid', dense_workers=workers)
            runner.assert_not_called()

    def test_failed_child_no_output_and_midrun_mutation_fail_closed(self):
        prepared = self.prepared()
        for name, child in [('failed', subprocess.CalledProcessError(1, ['test'])),
                            ('missing', lambda *a, **k: subprocess.CompletedProcess(a, 0))]:
            with patch.object(whole, 'proof_source_manifest', side_effect=self.mock_manifest), \
                 patch.object(whole.subprocess, 'run', side_effect=child), \
                 self.assertRaises((ValueError, subprocess.CalledProcessError)):
                whole.run(self.dp, self.sp, self.directory/name)
            self.assertFalse((self.directory/name/'whole.json').exists())
        ordinary = self.child(prepared)
        def changed(command, **kwargs):
            result = ordinary(command, **kwargs)
            self.map_path.write_text('changed map bytes')
            return result
        with patch.object(whole, 'proof_source_manifest', side_effect=self.mock_manifest), \
             patch.object(whole.subprocess, 'run', side_effect=changed) as runner, self.assertRaises(ValueError):
            whole.run(self.dp, self.sp, self.directory/'changed')
        self.assertEqual(runner.call_count, 1)
        self.assertFalse((self.directory/'changed/whole.json').exists())

    def test_manifest_detects_source_additions_and_input_changes(self):
        repo = self.directory/'repository'
        proof = repo/'proof'
        proof.mkdir(parents=True)
        (proof/'module.py').write_text('x=1\n')
        headers = repo/'spin/src/kernels/generated'
        headers.mkdir(parents=True)
        for name in ('BchCircuit.h', 'SelectedMaps.h'):
            (headers/name).write_text('// test\n')
        with patch.object(whole, 'PERMUTATION', proof), patch.object(whole, 'REPO', repo), \
             patch.object(whole.dense_t64.kernel_t64, 'SELECTED_MAP', self.map_path):
            initial = whole.proof_source_manifest([whole.source_metadata(self.dp)])
            paths = {row['path'] for row in initial['files']}
            self.assertIn(str(self.map_path.resolve()), paths)
            self.assertIn(str(self.dp.resolve()), paths)
            (proof/'new_module.py').write_text('x=2\n')
            self.assertNotEqual(initial, whole.proof_source_manifest([whole.source_metadata(self.dp)]))
            after = whole.proof_source_manifest([whole.source_metadata(self.dp)])
            self.dp.write_text('{}')
            self.assertNotEqual(after, whole.proof_source_manifest([whole.source_metadata(self.dp)]))

    def test_duplicate_json_keys_are_not_silently_discarded(self):
        self.dp.write_text('{"cover":{},"cover":{}}')
        with self.assertRaises(ValueError):
            whole.read_source(self.dp)


if __name__ == '__main__':
    unittest.main()
