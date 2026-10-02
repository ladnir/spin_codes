import copy
from fractions import Fraction as Q
import io
import json
from pathlib import Path
import sys
import subprocess
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from contextlib import redirect_stderr

from flint import arb, ctx
import hill_cover as hill


def scope(updates=2, ekr=False):
    return dict(schema=hill.CONTEXT_SCHEMA, ensemble=f'canonical-gl32-width8-shared4-r{updates}',
        K=1 << 20, N=1 << 21, updates=updates, block_width=8, distance='1/10', threshold=209715,
        minimum_groups=33, maximum_groups=2048, root=['0', '1'],
        comparison='direct-expected-shell-majorant', expected_cdf_sha256='a'*64,
        comparison_caps_sha256='b'*64, last_lp=104, refined=True, base_tilt='3/16',
        variance_bins=16, regional_count=True,
        outer_premises=dict(schema=hill.EKR_SCHEMA if ekr else hill.INCIDENCE_SCHEMA,
            canonical_cdf=[0, 1]), mixture=[dict(mass='1', activity='1/2')])


def record(updates=2, complete=False):
    cells = {p: dict(cell=list(map(str, hill.geometry.path_cell((0, 1), p)))) for p in ('0', '1')}
    leaves = {'0': dict(cells['0'], witness={'x': 1}, upper=[1, -10000])}
    pending = {'1': cells['1']}
    if complete:
        leaves['1'] = dict(cells['1'], witness={'x': 2}, upper=[1, -10000])
        pending = {}
    return dict(schema=hill.SCHEMA, scope=scope(updates), precision=256, cell_target_bits=52,
        whole_code_certificate=False, cover=dict(leaves=leaves, unresolved=pending, visited=7))


class Model:
    root = (Q(0), Q(1))

    def __init__(self, updates=2, exponent=-60):
        self.data = {'updates': updates}
        self.exponent = exponent
        self.calls = []

    def proposal(self, cell):
        return -80., {'proposal': True}

    def outward(self, cell, witness):
        self.calls.append((cell, copy.deepcopy(witness)))
        return arb(2)**self.exponent


class HillCoverTests(unittest.TestCase):
    def setUp(self):
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)

    def test_scope_separates_updates_and_outer_refinement(self):
        for updates in (2, 3, 4):
            self.assertFalse(hill.validate_scope(scope(updates)))
            self.assertTrue(hill.validate_scope(scope(updates, True)))
        for key, value in (('updates', True), ('updates', 5), ('ensemble', 'r2'),
                ('threshold', 209714), ('variance_bins', 0), ('minimum_groups', 0),
                ('refined', 1), ('regional_count', 1)):
            changed = scope()
            changed[key] = value
            with self.assertRaises(ValueError):
                hill.validate_scope(changed)

    def test_points_never_become_intervals(self):
        original = dict(schema=hill.POINT_SCHEMA, scope=scope(4, True),
            points=[dict(cell=['.1', '.1'], checked=dict(upper=[1, -900]))])
        before = copy.deepcopy(original)
        selected, cover = hill.from_source(original)
        self.assertEqual(selected['updates'], 4)
        self.assertFalse(cover['leaves'])
        self.assertEqual(cover['unresolved'], {'': {'cell': ['0', '1']}})
        self.assertEqual(original, before)

    def test_partition_gaps_overlaps_and_incomplete_replay_rejected(self):
        for change in ('gap', 'overlap', 'cell', 'witness'):
            bad = record(complete=True)
            if change == 'gap':
                del bad['cover']['leaves']['1']
            elif change == 'overlap':
                bad['cover']['unresolved'][''] = dict(cell=['0', '1'])
            elif change == 'cell':
                bad['cover']['leaves']['0']['cell'] = ['0', '1']
            else:
                del bad['cover']['leaves']['0']['witness']
            with self.assertRaises(ValueError):
                hill.validate_record(bad, complete=True)
        with patch.object(hill, 'fresh_model') as build, self.assertRaises(ValueError):
            hill.replay(record())
        build.assert_not_called()

    def test_search_preserves_prior_witnesses_without_claiming_replay(self):
        source = record(4)
        before = copy.deepcopy(source)
        model = Model(4)
        snapshots = []
        result = hill.search(model, source['scope'], source['cover'], max_cells=1,
            checkpoint=snapshots.append)
        self.assertFalse(result['unresolved'])
        self.assertEqual(result['visited'], 8)
        self.assertEqual(result['leaves']['0'], source['cover']['leaves']['0'])
        self.assertEqual(len(model.calls), 1)
        self.assertEqual(source, before)
        for snapshot in snapshots:
            hill.validate_cover(source['scope'], snapshot)

    def split_source(self):
        source = record(4)
        source['cover']['unresolved'] = {path: dict(
            cell=list(map(str, hill.geometry.path_cell((0, 1), path))), marker=path)
            for path in ('100', '101', '11')}
        return source

    def test_selection_guards_and_contiguous_shards(self):
        cover = self.split_source()['cover']
        shards = [hill.select_paths(cover, shard_index=i, shard_count=2) for i in range(2)]
        self.assertEqual(shards, [['100'], ['101', '11']])
        self.assertEqual(sum(shards, []), sorted(cover['unresolved']))
        for paths in ([], ['100', '100'], ['1'], ['1000'], ['0'], [True]):
            with self.assertRaises(ValueError):
                hill.select_paths(cover, paths)
        for options in (dict(shard_index=0), dict(shard_count=2),
                dict(shard_index=True, shard_count=2), dict(shard_index=2, shard_count=2),
                dict(shard_index=0, shard_count=4), dict(paths=['100'], shard_count=2)):
            with self.assertRaises(ValueError):
                hill.select_paths(cover, **options)

    def test_subset_search_retains_untouched_global_partition(self):
        source = self.split_source()
        before = copy.deepcopy(source)
        snapshots = []
        model = Model(4)
        result = hill.search(model, source['scope'], source['cover'], selected_paths=['101'],
            max_cells=1, checkpoint=snapshots.append)
        self.assertEqual(result['unresolved'], {p: source['cover']['unresolved'][p] for p in ('100', '11')})
        self.assertEqual(set(result['leaves']), {'0', '101'})
        self.assertEqual(model.calls[0][0], hill.geometry.path_cell((0, 1), '101'))
        self.assertEqual(len(model.calls), 1)
        self.assertEqual(source, before)
        for snapshot in snapshots:
            hill.validate_cover(source['scope'], snapshot)
        with patch.object(hill.cell_search, 'search', return_value=dict(
                leaves={'11': dict(cell=['3/4', '1'], witness={})}, unresolved={}, visited=1)), \
                self.assertRaises(ValueError):
            hill.search(Model(4), source['scope'], source['cover'], selected_paths=['101'])

    def test_cli_shard_preserves_pending_and_rejects_replay_selection(self):
        with TemporaryDirectory() as directory:
            source, output = Path(directory)/'source.json', Path(directory)/'new.json'
            source.write_text(json.dumps(self.split_source()))
            args = ['hill_cover', str(source), '--output', str(output), '--max-cells', '1',
                '--shard-index', '0', '--shard-count', '2', '--cell-target-bits', '46']
            with patch.object(sys, 'argv', args), patch.object(hill, 'fresh_model', return_value=Model(4)):
                hill.main()
            result = json.loads(output.read_text())
            self.assertEqual(result['search_selection']['assigned_roots'], ['100'])
            self.assertEqual(result['search_selection']['source_frontier'], ['100', '101', '11'])
            self.assertEqual(set(result['cover']['unresolved']), {'101', '11'})
            self.assertFalse(result['complete_search_partition'])
            self.assertFalse(result['whole_code_certificate'])
            args[args.index(str(output))] = str(Path(directory)/'replay.json')
            with patch.object(sys, 'argv', args+['--replay']), patch.object(hill, 'fresh_model') as build, \
                    redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                hill.main()
            build.assert_not_called()

    def merge_fixture(self):
        source = self.split_source()
        metadata = dict(path=str(Path(__file__).resolve()), sha256='c'*64)
        partials = []
        for path in ('100', '101', '11'):
            partial = copy.deepcopy(source)
            partial['cell_target_bits'] = 46
            partial['source'] = metadata
            partial['search_selection'] = dict(assigned_roots=[path],
                source_frontier=sorted(source['cover']['unresolved']))
            partial['cover'] = hill.search(Model(4), source['scope'], source['cover'],
                selected_paths=[path], max_cells=1, target_bits=46)
            partials.append(partial)
        return source, metadata, partials

    def test_merge_preserves_unassigned_and_discards_all_numeric_bounds(self):
        source, metadata, partials = self.merge_fixture()
        before = copy.deepcopy((source, metadata, partials))
        partial = hill.merge_records(source, metadata, partials[:1])
        self.assertEqual(set(partial['cover']['unresolved']), {'101', '11'})
        self.assertFalse(partial['complete_search_partition'])
        with patch.object(hill, 'fresh_model') as build, self.assertRaises(ValueError):
            hill.replay(partial)
        build.assert_not_called()
        merged = hill.merge_records(source, metadata, partials)
        hill.validate_record(merged, complete=True)
        self.assertEqual(merged['cover']['visited'], source['cover']['visited']+3)
        self.assertEqual(merged['cell_target_bits'], 46)
        self.assertEqual(merged['merge']['source_and_partial_cell_targets'], [52, 46, 46, 46])
        self.assertFalse(merged['whole_code_certificate'])
        self.assertTrue(merged['final_replay_required'])
        for row in merged['cover']['leaves'].values():
            self.assertEqual(set(row), {'cell', 'witness'})
        self.assertEqual((source, metadata, partials), before)
        with patch.object(hill, 'fresh_model', return_value=Model(4, -30)):
            self.assertFalse(hill.replay(merged, 256, 40)['passed'])

    def test_merge_rejects_source_scope_assignment_and_geometry_tampering(self):
        source, metadata, partials = self.merge_fixture()
        with self.assertRaises(ValueError):
            hill.merge_records(source, metadata, [partials[0], partials[0]])
        for change in ('hash', 'root', 'premises', 'frontier', 'assignment', 'outside', 'cell', 'work'):
            bad = copy.deepcopy(partials[0])
            if change == 'hash':
                bad['source']['sha256'] = 'd'*64
            elif change == 'root':
                bad['scope']['root'] = ['0', '2']
            elif change == 'premises':
                bad['scope']['outer_premises']['canonical_cdf'] = [0, 2]
            elif change == 'frontier':
                bad['search_selection']['source_frontier'] = ['100']
            elif change == 'assignment':
                bad['search_selection']['assigned_roots'] = ['101']
            elif change == 'outside':
                bad['cover']['leaves']['0']['witness'] = {'tampered': True}
            elif change == 'cell':
                bad['cover']['leaves']['100']['cell'] = ['1/2', '1']
            else:
                bad['cover']['visited'] = 6
            with self.subTest(change=change), self.assertRaises(ValueError):
                hill.merge_records(source, metadata, [bad])

    def test_cli_merge_records_all_input_hashes_and_never_builds_model(self):
        import hashlib
        with TemporaryDirectory() as directory:
            source_path, output = Path(directory)/'source.json', Path(directory)/'merged.json'
            source, _, partials = self.merge_fixture()
            source_path.write_text(json.dumps(source))
            metadata = dict(path=str(source_path.resolve()),
                sha256=hashlib.sha256(source_path.read_bytes()).hexdigest())
            paths = []
            for index, partial in enumerate(partials):
                partial['source'] = metadata
                path = Path(directory)/f'shard{index}.json'
                path.write_text(json.dumps(partial))
                paths.append(path)
            args = ['hill_cover', str(source_path), '--output', str(output),
                '--merge', *map(str, paths)]
            with patch.object(sys, 'argv', args), patch.object(hill, 'fresh_model') as build:
                hill.main()
            build.assert_not_called()
            result = json.loads(output.read_text())
            self.assertEqual(result['source'], metadata)
            self.assertEqual(result['partial_sources'], [dict(path=str(path.resolve()),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in paths])
            self.assertEqual(result['cell_target_bits'], 46)
            self.assertTrue(result['complete_search_partition'])
            self.assertFalse(result['whole_code_certificate'])
            self.assertNotIn('fresh_replay', result)
            self.assertFalse(output.with_name(output.name+'.tmp').exists())

    def test_cli_merge_detects_input_mutation_before_output(self):
        import hashlib
        with TemporaryDirectory() as directory:
            source_path, shard_path, output = (Path(directory)/name for name in
                ('source.json', 'shard.json', 'merged.json'))
            source, _, partials = self.merge_fixture()
            source_path.write_text(json.dumps(source))
            partials[0]['source'] = dict(path=str(source_path.resolve()),
                sha256=hashlib.sha256(source_path.read_bytes()).hexdigest())
            shard_path.write_text(json.dumps(partials[0]))
            original = hill.merge_records
            def changed(*args):
                result = original(*args)
                shard_path.write_text('{}')
                return result
            args = ['hill_cover', str(source_path), '--output', str(output), '--merge', str(shard_path)]
            with patch.object(sys, 'argv', args), patch.object(hill, 'merge_records', side_effect=changed), \
                    self.assertRaises(ValueError):
                hill.main()
            self.assertFalse(output.exists())

    def test_cli_hints_are_search_only_and_checked_at_checkpoints(self):
        helper = SimpleNamespace(wrap=unittest.mock.Mock(), check_sources=unittest.mock.Mock())
        helper.wrap.return_value = (Model(4), object(), [{'path': 'mock', 'sha256': 'c'*64}])
        with TemporaryDirectory() as directory:
            source, output = Path(directory)/'source.json', Path(directory)/'new.json'
            source.write_text(json.dumps(self.split_source()))
            hint = Path(directory)/'hint.json'
            args = ['hill_cover', str(source), '--output', str(output), '--max-cells', '1',
                '--paths', '100', '--hints', str(hint)]
            with patch.object(sys, 'argv', args), patch.dict(sys.modules, {'hill_reuse': helper}), \
                    patch.object(hill, 'fresh_model', return_value=Model(4)):
                hill.main()
            self.assertEqual(helper.wrap.call_args.args[2], [hint])
            self.assertGreater(helper.check_sources.call_count, 0)
            result = json.loads(output.read_text())
            self.assertEqual(result['hint_sources'], helper.wrap.return_value[2])
            self.assertFalse(result['whole_code_certificate'])
            args[args.index(str(output))] = str(Path(directory)/'replay.json')
            with patch.object(sys, 'argv', args+['--replay']), patch.object(hill, 'fresh_model') as build, \
                    redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                hill.main()
            build.assert_not_called()

    def test_replay_rechecks_all_leaves_and_sums_new_endpoints(self):
        source = record(4, complete=True)
        model = Model(4, -60)
        with patch.object(hill, 'fresh_model', return_value=model) as build:
            result = hill.replay(source, 384, 41)
        build.assert_called_once_with(source['scope'], 384)
        self.assertEqual(len(model.calls), 2)
        self.assertEqual(hill.cell_search.dense.dyadic(result['aggregate_upper']), Q(2)**-59)
        self.assertTrue(result['complete_dense'])
        self.assertTrue(result['passed'])
        self.assertFalse(result['whole_code_certificate'])
        with patch.object(hill, 'fresh_model', return_value=Model(4, -40)):
            self.assertFalse(hill.replay(source, 384, 41)['passed'])

    def test_precision_mutation_and_nonfinite_replay_rejected(self):
        class Broken(Model):
            def outward(self, cell, witness):
                ctx.prec = 128
                return arb(2)**-80
        with patch.object(hill, 'fresh_model', return_value=Broken()), self.assertRaises(ArithmeticError):
                hill.replay(record(complete=True))

    def test_complete_search_still_validates_work_limits(self):
        source = record(complete=True)
        for options in (dict(max_cells=True), dict(target_bits=19), dict(max_depth=65)):
            with self.assertRaises(ValueError):
                hill.search(Model(), source['scope'], source['cover'], **options)
        with patch.object(Model, 'outward', return_value=arb('nan')), \
                patch.object(hill, 'fresh_model', return_value=Model()), self.assertRaises(ArithmeticError):
            hill.replay(record(complete=True))

    def test_fresh_model_authenticates_mixture_and_actual_update_count(self):
        selected = scope(4, True)
        import monotone
        with patch.object(hill, '_authenticate_counts', return_value=([0, 1], selected['outer_premises'])) as authenticate, \
                patch.object(monotone, 'transport_shells', return_value=[0, 2]), \
                patch.object(hill.cell_search.dense, 'fingerprint', side_effect=['a'*64, 'b'*64]), \
                patch.object(hill.cell_search.dense, 'exact_mixture', return_value=([(Q(1), Q(1, 2))], {})) as mixture, \
                patch.object(hill, '_construct_model', return_value=Model(4)) as build:
            self.assertEqual(hill.fresh_model(selected, 384).data['updates'], 4)
            authenticate.assert_called_once_with(True)
            mixture.assert_called_once_with([0, 2], selected['mixture'])
            self.assertEqual(build.call_args.args[0]['updates'], 4)
        with patch.object(hill, '_authenticate_counts', return_value=([0, 1], {'canonical_cdf': [0, 1]})), \
                patch.object(monotone, 'transport_shells', return_value=[0, 2]), \
                patch.object(hill, '_construct_model') as build, self.assertRaises(ValueError):
            hill.fresh_model(selected, 384)
        build.assert_not_called()

    def test_cold_mixture_import_has_legacy_paths_before_inner_factory(self):
        script = '''
import builtins,json,sys
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,sys.argv[1])
import hill_cover as hill
import monotone
scope=json.loads(sys.argv[2])
class ReachedMixture(Exception): pass
original=builtins.__import__
def importing(name,*args,**kwargs):
    result=original(name,*args,**kwargs)
    if name=='shared_mixture':
        assert Path(result.__file__).name=='shared_mixture.py'
        raise ReachedMixture
    return result
with patch.object(hill,'_authenticate_counts',return_value=([0,1],scope['outer_premises'])), \\
     patch.object(monotone,'transport_shells',return_value=[0,2]), \\
     patch.object(hill.cell_search.dense,'fingerprint',side_effect=['a'*64,'b'*64]), \\
     patch.object(hill,'_construct_model') as build, \\
     patch('builtins.__import__',side_effect=importing):
    try: hill.fresh_model(scope,256)
    except ReachedMixture: print('cold legacy mixture import reached')
    else: raise AssertionError('real mixture import was not reached')
    build.assert_not_called()
'''
        result = subprocess.run([sys.executable, '-B', '-c', script,
            str(Path(__file__).resolve().parent), json.dumps(scope(4, True))],
            capture_output=True, text=True, timeout=45)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn('cold legacy mixture import reached', result.stdout)

    def test_cli_presplit_and_new_output_guard(self):
        with TemporaryDirectory() as directory:
            source, output = Path(directory)/'source.json', Path(directory)/'new.json'
            source.write_text(json.dumps(dict(schema=hill.POINT_SCHEMA, scope=scope(3))))
            args = ['hill_cover', str(source), '--output', str(output), '--max-cells', '0',
                '--initial-width', '1/4']
            with patch.object(sys, 'argv', args), patch.object(hill, 'fresh_model', return_value=Model(3)):
                hill.main()
            result = json.loads(output.read_text())
            self.assertEqual(len(result['cover']['unresolved']), 4)
            self.assertNotIn('fresh_replay', result)
            self.assertFalse(result['whole_code_certificate'])
            with patch.object(sys, 'argv', args), patch.object(hill, 'fresh_model') as build, \
                    redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                hill.main()
            build.assert_not_called()


if __name__ == '__main__':
    unittest.main()
