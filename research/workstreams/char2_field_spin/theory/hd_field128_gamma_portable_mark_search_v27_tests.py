"""Pure validation and mocked orchestration only; never generate coefficients."""
from contextlib import contextmanager
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import hd_field128_gamma_portable_mark_search_v27 as driver


def seed():
    return dict(active_groups=152, construction='mock-construction', state_dimension=30,
        target_distance='2/5', ancestral_bound=-999,
        witness=dict(z='1/8', alpha='1/128', markers=[['2', '3']], grid=128,
                     state_dimension=30, ancestral_bound=-1000))


def search_result(proposal, weight, *, has_best=True):
    witness = deepcopy(proposal['witness'])
    best = dict(active_groups=weight, construction=proposal['construction'], state_dimension=30,
        target_distance='2/5', complete_selected_weight_screen=True, all_profiles_retained=True,
        profile_pairs_evaluated=(weight//2+1)**2, directed=False, selected_weight_certificate=False,
        whole_code_certificate=False, old_numerical_bounds_reused=False,
        common_potential_used=False, log2_upper_screen=-88., evaluation_index=0,
        witness=witness, regional_logs_sha256='a'*64)
    return dict(schema='gamma-midpoint-fixed-marker-mark-search-FLOAT-v1', active_groups=weight,
        construction=proposal['construction'], state_dimension=30, target_distance='2/5',
        directed=False, selected_weight_certificate=False, whole_code_certificate=False,
        old_numerical_bounds_reused=False, common_potential_used=False,
        max_evaluations=32, max_seconds=75., initial_witness=deepcopy(witness),
        fixed_positive_marker_literals=deepcopy(witness['markers']),
        initial_z_preserved_exactly=witness['z'], initial_alpha_preserved_exactly=witness['alpha'],
        attempted_banks=1, completed_trials=int(has_best), best=best if has_best else None,
        history=[dict(evaluation_index=0, z=witness['z'], alpha=witness['alpha'],
                      log2_upper_screen=-88.)] if has_best else [], failures=[])


def replay_result(proposal, weight, *, complete=True, score=-88.):
    return dict(construction=proposal['construction'], state_dimension=30, target_distance='2/5',
        directed=False, whole_code_certificate=False, old_numerical_bounds_reused=False,
        no_weight_interpolation=True, witness=deepcopy(proposal['witness']),
        requested_weights=[weight], status='complete' if complete else 'timeout',
        pending_weights=[] if complete else [weight],
        completed_weights=[dict(active_groups=weight, profile_pairs_evaluated=(weight//2+1)**2,
            complete_selected_weight_screen=True, all_profiles_retained=True, directed=False,
            whole_code_certificate=False, log2_upper_screen=score,
            regional_logs_sha256='a'*64)] if complete else [])


def normalize(value, root):
    if isinstance(value, list):
        return [normalize(v, root) for v in value]
    if not isinstance(value, dict):
        return value
    result = {k: normalize(v, root) for k, v in value.items()}
    if 'path' in result and 'sha256' in result and Path(result['path']).is_absolute():
        result['path'] = Path(result['path']).relative_to(root).as_posix()
    return result


class Tests(unittest.TestCase):
    def test_target_and_exact_retargeted_projection(self):
        for value in (True, 1, 4096, 172., '172'):
            with self.assertRaises(ValueError):
                driver.validate_target(value)
        parent = seed()
        projected = driver.calibration(parent, dict(path='seed.json', sha256='a'*64), 172)
        self.assertEqual(projected['active_groups'], 172)
        self.assertEqual(projected['source_calibration_active_groups'], 152)
        self.assertEqual(projected['source_sha256'], {})
        self.assertNotIn('ancestral_bound', projected)
        self.assertNotIn('ancestral_bound', projected['witness'])
        self.assertFalse(projected['old_numerical_bounds_reused'])
        projected['witness']['markers'][0][0] = '7'
        self.assertEqual(parent['witness']['markers'], [['2', '3']])
        for value in ('1/3', '1', '1/4294967297'):
            bad = seed(); bad['witness']['z'] = value
            with self.assertRaises(ValueError):
                driver.calibration(bad, {}, 172)

    def test_best_validation_and_missing_complete_trial(self):
        proposal = driver.calibration(seed(), {}, 172)
        valid = search_result(proposal, 172)
        self.assertEqual(driver.validate_search(valid, proposal, 172)['evaluation_index'], 0)
        self.assertIsNone(driver.validate_search(search_result(proposal, 172, has_best=False), proposal, 172))
        for mutate in (
                lambda r: r.update(active_groups=152),
                lambda r: r.update(attempted_banks=33),
                lambda r: r.update(completed_trials=2),
                lambda r: r['history'][0].update(log2_upper_screen=float('nan')),
                lambda r: r['best'].update(profile_pairs_evaluated=1),
                lambda r: r['best'].update(all_profiles_retained=False),
                lambda r: r['best'].update(log2_upper_screen=-90),
                lambda r: r['best']['witness'].update(markers=[['9', '3']]),
                lambda r: r['best'].update(evaluation_index=4),
                lambda r: r.update(directed=True)):
            bad = deepcopy(valid); mutate(bad)
            with self.assertRaises(ArithmeticError):
                driver.validate_search(bad, proposal, 172)

    def test_best_projection_carries_only_literals_and_lineage(self):
        proposal = driver.calibration(seed(), {}, 172)
        best = search_result(proposal, 172)['best']
        projected = driver.best_proposal(proposal, best, dict(path='run/search.json', sha256='b'*64))
        self.assertEqual(projected['source_sha256'], {})
        self.assertNotIn('log2_upper_screen', projected)
        self.assertNotIn('regional_logs_sha256', projected)
        self.assertEqual(projected['source']['path'], 'run/search.json')
        self.assertEqual(projected['witness']['markers'], proposal['witness']['markers'])

    def _run(self, directory, *, has_best=True, complete=True, score=-88., mutate=False,
             maximum_bytes=None):
        root, parent, inside = Path(directory), seed(), []
        raw = json.dumps(parent).encode(); (root/'seed.json').write_bytes(raw)
        source = dict(path='seed.json', sha256=sha256(raw).hexdigest())
        @contextmanager
        def limits():
            inside.append(True)
            try:
                yield
            finally:
                inside.pop()
        def load(name, expected):
            self.assertEqual((name, expected), (source['path'], source['sha256']))
            return parent, None, source
        calls = []
        def search(path, weight, **options):
            self.assertTrue(inside)
            self.assertEqual(options, dict(max_evaluations=32, max_seconds=75.))
            self.assertEqual(weight, 172)
            proposal = json.loads(path.read_bytes())
            self.assertEqual(proposal['active_groups'], 172)
            self.assertEqual(proposal['source_calibration_active_groups'], 152)
            result = search_result(proposal, weight, has_best=has_best)
            result['source'] = dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest())
            if mutate:
                (root/'seed.json').write_bytes(b'mutated')
            calls.append('search')
            return result
        def replay(path, weights, **options):
            self.assertTrue(inside)
            self.assertEqual(weights, (172,))
            self.assertEqual(options, dict(max_seconds=15.))
            calls.append('replay')
            result = replay_result(json.loads(path.read_bytes()), 172, complete=complete, score=score)
            result['source'] = dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest())
            return result
        loaded = (SimpleNamespace(load_seed=load, relative_name=lambda value: value),
            SimpleNamespace(normalize_new_metadata=normalize), SimpleNamespace(run=replay))
        def helpers():
            self.assertTrue(inside)
            return loaded, SimpleNamespace(run_search=search)
        with patch.object(driver, 'HERE', root), patch.object(driver.portable, 'HERE', root), \
                patch.object(driver.portable, 'remote_limits', limits), \
                patch.object(driver, 'helpers', helpers), patch.object(driver, 'source_pins', return_value={}), \
                patch.object(driver.portable, 'reserve_guard'), \
                patch.object(driver.portable, 'MAX_OUTPUT_BYTES', maximum_bytes or (1<<20)):
            result = driver.run(source['path'], source['sha256'], 172, 'run')
        return result, calls

    def test_mock_complete_replay_portability_and_no_trial_files(self):
        with tempfile.TemporaryDirectory() as directory:
            result, calls = self._run(directory)
            self.assertEqual(calls, ['search', 'replay'])
            self.assertEqual(result['status'], 'complete')
            root = Path(directory)/'run'
            self.assertEqual({p.name for p in root.iterdir()},
                {'driver.json','calibration.json','search.json','best_proposal.json','replay.json','summary.json'})
            saved = json.loads((root/'search.json').read_bytes())
            self.assertEqual(saved['source']['path'], 'run/calibration.json')
            summary = json.loads((root/'summary.json').read_bytes())
            self.assertEqual(summary['profile_pairs_evaluated'], 87**2)
            self.assertTrue(summary['exact_score_and_regional_hash_match_search'])
            self.assertFalse(summary['directed'])
            self.assertFalse(summary['whole_code_certificate'])

    def test_mock_incomplete_and_mismatch_are_not_complete(self):
        for options, status in ((dict(has_best=False),'search-incomplete'),
                (dict(complete=False),'replay-incomplete'), (dict(score=-87.),'replay-mismatch')):
            with tempfile.TemporaryDirectory() as directory:
                result, calls = self._run(directory, **options)
                self.assertEqual(result['status'], status)
                self.assertIsNone(result['log2_upper_screen'])
                self.assertEqual(calls, ['search'] if options.get('has_best') is False else ['search','replay'])

    def test_mutation_and_cumulative_budget_fail_without_complete_summary(self):
        for options, error in ((dict(mutate=True),RuntimeError), (dict(maximum_bytes=128),OSError)):
            with tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(error):
                    self._run(directory, **options)
                self.assertFalse((Path(directory)/'run/summary.json').exists())

    def test_reused_resource_limits_and_frozen_helper_identity(self):
        self.assertEqual((driver.portable.WALL_SECONDS,driver.portable.MEMORY_BYTES,
                          driver.portable.RESERVE,driver.portable.MAX_OUTPUT_BYTES),
                         (120,4<<30,20<<30,1<<20))
        self.assertEqual(sha256(Path(driver.portable.__file__).read_bytes()).hexdigest(), driver.PORTABLE_SHA256)
        self.assertEqual(sha256((driver.HERE/(driver.SEARCH_NAME+'.py')).read_bytes()).hexdigest(),
                         driver.SEARCH_SHA256)


if __name__ == '__main__':
    unittest.main()
