"""Pure and mocked orchestration tests; no numerical banks or objectives."""
from contextlib import contextmanager
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import hd_field128_gamma_portable_neighbors_v27 as driver


def seed():
    return dict(active_groups=152, construction='mock-construction', state_dimension=30,
        target_distance='2/5', historical_bound=-999,
        witness=dict(z='1/8', alpha='1/128', markers=[['2', '3']], grid=128,
                     state_dimension=30, historical_bound=-1000))


def result(source, weights, completed=None):
    completed = len(weights) if completed is None else completed
    return dict(construction=source['construction'], state_dimension=30, target_distance='2/5',
        witness=deepcopy(source['witness']), requested_weights=list(weights),
        status='complete' if completed == len(weights) else 'timeout',
        completed_weights=[dict(active_groups=r, profile_pairs_evaluated=(r//2+1)**2,
            complete_selected_weight_screen=True, all_profiles_retained=True,
            directed=False, whole_code_certificate=False, log2_upper_screen=-1.,
            regional_alias_diagnostics=dict(selected_positive_fallback_entries=3))
            for r in weights[:completed]], pending_weights=list(weights[completed:]),
        directed=False, whole_code_certificate=False, old_numerical_bounds_reused=False,
        no_weight_interpolation=True)


def normalize(value, root):
    if isinstance(value, list):
        return [normalize(v, root) for v in value]
    if not isinstance(value, dict):
        return value
    out = {k: normalize(v, root) for k, v in value.items()}
    if 'path' in out and 'sha256' in out and Path(out['path']).is_absolute():
        out['path'] = Path(out['path']).relative_to(root).as_posix()
    return out


class Tests(unittest.TestCase):
    def test_weight_guards_and_literal_projection(self):
        for bad in ([], [1], [4096], [165, 165], [True], list(range(2, 11))):
            with self.assertRaises(ValueError):
                driver.validate_weights(bad)
        parent = seed()
        projected = driver.calibration(parent, {'path': 'parent.json', 'sha256': 'a'*64})
        self.assertEqual(projected['source_sha256'], {})
        self.assertNotIn('historical_bound', projected)
        self.assertNotIn('historical_bound', projected['witness'])
        self.assertEqual(projected['witness']['markers'], parent['witness']['markers'])
        projected['witness']['markers'][0][0] = '9'
        self.assertEqual(parent['witness']['markers'][0][0], '2')

    def test_complete_original_profiles_and_partial_fail_closed(self):
        s, weights = seed(), (165, 172, 179)
        valid = result(s, weights)
        self.assertTrue(driver.validate_result(valid, s, weights))
        partial = result(s, weights, 1)
        self.assertFalse(driver.validate_result(partial, s, weights))
        for mutate in (
                lambda x: x['completed_weights'][0].update(profile_pairs_evaluated=1),
                lambda x: x['completed_weights'][0].update(active_groups=172),
                lambda x: x['completed_weights'][0].update(log2_upper_screen=float('nan')),
                lambda x: x.update(directed=True),
                lambda x: x['witness'].update(alpha='1/64')):
            broken = deepcopy(valid)
            mutate(broken)
            with self.assertRaises(ArithmeticError):
                driver.validate_result(broken, s, weights)
        partial['status'] = 'complete'
        with self.assertRaises(ArithmeticError):
            driver.validate_result(partial, s, weights)

    def test_outputs_exclusive_budget_and_final_integrity(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(driver, 'HERE', Path(directory)), \
                patch.object(driver, 'reserve_guard'), patch.object(driver, 'MAX_OUTPUT_BYTES', 128):
            output = driver.Outputs(Path(directory), normalize)
            output.write('first.json', {'value': 1})
            with self.assertRaises(FileExistsError):
                output.write('first.json', {})
            with self.assertRaises(OSError):
                output.write('huge.json', {'value': 'x'*129})
            self.assertFalse((Path(directory)/'huge.json').exists())
            (Path(directory)/'first.json').write_bytes(b'changed')
            with self.assertRaises(RuntimeError):
                output.check()

    def test_reserve_includes_planned_bytes(self):
        with patch.object(driver.shutil, 'disk_usage', return_value=SimpleNamespace(free=driver.RESERVE+100)):
            driver.reserve_guard('.', 100)
            with self.assertRaises(OSError):
                driver.reserve_guard('.', 101)

    def test_linux_limits_precede_heavy_import_and_hard_exit(self):
        resource = SimpleNamespace(RLIMIT_AS=9, RLIM_INFINITY=-1,
            getrlimit=lambda unused: (-1, -1), setrlimit=lambda *args: None)
        signals = SimpleNamespace(SIGALRM=14, ITIMER_REAL=0, getsignal=lambda unused: 'old',
                                  signal=lambda *args: None, setitimer=lambda *args: None)
        with patch.object(driver.sys, 'platform', 'linux'), patch.dict('sys.modules', resource=resource), \
                patch.object(driver, 'signal', signals), \
                patch.object(signals, 'signal') as handler, \
                patch.object(signals, 'setitimer') as timer, \
                patch.object(driver.os, '_exit') as exit_process, \
                patch.object(resource, 'setrlimit') as limits:
            with driver.remote_limits():
                limits.assert_called_once_with(9, (4 << 30, 4 << 30))
                handler.call_args_list[0].args[1](None, None)
                exit_process.assert_called_once_with(124)
                self.assertEqual(timer.call_args_list[0].args[1], 120)
            self.assertEqual(timer.call_args_list[-1].args[1], 0)

    def _mock_run(self, directory, *, completed=3, mutate_seed=False):
        root, s, weights, active = Path(directory), seed(), (165, 172, 179), []
        raw = json.dumps(s).encode()
        (root/'seed.json').write_bytes(raw)
        digest = sha256(raw).hexdigest()
        @contextmanager
        def limits():
            active.append(True)
            try:
                yield
            finally:
                active.pop()
        def load(name, expected):
            self.assertEqual((name, expected), ('seed.json', digest))
            return s, None, dict(path=name, sha256=expected)
        def screen(path, requested, max_seconds):
            self.assertTrue(active)
            projection = json.loads(path.read_bytes())
            self.assertEqual(projection['source_sha256'], {})
            self.assertNotIn('historical_bound', projection['witness'])
            self.assertEqual(max_seconds, 90.)
            out = result(s, requested, completed)
            out['source'] = dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest())
            if mutate_seed:
                (root/'seed.json').write_bytes(b'mutated')
            return out
        loaded = (SimpleNamespace(load_seed=load, relative_name=lambda x: x),
                  SimpleNamespace(normalize_new_metadata=normalize), SimpleNamespace(run=screen))
        with patch.object(driver, 'HERE', root), patch.object(driver, 'remote_limits', limits), \
                patch.object(driver, 'helpers', side_effect=lambda: (self.assertTrue(active), loaded)[1]), \
                patch.object(driver, 'source_pins', return_value={}), patch.object(driver, 'reserve_guard'):
            return driver.run('seed.json', digest, weights, 'output')

    def test_mock_complete_and_partial_metadata(self):
        for completed in (1, 3):
            with tempfile.TemporaryDirectory() as directory:
                report = self._mock_run(directory, completed=completed)
                saved = json.loads((Path(directory)/'output/neighbors.json').read_bytes())
                self.assertEqual(saved['all_requested_weights_completed'], completed == 3)
                self.assertEqual(saved['source']['path'], 'output/calibration.json')
                self.assertEqual(saved['completed_weights'][0]['regional_alias_diagnostics']
                                 ['selected_positive_fallback_entries'], 3)
                self.assertFalse(report['whole_code_certificate'])
                self.assertEqual(len(report['complete_weights']), completed)

    def test_mutation_rejects_complete_output(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(RuntimeError):
                self._mock_run(directory, mutate_seed=True)
            self.assertFalse((Path(directory)/'output/neighbors.json').exists())


if __name__ == '__main__':
    unittest.main()
