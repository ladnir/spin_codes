"""Memory-only metadata/projection/prefactor/guard tests; no coefficient jobs."""
from fractions import Fraction as F
import json
from pathlib import Path
import unittest
from unittest.mock import Mock, patch

import hd_field128_gamma_state_ladder_v33 as pilot


def source():
    return dict(construction=pilot.FIXED30, state_dimension=30, active_groups=4096,
        target_distance='2/5', log2_upper_screen=-999, source_sha256={'old.py': 'unused'},
        witness=dict(z='1/4', alpha='1/128', markers=[['2', '3'], ['5/2', '7']], grid=128,
                     forbidden_weight=838860, target_distance='2/5'))


class Checks(unittest.TestCase):
    def test_state_semantics(self):
        for s in (30, 31):
            m = pilot.state_metadata(s)
            self.assertFalse(m['conditional_on_global_principal_nonsingularity'])
            self.assertEqual(m['setup_failure_probability_upper_rational'], '0')
        m = pilot.state_metadata(32)
        self.assertTrue(m['conditional_on_global_principal_nonsingularity'])
        self.assertEqual(m['emission_minimum'], 33)
        self.assertEqual(F(m['setup_failure_probability_upper_rational']), F(16*((1 << 62)-1), (1 << 128)-(1 << 64)))
        self.assertLess(F(m['setup_failure_probability_upper_rational']), F(1, 1 << 62))
        self.assertFalse(m['setup_failure_included_in_selected_FLOAT_score'])
        self.assertFalse(m['whole_code_certificate'])
        for s in (True, 32., 33):
            with self.assertRaises(ValueError):
                pilot.state_metadata(s)

    def test_exact_literal_projection(self):
        original = source()
        before = json.dumps(original, sort_keys=True)
        projected = pilot.project_seed(original, 5120, F(9, 20))
        self.assertEqual(projected['active_groups'], 5120)
        self.assertEqual(projected['target_distance'], '9/20')
        self.assertEqual(projected['witness']['markers'], original['witness']['markers'])
        self.assertNotIn('forbidden_weight', projected['witness'])
        self.assertNotIn('source_sha256', projected)
        self.assertNotIn('log2_upper_screen', projected)
        self.assertEqual(json.dumps(original, sort_keys=True), before)
        for key, value in (('z', '1/3'), ('alpha', '0'), ('markers', [['0', '1']])):
            bad = source(); bad['witness'][key] = value
            with self.assertRaises(ValueError):
                pilot.project_seed(bad, 4096, F(9, 20))

    def test_prefactor_distance_only(self):
        # Stub only the exact projective log correction, not numerical kernels.
        import sys
        fake = Mock(log1mexp2=lambda bits: -0.125)
        with patch.dict(sys.modules, {'hd_field128_lift_alias_screen': fake}):
            a = pilot.prefactor(.25, .125, F(2, 5), 913, 128.)
            b = pilot.prefactor(.25, .125, F(9, 20), 913, 128.)
        self.assertEqual(int(F(9, 20)*pilot.N), 943718)
        self.assertEqual(b-a, F(104858, 4))
        self.assertEqual(a, .125*(128.*913-128+.125+2*838860))

    def test_scope_and_limits(self):
        pilot.validate_options(4096, (30, 31, 32), F(9, 20), 'search', 32, 24, 150.)
        for index, replacement in ((0, 8192), (1, (30, 30)), (2, F(1, 2)), (4, 29), (5, 25), (6, 151.)):
            values = [4096, (30, 31, 32), F(9, 20), 'search', 32, 24, 150.]
            values[index] = replacement
            with self.assertRaises(ValueError):
                pilot.validate_options(*values)
        for name in ('../out', '/tmp/out', 'C:/out', 'a\\b', '.'):
            with self.assertRaises(ValueError):
                pilot.inside(name)
        with patch.object(pilot.shutil, 'disk_usage', return_value=Mock(free=pilot.RESERVE+pilot.OUTPUT_LIMIT-1)):
            with patch.object(Path, 'mkdir') as mkdir:
                with self.assertRaises(RuntimeError):
                    pilot.Writer(Path('unused'))
                mkdir.assert_not_called()
        writer = object.__new__(pilot.Writer)
        writer.directory, writer.written, writer.pins = Path('unused'), pilot.OUTPUT_LIMIT, {}
        with self.assertRaises(RuntimeError):
            writer.write('x.json', {})

    def test_run_installs_hard_guard(self):
        guard = Mock()
        guard.__enter__ = Mock(return_value=None); guard.__exit__ = Mock(return_value=False)
        with patch.object(pilot, 'remote_limits', return_value=guard), patch.object(pilot, '_run', return_value={'ok': True}) as run:
            self.assertEqual(pilot.run('seed', 'pin', 4096, 'out'), {'ok': True})
            guard.__enter__.assert_called_once()
            run.assert_called_once_with('seed', 'pin', 4096, 'out')
            guard.__exit__.assert_called_once()


if __name__ == '__main__':
    unittest.main()
