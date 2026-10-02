import copy
from contextlib import redirect_stderr
from fractions import Fraction as Q
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import dense_cover as cover


class DenseScope(unittest.TestCase):
    def record(self):
        return dict(schema=cover.SCHEMA, ensemble=cover.ENSEMBLE,
            K=1 << 20, N=1 << 21, updates=2, block_width=8,
            minimum_groups=33, maximum_groups=2048, distance='19/200', threshold=199229,
            comparison='expected-cdf-shell-majorant', last_lp=104, refined=True,
            base_tilt='3/16', variance_bins=16, regional_count=True,
            cover=dict(leaves={}, unresolved={'': {}}))

    def test_scope_acceptance(self):
        record = self.record()
        self.assertEqual(cover.validate_record(record), record['cover'])

    def test_other_constructions_rejected(self):
        for key, value in [('ensemble', 'shared4-gf16-r2'), ('K', 1 << 18),
                           ('block_width', 4), ('updates', 4), ('threshold', 199230),
                           ('minimum_groups', True), ('maximum_groups', 1024), ('refined', False),
                           ('comparison', 'actual-shells'), ('base_tilt', '1/4')]:
            record = copy.deepcopy(self.record())
            record[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                cover.validate_record(record)

    def test_cutoff_is_exact(self):
        self.assertEqual(cover.cutoff('.095'), 199229)
        self.assertEqual(cover.cutoff(Q(19, 200)), 199229)
        for value in ('0', '-.1', '.5', '1'):
            with self.assertRaises(ValueError):
                cover.cutoff(value)

    def test_dyadic_endpoints(self):
        self.assertEqual(cover.dyadic([3, -7]), Q(3, 128))
        self.assertEqual(cover.dyadic([0, -9]), 0)
        for value in ([-1, 0], [True, 0], [1, .5], [1], '1'):
            with self.assertRaises(ValueError):
                cover.dyadic(value)

    def test_cap_fingerprints_are_rational(self):
        self.assertEqual(cover.fingerprint([0, Q(1, 2), 4]),
                         cover.fingerprint([Q(0), Q(2, 4), Q(4)]))
        self.assertNotEqual(cover.fingerprint([0, Q(1, 2), 4]),
                            cover.fingerprint([0, Q(1, 3), 4]))

    def test_cli_rejects_search_flags_before_replay(self):
        with tempfile.TemporaryDirectory() as folder:
            output = str(Path(folder)/'unused.json')
            for flags in (['--regional-presplit'], ['--point-seeds', 'missing.json'],
                          ['--region-cache'], ['--seed-width', '1/256']):
                argv = ['dense_cover.py', '--replay', 'missing.json', '--output', output, *flags]
                with self.subTest(flags=flags), patch('sys.argv', argv), redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as caught:
                        cover.main()
                    self.assertEqual(caught.exception.code, 2)
                self.assertFalse(Path(output).exists())

    def test_cli_seed_width_constraints(self):
        with tempfile.TemporaryDirectory() as folder:
            output = str(Path(folder)/'unused.json')
            for flags in (['--seed-width', '1/256'],
                          ['--seed-width', '1/2048', '--point-seeds', 'missing.json'],
                          ['--seed-width', '2', '--point-seeds', 'missing.json'],
                          ['--seed-width', '1/256', '--point-seeds', 'missing.json', '--regional-presplit']):
                argv = ['dense_cover.py', '--output', output, *flags]
                with self.subTest(flags=flags), patch('sys.argv', argv), redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as caught:
                        cover.main()
                    self.assertEqual(caught.exception.code, 2)
            argv = ['dense_cover.py', '--output', output, '--seed-width', '1/256',
                    '--point-seeds', 'missing.json', '--region-cache', '--max-cells', '0']
            with patch('sys.argv', argv), patch.object(cover, 'authenticated_bch_cdf',
                    side_effect=RuntimeError('validated CLI; no numeric work')) as authenticate:
                with self.assertRaisesRegex(RuntimeError, 'validated CLI'):
                    cover.main()
                authenticate.assert_called_once_with()
            self.assertFalse(Path(output).exists())

    def test_parallel_workers_are_bounded_and_replay_only(self):
        with tempfile.TemporaryDirectory() as folder:
            output = str(Path(folder)/'unused.json')
            for flags in (['--replay-workers', '3'],
                          ['--replay-workers', '0', '--replay', 'missing.json'],
                          ['--replay-workers', '5', '--replay', 'missing.json']):
                argv = ['dense_cover.py', '--output', output, *flags]
                with self.subTest(flags=flags), patch('sys.argv', argv), redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as caught:
                        cover.main()
                    self.assertEqual(caught.exception.code, 2)
            self.assertFalse(Path(output).exists())


if __name__ == '__main__':
    unittest.main()
