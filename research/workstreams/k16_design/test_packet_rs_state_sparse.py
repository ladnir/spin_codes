"""Small-kernel checks for state-aware RS16 sparse receipts; no s20 census."""
from contextlib import redirect_stdout
from fractions import Fraction as Q
import hashlib
from io import StringIO
import json
from math import comb
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx
import packet_rs_state_sparse as screen


class BasisImages:
    """Minimal test double for indexed basis images, never an accepted census."""
    def __getitem__(self, index):
        return index


class StateSparseTests(unittest.TestCase):
    def setUp(self):
        self.old_precision = ctx.prec
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.source = Path(self.folder.name) / 'base16.json'
        self.source.write_text('{"base_state_bits":16}\n')
        self.digest = hashlib.sha256(self.source.read_bytes()).hexdigest()
        self.physical = dict(bits=20, distribution='uniform_gl', map_images=BasisImages(),
                             columns=[1 << (j % 20) for j in range(64)])
        self.data = dict(bits=20, map_sha256='fresh-derived-map', physical_step_bits=64,
                         macro_windows=32, macro_step_bits=128)
        self.map_record = dict(t=64, s=20, distribution='uniform_gl',
            map_sha256='fresh-derived-map', expansion_rows_hex=[hex(1 << j) for j in range(20)],
            feedback_columns=self.physical['columns'][:],
            source=dict(path=str(self.source), sha256=self.digest),
            source_role='base16 only; actual derived map is in the record')
        self.geometry = screen.q1.Geometry(32, 64, 128)

    def tearDown(self):
        ctx.prec = self.old_precision

    def test_actual_dimension_and_base_source_are_bound(self):
        with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical):
            self.assertEqual(screen.validated(self.data, self.map_record), 20)
        sources = screen.source_snapshot(self.map_record)
        self.assertEqual(sources[str(self.source.resolve())], self.digest)
        self.assertIn(str(Path(screen.__file__).resolve()), sources)

    def test_stale_s16_record_is_rejected_before_census(self):
        record = dict(self.map_record, s=16)
        with patch.object(screen.q1.kernel_t64, 'authenticate') as authenticate, \
                self.assertRaises(ValueError):
            screen.validated(self.data, record)
        authenticate.assert_not_called()

    def test_derived_rows_columns_hash_and_distribution_must_match(self):
        bad_records = [dict(self.map_record, expansion_rows_hex=['0x0'] * 20),
                       dict(self.map_record, feedback_columns=[0] * 64),
                       dict(self.map_record, map_sha256='stale'),
                       dict(self.map_record, full_state_census=False),
                       dict(self.map_record, state_count=1 << 16),
                       dict(self.map_record, distribution='transvections'),
                       dict(self.map_record, t=128)]
        with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical):
            for record in bad_records:
                with self.subTest(record=record), self.assertRaises(ValueError):
                    screen.validated(self.data, record)

    def test_changed_base_or_constructor_file_is_rejected(self):
        bad = dict(self.map_record, source=dict(path=str(self.source), sha256='0' * 64))
        constructor = dict(self.map_record, constructor_source=dict(path=str(self.source), sha256='0' * 64))
        with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical):
            for record in (bad, constructor):
                with self.assertRaises(ValueError):
                    screen.validated(self.data, record)

    def test_identity_kernel_counts_every_message_and_uses_no_old_evaluator(self):
        local = [arb_mat([[1]]) for _ in range(33)]
        for occupancy in (1, 2):
            with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical), \
                    patch.object(screen.q1.kernel_t64, 'local_operators', return_value=local), \
                    patch.object(screen.q1, 'evaluate_q1', side_effect=AssertionError('old q1 evaluator')), \
                    patch.object(screen.q2, 'evaluate_q2', side_effect=AssertionError('old q2 evaluator')), \
                    redirect_stdout(StringIO()):
                record = screen.run_low(occupancy, data=self.data, map_record=self.map_record,
                    geometry=self.geometry, tilts=['.001'], precision=192)
            mantissa, exponent = record[f'q{occupancy}_upper']
            exact = comb(32, occupancy) * (2**128 - 1)**occupancy
            upper = Q(mantissa) * Q(2)**exponent
            self.assertGreaterEqual(upper, exact)
            self.assertLessEqual(upper - exact, Q(exact) * Q(2)**-180)
            self.assertEqual(record['state_bits'], 20)
            self.assertEqual(record['inner_group'], 'GL(20,2)')
            self.assertIn('inner GL(20,2)', record['scope'])
            self.assertEqual(record['macros_per_region'], 1)
            self.assertEqual(record['threshold'], self.geometry.N // 10)
            self.assertFalse(record['whole_code_certificate'])
            self.assertEqual(record['occupancy_covered'], [occupancy])
            self.assertEqual(record['count_kind'], 'exact_expected_shells')
            self.assertEqual(record['count_sha256'], screen.shell_hash(
                screen.expected_group_support_counts(16, 8, 4, 4)))

    def test_k16_and_k20_geometry_remain_explicit(self):
        for group_count, message_bits, macros in ((512, 65536, 16), (8192, 1048576, 256)):
            geometry = screen.q1.Geometry(group_count, 64, 128)
            self.assertEqual(geometry.K, message_bits)
            self.assertEqual(geometry.macros_per_region, macros)
            screen._options(1, geometry, ['.001'], 192, None)
        with self.assertRaises(ValueError):
            screen._options(1, screen.q1.Geometry(32, 2, 4), ['.001'], 192, None)

    def test_invalid_options_and_existing_output_fail_before_operators(self):
        with patch.object(screen.q1.kernel_t64, 'local_operators') as operators:
            for occupancy in (0, 3, True):
                with self.assertRaises(ValueError):
                    screen.run_low(occupancy, data=self.data, map_record=self.map_record, tilts=['.001'])
            for tilts in ([], ['0'], ['nan'], ['.001', '1/1000'], '.001'):
                with self.assertRaises(ValueError):
                    screen.run_low(1, data=self.data, map_record=self.map_record, tilts=tilts)
            with self.assertRaises(ValueError):
                screen.run_low(1, data=self.data, map_record=self.map_record,
                               tilts=['.001'], output=self.source)
            operators.assert_not_called()

    def test_supportwise_tilt_minimum_and_fresh_json(self):
        # Uniform attenuation gives a positive exact moment independent of support.
        # Both tilts improve the trivial probability-one bound; the smaller tilt
        # wins because the synthetic operator does not depend on the tilt.
        local = [arb_mat([[arb(1)/2]]) for _ in range(33)]
        output = Path(self.folder.name) / 'fresh-component.json'
        with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical), \
                patch.object(screen.q1.kernel_t64, 'local_operators', return_value=local), \
                redirect_stdout(StringIO()):
            record = screen.run_low(1, data=self.data, map_record=self.map_record,
                geometry=self.geometry, tilts=['.002', '.001'], output=output)
        self.assertEqual(set(record['support_choices']), {'1/1000'})
        self.assertEqual(json.loads(output.read_text()), record)
        self.assertEqual(len(record['trials']), 2)

    def test_source_change_is_detected_before_output(self):
        local = [arb_mat([[1]]) for _ in range(33)]
        output = Path(self.folder.name) / 'not-created.json'
        with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical), \
                patch.object(screen.q1.kernel_t64, 'local_operators', return_value=local), \
                patch.object(screen, 'source_snapshot', side_effect=[{'a': 'b'}, {'a': 'changed'}]), \
                redirect_stdout(StringIO()), self.assertRaises(RuntimeError):
            screen.run_low(1, data=self.data, map_record=self.map_record,
                           geometry=self.geometry, tilts=['.001'], output=output)
        self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
