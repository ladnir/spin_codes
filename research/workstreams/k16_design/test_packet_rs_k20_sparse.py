"""Lightweight K20 geometry, delegation, and outward-expression checks."""
from contextlib import redirect_stdout
from dataclasses import asdict
from fractions import Fraction as Q
from io import StringIO
import json
from math import comb
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx
import packet_rs_k20_sparse as screen


class K20SparseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beta, cls.counts = screen.exact_outer()

    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 192

    def tearDown(self):
        ctx.prec = self.precision

    def test_geometry_and_outer_are_explicit(self):
        geometry = screen.GEOMETRY
        self.assertEqual((geometry.K, geometry.N), (2**20, 2**21))
        self.assertEqual((geometry.group_count, geometry.regions,
                          geometry.group_dimension, geometry.macros_per_region), (8192, 64, 128, 256))
        self.assertEqual(geometry.regions * geometry.macros_per_region * 128, geometry.N)
        self.assertEqual(screen.THRESHOLD, 209715)
        self.assertEqual(self.beta, Q(2**256, 65535**8))
        self.assertEqual(sum(self.counts), 2**128 - 1)
        self.assertEqual(self.counts[0], 0)
        self.assertNotEqual(screen.count_hash(self.counts), screen.count_hash(self.counts, cumulative=True))

    def test_options_reject_nonfinite_duplicate_nonpositive_and_bad_precision(self):
        for tilts in ([], ['0'], ['-1'], ['nan'], ['inf'], ['.01', '1/100'],
                      [float('inf')], None, '.01'):
            with self.subTest(tilts=tilts), self.assertRaises(ValueError):
                screen.checked_options(tilts, 192, None)
        for precision in (True, 127, 192.0):
            with self.assertRaises(ValueError):
                screen.checked_options(['.01'], precision, None)
        self.assertEqual(screen.checked_options(['.01'], 192, None), (('1/100',), None))

    def test_bad_occupancy_and_existing_outputs_fail_before_preparation(self):
        with patch.object(screen.q1.kernel_t64, 'prepare') as prepare:
            for first, last in ((2, 3), (3, 8193), (4, 3), (True, 3), (3, 3.0)):
                with self.assertRaises(ValueError):
                    screen.run_tail(q_min=first, q_max=last)
            with tempfile.TemporaryDirectory() as folder:
                output = Path(folder) / 'occupied.json'
                output.touch()
                for runner in (screen.run_q1, screen.run_q2, screen.run_tail):
                    with self.assertRaises(ValueError):
                        runner(output=output)
            prepare.assert_not_called()

    def _receipt(self, occupancy, map_record):
        return dict(schema=f'finite-packet-q{occupancy}-diagnostic-1',
            geometry=asdict(screen.GEOMETRY), K=2**20, N=2**21, threshold=209715,
            distance='1/10', map_record=map_record, zero_initial_state=True,
            final_flush=False, occupancy_covered=[occupancy],
            whole_code_certificate=False,
            source_sha256=screen.q1.source_snapshot(),
            count_sha256=screen.count_hash(self.counts, cumulative=occupancy == 1),
            **{f'q{occupancy}_upper': [1, -50]})

    def test_q1_and_q2_delegate_exact_shells_and_k20_geometry(self):
        for occupancy, runner, module in ((1, screen.run_q1, screen.q1),
                                          (2, screen.run_q2, screen.q2)):
            map_record = {'fresh_test_map': 1}
            result = self._receipt(occupancy, map_record)
            with patch.object(screen.q1.kernel_t64, 'prepare', return_value=({}, map_record)), \
                    patch.object(module, f'evaluate_q{occupancy}', return_value=result) as evaluate:
                receipt = runner(tilts=['.00032'])
            args, kwargs = evaluate.call_args
            self.assertEqual(args, (self.counts,))
            self.assertEqual(kwargs['tilts'], ('1/3125',))
            self.assertIsNone(kwargs['output'])
            self.assertTrue(kwargs['metadata']['independent_setups_between_groups'])
            self.assertEqual(kwargs['map_record'], map_record)
            if occupancy == 1:
                self.assertEqual(kwargs['group_count'], 8192)
                self.assertEqual(kwargs['epochs_per_region'], 256)
                self.assertEqual(kwargs['count_kind'], 'shells')
            else:
                self.assertEqual(kwargs['geometry'], screen.GEOMETRY)
            self.assertEqual(receipt['outer'], 'rs16')
            self.assertTrue(receipt['fresh_computation'])
            self.assertFalse(receipt['whole_code_certificate'])
            self.assertIn(str(Path(screen.__file__).resolve()), receipt['source_sha256'])

    def test_wrapper_rejects_wrong_geometry_or_nonpositive_endpoint(self):
        for key, bad in (('threshold', 13107), ('q1_upper', [0, 0]),
                         ('whole_code_certificate', True)):
            result = self._receipt(1, {'map': 1})
            result[key] = bad
            with patch.object(screen.q1.kernel_t64, 'prepare', return_value=({}, {'map': 1})), \
                    patch.object(screen.q1, 'evaluate_q1', return_value=result), \
                    self.assertRaises((ValueError, ArithmeticError)):
                screen.run_q1(tilts=['.001'])

    def test_uniform_expression_contains_choose_beta_and_all_terminal_mass(self):
        # Each regional operator has e_zero*R=[1/4,1/4]; subsequently total
        # mass halves in either state. Thus e_zero*R^64*1=2^-64, not 2^-65.
        regional = [arb_mat([[arb(1)/4, arb(1)/4], [arb(1)/4, arb(1)/4]]) for _ in range(4)]
        result = screen.occupancy_upper(regional, occupancy=3, beta=Q(7, 3), tilt='1/1000')
        expected = (arb(comb(8192, 3)) * screen.q1.kernel_t64.aq(Q(7, 3))**3 *
                    (arb(209715) / 1000).exp() * arb(2)**-64)
        self.assertTrue(result >= expected.lower())
        self.assertLess(float(result / expected), 1.000000001)

    def test_endpoint_rounds_up_and_rejects_invalid_values(self):
        endpoint = screen.positive_endpoint(arb(1) / 3)
        self.assertGreaterEqual(screen.endpoint_value(endpoint), Q(1, 3))
        for value in (arb(0), arb(-1), arb('nan'), arb('inf')):
            with self.assertRaises(ArithmeticError):
                screen.positive_endpoint(value)
        for pair in ([0, 0], [-1, 0], [True, 0], [1, 1000001], (1, 0)):
            with self.assertRaises(ValueError):
                screen.endpoint_value(pair)

    def test_tail_uses_complete_placement_and_best_tilt_per_occupancy(self):
        map_record = {'synthetic_map': 1}
        local = [arb_mat([[1]])] * 33
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / 'fresh.json'
            with patch.object(screen.q1.kernel_t64, 'prepare', return_value=({}, map_record)), \
                    patch.object(screen.q1.kernel_t64, 'local_operators', return_value=local), \
                    patch.object(screen.q1, 'placement', return_value=local) as placement, \
                    redirect_stdout(StringIO()):
                record = screen.run_tail(q_min=3, q_max=4, tilts=['.002', '.001'], output=output)
            self.assertEqual(record['schema'], screen.TAIL_SCHEMA)
            self.assertEqual(record['occupancy_covered'], [3, 4])
            self.assertEqual(set(record['occupancy_uppers']), {'3', '4'})
            self.assertEqual(record['occupancy_choices'], {'3': '1/1000', '4': '1/1000'})
            self.assertFalse(record['whole_code_certificate'])
            for call in placement.call_args_list:
                self.assertEqual(call.kwargs['epochs'], 256)
                self.assertEqual(call.kwargs['windows'], 32)
                self.assertEqual(call.kwargs['maximum_groups'], 4)
            saved = json.loads(output.read_text())
            self.assertEqual(saved, record)
            summed = sum(map(screen.endpoint_value, record['occupancy_uppers'].values()))
            self.assertGreaterEqual(screen.endpoint_value(record['union_upper']), summed)

    def test_source_change_rejects_component_and_creates_no_output(self):
        record = self._receipt(1, {'map': 1})
        record['source_sha256'] = {'source': 'a'}
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / 'never.json'
            with patch.object(screen.q1.kernel_t64, 'prepare', return_value=({}, {'map': 1})), \
                    patch.object(screen.q1, 'evaluate_q1', return_value=record), \
                    patch.object(screen.q1, 'source_snapshot', side_effect=[{'source': 'a'}, {'source': 'b'}]), \
                    self.assertRaises(RuntimeError):
                screen.run_q1(tilts=['.001'], output=output)
            self.assertFalse(output.exists())

    def test_shared_preparation_is_reused_and_partial_inputs_rejected(self):
        map_record = {'fresh_shared_map': 1}
        data = {'fresh_shared_data': 1}
        with patch.object(screen.q1.kernel_t64, 'prepare') as prepare, \
                patch.object(screen.q1, 'evaluate_q1', return_value=self._receipt(1, map_record)) as evaluate:
            screen.run_q1(tilts=['.001'], data=data, map_record=map_record)
            self.assertIs(evaluate.call_args.kwargs['data'], data)
            prepare.assert_not_called()
        for runner in (screen.run_q1, screen.run_q2, screen.run_tail):
            with self.assertRaises(ValueError):
                runner(data=data)
            with self.assertRaises(ValueError):
                runner(map_record=map_record)


if __name__ == '__main__':
    unittest.main()
