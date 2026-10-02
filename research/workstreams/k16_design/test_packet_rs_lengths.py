"""Lightweight length, selected-map, delegation, and positive-envelope checks."""
from contextlib import redirect_stdout
from copy import deepcopy
from dataclasses import asdict
from fractions import Fraction as Q
from io import StringIO
import hashlib
import json
from math import comb
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx
import packet_rs_lengths as screen


def selected_fixture():
    """Only read basis rows; no census or prepared transition data is generated."""
    path = screen.q1.kernel_t64.SELECTED_MAP.resolve()
    raw = path.read_bytes()
    declared = json.loads(raw)
    rows = [int(value, 16) for value in declared['generator_rows_hex']]
    columns = declared['columns']
    identity = dict(bits=16, width=64, expansion_rows=list(map(hex, rows)), feedback_columns=columns)
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    macro = dict(physical_steps=2, physical_step_bits=64, macro_step_bits=128,
        physical_windows=16, macro_windows=32, state_continuity='retained_between_halves',
        initial_state='zero', flush=False)
    record = dict(schema='s16-selected-t64-fixed-maps-1', t=64, s=16,
        expansion_rows_hex=list(map(hex, rows)), feedback_columns=columns,
        source=dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest()),
        macro=macro, distribution='uniform_gl', map_sha256=digest,
        sampling='independent uniform GL16 for every physical t64 step', whole_code_certificate=False)
    physical = dict(bits=16, birth_density='capped', columns=columns,
                    map_images={1 << j: row for j, row in enumerate(rows)})
    data = dict(bits=16, physical_step_bits=64, windows=32, macro_step_bits=128,
                birth_density='capped', map_sha256=digest)
    return data, record, physical


class LengthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.beta, cls.counts = screen.exact_outer()

    def setUp(self):
        self.previous_precision = ctx.prec
        ctx.prec = 192
        self.data, self.map_record, self.physical = selected_fixture()

    def tearDown(self):
        ctx.prec = self.previous_precision

    def test_natural_geometry_and_exact_outer(self):
        for K in (4096, 65536, 69632, 98304, 131072, 1048576):
            geom = screen.geometry(K)
            self.assertEqual((geom.K, geom.N, geom.group_count, geom.regions,
                              geom.macros_per_region), (K, 2*K, K//128, 64, K//4096))
            self.assertEqual(geom.regions * geom.macros_per_region * 128, 2*K)
        self.assertEqual(self.beta, Q(2**256, 65535**8))
        self.assertEqual(sum(self.counts), 2**128 - 1)
        self.assertEqual(self.counts[0], 0)
        for invalid in (True, 0, -4096, 65537, 65536.0, '65536'):
            with self.assertRaises(ValueError):
                screen.geometry(invalid)

    def test_selected_record_authentication_needs_no_census(self):
        with patch.object(screen.q1.kernel_t64, 'prepare') as census:
            self.assertEqual(screen.authenticate_record(self.map_record), self.map_record['map_sha256'])
            census.assert_not_called()
        for key, value in (('s', 20), ('t', 128), ('map_sha256', 'stale'),
                           ('sampling', 'transvections'), ('whole_code_certificate', True)):
            bad = deepcopy(self.map_record)
            bad[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                screen.authenticate_record(bad)
        for entry, key, value in (('source', 'sha256', 'stale'), ('macro', 'flush', True)):
            bad = deepcopy(self.map_record)
            bad[entry][key] = value
            with self.assertRaises(ValueError):
                screen.authenticate_record(bad)

    def test_shared_prepare_binds_selected_basis_and_capped_density(self):
        with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical), \
                patch.object(screen.q1.kernel_t64, 'prepare') as census:
            data, record, sources = screen.prepare(192, self.data, self.map_record)
            self.assertIs(data, self.data)
            self.assertIs(record, self.map_record)
            self.assertIn(str(Path(screen.__file__).resolve()), sources)
            census.assert_not_called()
            for key, value in (('bits', 20), ('birth_density', 'classes'), ('map_sha256', 'stale')):
                bad = dict(self.data, **{key: value})
                with self.assertRaises(ValueError):
                    screen.prepare(192, bad, self.map_record)
            bad_physical = deepcopy(self.physical)
            bad_physical['map_images'][1] ^= 1
            with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=bad_physical), \
                    self.assertRaises(ValueError):
                screen.prepare(192, self.data, self.map_record)
        for data, record in ((self.data, None), (None, self.map_record)):
            with self.assertRaises(ValueError):
                screen.prepare(192, data, record)

    def test_options_fail_before_numerical_preparation(self):
        with patch.object(screen, 'prepare') as prepare:
            for occupancies in ([], [2], [769], [3, 3], [True], [3.0], None, '3'):
                with self.subTest(occupancies=occupancies), self.assertRaises(ValueError):
                    screen.run_tail(K=98304, occupancies=occupancies)
            for tilts in ([], ['0'], ['nan'], ['inf'], ['.01', '1/100']):
                with self.assertRaises(ValueError):
                    screen.run_tail(K=98304, occupancies=[3], tilts=tilts)
            for q in (True, 0, 3):
                with self.assertRaises(ValueError):
                    screen.run_low(q, K=98304)
            with tempfile.TemporaryDirectory() as folder:
                existing = Path(folder) / 'existing.json'
                existing.touch()
                with self.assertRaises(ValueError):
                    screen.run_tail(K=98304, occupancies=[3], output=existing)
            prepare.assert_not_called()

    def receipt(self, q, K=98304):
        geom = screen.geometry(K)
        return dict(schema=f'finite-packet-q{q}-diagnostic-1', geometry=asdict(geom),
            K=K, N=2*K, threshold=2*K//10, distance='1/10', map_record=self.map_record,
            zero_initial_state=True, final_flush=False, occupancy_covered=[q],
            count_sha256=screen.count_hash(self.counts, cumulative=q == 1),
            source_sha256=screen.source_snapshot(), whole_code_certificate=False,
            **{f'q{q}_upper': [1, -50]})

    def test_low_delegation_uses_actual_length_and_exact_counts(self):
        for q, evaluator in ((1, screen.q1), (2, screen.q2)):
            receipt = self.receipt(q)
            with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical), \
                    patch.object(evaluator, f'evaluate_q{q}', return_value=receipt) as evaluate:
                result = screen.run_low(q, K=98304, tilts=['.001'], data=self.data, map_record=self.map_record)
            args, kwargs = evaluate.call_args
            self.assertEqual(args, (self.counts,))
            self.assertEqual(kwargs['tilts'], ('1/1000',))
            self.assertIs(kwargs['data'], self.data)
            self.assertIsNone(kwargs['output'])
            if q == 1:
                self.assertEqual(kwargs['group_count'], 768)
                self.assertEqual(kwargs['epochs_per_region'], 24)
                self.assertEqual(kwargs['count_kind'], 'shells')
            else:
                self.assertEqual(kwargs['geometry'], screen.geometry(98304))
            self.assertEqual(result['length_component'], f'q{q}')
            self.assertFalse(result['whole_code_certificate'])

    def test_low_rejects_inconsistent_receipts(self):
        for key, bad in (('K', 65536), ('threshold', 13107), ('q1_upper', [0, 0]),
                         ('whole_code_certificate', True), ('source_sha256', {'stale': 'source'})):
            receipt = self.receipt(1)
            receipt[key] = bad
            with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical), \
                    patch.object(screen.q1, 'evaluate_q1', return_value=receipt), \
                    self.assertRaises((ValueError, ArithmeticError)):
                screen.run_low(1, K=98304, data=self.data, map_record=self.map_record)

    def test_uniform_formula_retains_choose_beta_cutoff_and_terminal_mass(self):
        # The mass halves in either state at every region; both terminal
        # coordinates count, hence2^-64, not2^-65.
        regional = [arb_mat([[arb(1)/4, arb(1)/4], [arb(1)/4, arb(1)/4]])] * 4
        for K in (4096, 98304, 1048576):
            upper = screen.occupancy_upper(regional, K=K, occupancy=3, beta=Q(7, 3), tilt='1/1000')
            expected = (arb(comb(K//128, 3)) * screen.q1.kernel_t64.aq(Q(7, 3))**3 *
                        (arb(2*K//10)/1000).exp() * arb(2)**-64)
            self.assertGreaterEqual(upper, expected.lower())
            self.assertLess(float(upper/expected), 1.000000001)

    def test_tail_records_exact_nonconsecutive_coverage_and_best_tilts(self):
        local = [arb_mat([[1]])] * 33
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / 'fresh.json'
            with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical), \
                    patch.object(screen.q1.kernel_t64, 'local_operators', return_value=local), \
                    patch.object(screen.q1, 'placement', return_value=local) as placement, \
                    redirect_stdout(StringIO()):
                result = screen.run_tail(K=98304, occupancies=[5, 3], tilts=['.002', '.001'],
                    data=self.data, map_record=self.map_record, output=output)
            self.assertEqual(result['schema'], 'rs16-length-tail-1')
            self.assertEqual(result['occupancy_covered'], [3, 5])
            self.assertEqual(result['occupancy_choices'], {'3': '1/1000', '5': '1/1000'})
            self.assertEqual(set(result['occupancy_uppers']), {'3', '5'})
            for call in placement.call_args_list:
                self.assertEqual(call.kwargs['epochs'], 24)
                self.assertEqual(call.kwargs['windows'], 32)
                self.assertEqual(call.kwargs['maximum_groups'], 5)
            self.assertFalse(result['whole_code_certificate'])
            self.assertEqual(json.loads(output.read_text()), result)
            self.assertGreaterEqual(screen.endpoint_value(result['union_upper']),
                sum(map(screen.endpoint_value, result['occupancy_uppers'].values())))

    def test_changed_sources_reject_without_output(self):
        local = [arb_mat([[1]])] * 33
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / 'never.json'
            with patch.object(screen.q1.kernel_t64, 'authenticate', return_value=self.physical), \
                    patch.object(screen.q1.kernel_t64, 'local_operators', return_value=local), \
                    patch.object(screen.q1, 'placement', return_value=local), \
                    patch.object(screen, 'source_snapshot', side_effect=[{'a': 'old'}, {'a': 'changed'}]), \
                    self.assertRaises(RuntimeError):
                screen.run_tail(K=98304, occupancies=[3], tilts=['.001'],
                    data=self.data, map_record=self.map_record, output=output)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
