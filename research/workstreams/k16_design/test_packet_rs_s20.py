"""Light driver tests using synthetic matrices, never the S20 state census."""
from contextlib import ExitStack, redirect_stdout
from copy import deepcopy
from fractions import Fraction as Q
from io import StringIO
import json
from math import comb
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat
import packet_rs_s20 as screen


def rational(value):
    value = Q(value)
    return fmpq(value.numerator, value.denominator)


def endpoint_fraction(pair):
    return Q(pair[0]) * Q(2) ** pair[1]


class S20DriverTests(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 192
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        _, _, record = screen.inner.construction()
        record.update(full_state_census=True, state_count=2**20,
                      exact_profiles_regenerated=True, birth_density='classes')
        self.record = record
        self.data = dict(bits=20, physical_step_bits=64, macro_windows=32,
                         map_sha256=record['map_sha256'], birth_density='classes')
        self.physical = dict(bits=20, distribution='uniform_gl')

    def tearDown(self):
        ctx.prec = self.precision

    def options(self, **updates):
        args = dict(occupancies=[3, 5], tilts=['.01'], precision=192,
                    output=None, method='exact', markers=())
        return screen.options(**(args | updates))

    def test_options_require_explicit_increasing_coverage(self):
        self.assertEqual(self.options()[0], (3, 5))
        self.assertEqual(self.options(occupancies=[8192])[0], (8192,))
        for values in ([], [2], [8193], [5, 3], [3, 3], [True], [3.0]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                self.options(occupancies=values)
        for precision in (True, 127, 192.0):
            with self.assertRaises(ValueError):
                self.options(precision=precision)

    def test_options_reject_invalid_tilts_and_marker_laws(self):
        for values in ([], ['0'], ['-1'], ['.5', '1/2'], '0.1', ['nan']):
            with self.subTest(values=values), self.assertRaises(ValueError):
                self.options(tilts=values)
        for values in ([], ['0'], ['-1'], ['2'], ['1/2', '.5'], ['1'], '0.1'):
            with self.subTest(values=values), self.assertRaises(ValueError):
                self.options(method='fugacity', markers=values)
        self.assertEqual(self.options(occupancies=[8192], method='fugacity', markers=['1'])[2], (Q(1),))
        with self.assertRaises(ValueError):
            self.options(method='unrecognized')

    def test_output_and_partial_preparation_fail_before_census(self):
        output = self.folder / 'existing.json'
        output.touch()
        with patch.object(screen.inner, 'prepare') as prepare:
            for kwargs in (dict(output=output), dict(data=self.data), dict(map_record=self.record)):
                with self.subTest(kwargs=tuple(kwargs)), self.assertRaises(ValueError):
                    screen.run_tail(occupancies=[3], tilts=['.01'], **kwargs)
            prepare.assert_not_called()

    def test_authentication_calls_exact_constructor_and_rejects_s16(self):
        with patch.object(screen.inner, 'authenticate', return_value=self.physical) as authenticate:
            self.assertIs(screen.authenticate(self.data, self.record), self.physical)
            authenticate.assert_called_once_with(self.data)
        bad_cases = [('data', 'bits', 16), ('data', 'physical_step_bits', 128),
                     ('data', 'macro_windows', 16), ('data', 'map_sha256', 'foreign'),
                     ('record', 's', 16), ('record', 't', 128),
                     ('record', 'map_sha256', 'foreign'), ('record', 'distribution', 'transvections'),
                     ('record', 'full_state_census', False), ('record', 'state_count', 2**16)]
        for where, key, value in bad_cases:
            data, record = deepcopy(self.data), deepcopy(self.record)
            (data if where == 'data' else record)[key] = value
            with self.subTest(where=where, key=key), \
                    patch.object(screen.inner, 'authenticate', return_value=self.physical), \
                    self.assertRaises(ValueError):
                screen.authenticate(data, record)
        with patch.object(screen.inner, 'authenticate', side_effect=ValueError('foreign map')), \
                self.assertRaises(ValueError):
            screen.authenticate(self.data, self.record)

    def test_authentication_rejects_changed_base_or_constructor_pin(self):
        record = deepcopy(self.record)
        record['source']['sha256'] = '0' * 64
        with patch.object(screen.inner, 'authenticate', return_value=self.physical), \
                self.assertRaises(ValueError):
            screen.authenticate(self.data, record)
        record = deepcopy(self.record)
        record['source_sha256'][str(Path(screen.inner.__file__).resolve())] = '0' * 64
        with patch.object(screen.inner, 'authenticate', return_value=self.physical), \
                self.assertRaises(RuntimeError):
            screen.authenticate(self.data, record)

    def test_generic_return_atom_uses_twenty_bits_without_state_census(self):
        # Synthetic unit-moment profiles isolate the refresh denominator.
        # No claim is made that this one-profile fixture is the selected map.
        kernel = screen.q1.kernel_t64.sparse_kernel
        states = 2**20 - 1
        data = dict(distribution='uniform_gl', bits=20, windows=1,
            birth_class_levels=[1], image_histogram_weights=[1],
            histograms=[(0, 1, 0, 0, 0)], histogram_multiplicities=[states],
            zero_probabilities=[Q(1), Q(0)])
        with patch.object(kernel.births, 'class_masses',
                return_value=[[arb(1), arb(0)], [arb(0), arb(1)]]), \
                patch.object(kernel.occupancy, 'polynomial', return_value=[arb(1), arb(1)]):
            local = kernel.outward_at_z(data, arb(1) / 2)
        for source in (1, 2, 3):
            self.assertEqual(local[0][source, 0], arb(0))
            returned = endpoint_fraction(screen.q1.endpoint(local[1][source, 0]))
            self.assertGreaterEqual(returned, Q(1, states))
            self.assertLess(returned - Q(1, states), Q(1, 2**180))
            self.assertLess(returned, Q(1, 2**16 - 1))

    def toy_local(self):
        # Positive, noncommuting occupancy operators exercise state continuity.
        exact = [fmpq_mat([[rational(Q(j + 2, 16)), rational(Q(1, 16))],
                           [rational(Q(1, 16)), rational(Q(3, 16))]]) for j in range(33)]
        outward = [arb_mat([[screen.q1.kernel_t64.aq(Q(str(value[i, j])))
                             for j in range(2)] for i in range(2)]) for value in exact]
        self.assertNotEqual(exact[0] * exact[1], exact[1] * exact[0])
        return exact, outward

    def run_toy(self, *, occupancies, method='exact', markers=(), tilts=('.01',),
                output=None, fresh=False, source_change=False):
        exact, outward = self.toy_local()
        geometry = screen.q1.Geometry(32, 64, 128)
        counts = (Q(0), Q(1))  # Synthetic outer count metadata, not a code claim.
        sources = {'synthetic-source': 'fixed'}
        with ExitStack() as stack:
            stack.enter_context(patch.object(screen, 'GEOMETRY', geometry))
            stack.enter_context(patch.object(screen, 'authenticate', return_value=self.physical))
            prepared = stack.enter_context(patch.object(screen.inner, 'prepare',
                return_value=(self.data, self.record)))
            stack.enter_context(patch.object(screen.uniform, 'uniform_envelope', return_value=(Q(7, 3), counts)))
            stack.enter_context(patch.object(screen.q1.kernel_t64, 'local_operators', return_value=outward))
            stack.enter_context(patch.object(screen.q1, 'source_snapshot', return_value=sources,
                side_effect=[sources, {'changed': 'source'}] if source_change else None))
            stack.enter_context(redirect_stdout(StringIO()))
            shared = {} if fresh else dict(data=self.data, map_record=self.record)
            result = screen.run_tail(occupancies=occupancies, method=method,
                marker_probabilities=markers, tilts=tilts, output=output, **shared)
        return result, geometry, exact, prepared

    def assert_expected_upper(self, endpoint, expected):
        upper = screen.q1.kernel_t64.aq(endpoint_fraction(endpoint))
        self.assertGreaterEqual(endpoint_fraction(endpoint),
                                endpoint_fraction(screen.q1.endpoint(expected.lower())))
        self.assertTrue((upper / expected - 1).abs_upper() < arb(2)**-150)

    def test_exact_formula_preserves_state_and_all_terminal_mass(self):
        result, geometry, local, prepared = self.run_toy(occupancies=[3, 5])
        prepared.assert_not_called()
        for q in (3, 5):
            regional = screen.uniform.regional_uniform(local, q, matrix=fmpq_mat,
                rational=rational, rounding=lambda value: value)
            power = regional**geometry.regions
            moment = Q(str(sum(power[0, j] for j in range(2))))
            reset = Q(str(sum(regional[0, j] for j in range(2))))**geometry.regions
            self.assertNotEqual(moment, reset)
            expected = (screen.q1.kernel_t64.aq(comb(32, q) * Q(7, 3)**q * moment)
                        * (screen.q1.kernel_t64.aq(Q(1, 100)) * (geometry.N // 10)).exp())
            self.assert_expected_upper(result['occupancy_uppers'][str(q)], expected)
        self.assertEqual(result['occupancy_values'], [3, 5])
        self.assertFalse(result['contiguous_occupancy_interval'])
        self.assertIsNone(result['occupancy_covered'])
        self.assertEqual(set(result['occupancy_uppers']), {'3', '5'})

    def test_fugacity_formula_and_certain_event_endpoint(self):
        result, geometry, local, _ = self.run_toy(occupancies=[3, 32], method='fugacity', markers=['1/4', '1'])
        for q in (3, 32):
            p = Q(result['occupancy_choices'][str(q)]['marker_probability'])
            if q == 3:
                self.assertEqual(p, Q(1, 4))
            mixed = screen.density.iid_macro(local, Q(15, 16) * p,
                matrix=fmpq_mat, rational=rational, rounding=lambda value: value)
            power = mixed**(geometry.N // 128)
            moment = Q(str(sum(power[0, j] for j in range(2))))
            loss = 1 / (comb(32, q) * p**q * (1 - p)**(32 - q))
            expected = (screen.q1.kernel_t64.aq(comb(32, q) * Q(7, 3)**q * loss**64 * moment)
                        * (screen.q1.kernel_t64.aq(Q(1, 100)) * (geometry.N // 10)).exp())
            self.assert_expected_upper(result['occupancy_uppers'][str(q)], expected)
        # The endpoint-only law has exactly D_L=1 and cannot evaluate q<L.
        endpoint, _, _, _ = self.run_toy(occupancies=[32], method='fugacity', markers=['1'])
        self.assertEqual(endpoint['occupancy_choices']['32']['marker_probability'], '1')

    def test_receipt_keeps_s20_scope_and_only_explicit_coverage(self):
        output = self.folder / 'toy.json'
        result, _, _, prepare = self.run_toy(occupancies=[3, 4], output=output, fresh=True)
        prepare.assert_called_once_with(birth_density='classes')
        self.assertEqual(result['schema'], 'rs16-s20-occupancy-component-1')
        self.assertEqual(result['state_bits'], 20)
        self.assertEqual(result['return_denominator'], 2**20 - 1)
        self.assertEqual(result['birth_density'], 'classes')
        self.assertEqual(result['outer_symbol_update'], 'independent uniform GL16')
        self.assertEqual(result['inner_state_update'], 'independent uniform GL20')
        self.assertEqual(result['occupancy_covered'], [3, 4])
        self.assertTrue(result['contiguous_occupancy_interval'])
        self.assertTrue(result['every_requested_occupancy_checked'])
        self.assertTrue(result['fresh_computation'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['final_flush'])
        self.assertEqual(json.loads(output.read_text()), result)
        self.assertGreaterEqual(endpoint_fraction(result['union_upper']),
            sum(endpoint_fraction(value) for value in result['occupancy_uppers'].values()))

    def test_per_occupancy_minimum_uses_only_valid_tilts(self):
        result, _, _, _ = self.run_toy(occupancies=[3, 5], tilts=('.02', '.01'))
        self.assertEqual(result['occupancy_choices'],
                         {'3': {'tilt': '1/100'}, '5': {'tilt': '1/100'}})
        self.assertEqual(len(result['trials']), 2)

    def test_source_change_rejects_receipt(self):
        output = self.folder / 'must-not-exist.json'
        with self.assertRaises(RuntimeError):
            self.run_toy(occupancies=[3], output=output, source_change=True)
        self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
