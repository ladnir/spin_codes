"""Tiny controls; root executes remotely. This file generates no artifacts."""
from fractions import Fraction as F
from math import log2
import unittest
from unittest.mock import patch

import numpy as np
import hd_field128_wide_state_logs_v62 as fast


class Tests(unittest.TestCase):
    def test_integer_logs_and_exact_zeros(self):
        self.assertEqual(fast._log_exact(0), -np.inf)
        self.assertEqual(fast._log_exact(F(1, 1 << 50000)), -50000.)
        self.assertEqual(fast._log_exact(F(1 << 50000, 1)), 50000.)
        for value in (-1, .5, True):
            with self.assertRaises(ValueError):
                fast._log_exact(value)

    def test_all_tiny_rows_match_exact_packets(self):
        # Both sides use the SAME caps, not a claim that the caps are attained.
        # Includes s>n/2, nontrivial clipping/tails and both endpoint regimes.
        for q, n, state in ((16, 3, 2), (16, 5, 3), (16, 6, 2)):
            table = fast._compile(q, n, state)
            for z in (F(1, 16), F(1, 2), F(15, 16)):
                got = fast._evaluate(table, float(z)).reshape(n+1, 4)
                d = q**state-1
                rho = F(d, (q-1)**state)
                for j in range(n+1):
                    quiet, birth, ret, live = fast.exact._moments(j, z, q, n, state, True)
                    want = (quiet, rho*birth, ret/d,
                            min(rho*live, z**(n-state+1)) if j == 0 else rho*live)
                    for entry, value in enumerate(want):
                        expected = fast._log_exact(value)
                        if value:
                            self.assertAlmostEqual(got[j, entry], expected, delta=2e-11)
                        else:
                            self.assertEqual(got[j, entry], -np.inf)
                self.assertFalse(got.flags.writeable)

    def test_tail_identity_and_bernstein_shapes(self):
        table = fast._compile(16, 5, 3)
        for j in range(6):
            bounds = fast.exact._factorial_bounds(j, 16, 5, 3)
            for z in (F(1, 4), F(3, 4)):
                result = fast.exact.tail.evaluate_envelope(bounds, z)
                lz, lw = log2(z), log2(1-z)
                basis = (5-np.arange(6))*lz+np.arange(6)*lw
                actual = (np.logaddexp2.reduce(table.original[j]+basis),
                          np.logaddexp2.reduce(table.clipped[j]+basis),
                          np.logaddexp2(5*lz+table.original[j, 0], lw+
                              np.logaddexp2.reduce(table.tails[j, 1:]
                                                  +(5-np.arange(1, 6))*lz)))
                for got, key in zip(actual, ('original_bernstein', 'clipped_bernstein', 'tail')):
                    expected = fast._log_exact(result[key])
                    if result[key]:
                        self.assertAlmostEqual(got, expected, delta=2e-11)
                    else:
                        self.assertEqual(got, -np.inf)

    def test_extreme_float_marks_do_not_underflow_positive_entries(self):
        table = fast._compile(16, 3, 2)
        for z in (float.fromhex('0x0.0000000000001p-1022'),
                  float.fromhex('0x1.fffffffffffffp-1')):
            got = fast._evaluate(table, z).reshape(4, 4)
            self.assertEqual(got[0, 1], -np.inf)
            self.assertEqual(got[0, 2], -np.inf)
            self.assertTrue(np.all(np.isfinite(got[1:])))
            self.assertTrue(np.all(np.isfinite(got[0, (0, 3)])))

    def test_public_geometry_and_warm_path(self):
        table = fast._compile(16, 3, 2)
        with patch.object(fast, '_prepared', return_value=table) as prepared, \
             patch.object(fast.exact, 'exact_packet', side_effect=AssertionError('full packet')), \
             patch.object(fast.exact, 'exact_moments', side_effect=AssertionError('exact at z')), \
             patch.object(fast, 'source_pins', side_effect=AssertionError('hot source IO')):
            for state in fast.STATES:
                got = fast.local_logs(.5, state)
                self.assertEqual(got.shape, (4, 2, 2))  # injected tiny table only
            self.assertEqual(prepared.call_count, 3)
            for state in (True, 31, 33, 64):
                with self.assertRaises(ValueError):
                    fast.local_logs(.5, state)
            for z in (0, 1, -1, np.nan, np.inf, True, np.bool_(False)):
                with self.assertRaises(ValueError):
                    fast.local_logs(z, 40)

    def test_preparation_closure_guards_and_scope(self):
        fast._prepared.cache_clear()
        fake_table = object()
        with patch.object(fast, 'source_pins', return_value={'test': 'pin'}) as pins, \
             patch.object(fast, 'check_source_pins') as check, \
             patch.object(fast, '_compile', return_value=fake_table) as compile_table:
            self.assertIs(fast._prepared(40), fake_table)
            self.assertIs(fast._prepared(40), fake_table)
            self.assertEqual(pins.call_count, 1)
            check.assert_called_once_with({'test': 'pin'})
            compile_table.assert_called_once_with(1 << 128, 64, 40)
        fast._prepared.cache_clear()
        with patch.object(fast, 'source_pins', return_value={'test': 'pin'}), \
             patch.object(fast, 'check_source_pins') as check, \
             patch.object(fast, '_compile', side_effect=ValueError('preparation failed')):
            with self.assertRaises(ValueError):
                fast._prepared(40)
            check.assert_called_once_with({'test': 'pin'})
        with patch.object(fast, 'source_pins', return_value={'test': 'pin'}):
            meta = fast.metadata()
            with self.assertRaises(RuntimeError):
                fast.check_source_pins({'changed': 'pin'})
        self.assertFalse(meta['directed'])
        self.assertFalse(meta['whole_code_certificate'])
        self.assertFalse(meta['construction_changed'])


if __name__ == '__main__':
    unittest.main()
