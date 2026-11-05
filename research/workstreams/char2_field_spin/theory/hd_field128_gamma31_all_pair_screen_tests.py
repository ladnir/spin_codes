"""Tiny map, local-oracle, and state30 replay checks; no distance gate."""

from fractions import Fraction as F
import unittest
from unittest.mock import patch

import numpy as np

import hd_field128_gamma31_all_pair_screen as candidate
from hd_blockwise_maps import matmul, mul, power


def xor_sum(values):
    result = 0
    for value in values:
        result ^= value
    return result


def inverse(matrix):
    size = len(matrix)
    rows = [list(row)+[int(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(row for row in range(column, size) if rows[row][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = power(rows[column][column], 254, 256)
        rows[column] = [mul(value, scale, 256) for value in rows[column]]
        for row in range(size):
            if row != column:
                scale = rows[row][column]
                rows[row] = [a^mul(scale, b, 256) for a, b in zip(rows[row], rows[column])]
    return tuple(tuple(row[size:]) for row in rows)


class Gamma31Tests(unittest.TestCase):
    def test_actual_byte_maps_ca_zero_unit_diagonal_and_refresh_inverse(self):
        a, c, h = candidate.fixed_maps()
        self.assertEqual((len(a), len(a[0]), len(h)), (64, 31, 31))
        self.assertFalse(any(value for row in matmul(c, a, 256) for value in row))
        self.assertTrue(all(xor_sum(row) == 1 for row in a))
        self.assertTrue(all(xor_sum(mul(x, x, 256) for x in row) == 1 for row in a))
        identity = tuple(tuple(int(i == j) for j in range(31)) for i in range(31))
        inv = inverse(h)
        self.assertEqual(tuple(map(tuple, matmul(h, inv, 256))), identity)
        self.assertEqual(tuple(map(tuple, matmul(inv, h, 256))), identity)
        for size in (1, 2, 3, 7, 15):
            minor = tuple(tuple(h[i][j] for j in range(31-size, 31)) for i in range(size))
            inverse(minor)
        self.assertEqual(candidate.fixed_maps(30), candidate.gamma.fixed_maps(30))

    def test_state31_exact_local_oracle_and_injective_quiet_cases(self):
        q = candidate.FIELD_ORDER
        for z in (F(1, 2), F(1, 8)):
            logs = candidate.local_logs(float(z))
            self.assertEqual(logs.shape, (65, 2, 2))
            for j in (0, 1, 30, 31, 32, 64):
                exact = candidate.exact_packet(j, z)
                np.testing.assert_allclose(logs[j].ravel(),
                    [candidate.old.log_fraction(value) for value in exact], atol=2e-9, rtol=0)
                if j <= 31:
                    self.assertEqual(exact[0], F(1, q**j))
            for j in (1, 31, 64):
                self.assertEqual(candidate.gamma.joint_probability(j, j, 0, q, 31), 0)
                self.assertEqual(candidate.gamma.joint_probability(j, 0, 31, q, 31), 0)

    def test_state30_local_coefficients_and_alias_replay(self):
        z, alpha, odds = .5, .125, ((-2., -3.),)
        old = candidate.frozen.AliasCoefficients(z, alpha, state=30)
        control = candidate.AliasCoefficients(z, alpha, state=30)
        np.testing.assert_array_equal(candidate.local_logs(z, state=30),
                                      candidate.gamma.local_logs(z, state=30))
        np.testing.assert_array_equal(control.single, old.single)
        np.testing.assert_array_equal(control.composed, old.composed)
        np.testing.assert_array_equal(control.pair_logs(odds[0]), old.pair_logs(odds[0]))
        previous = candidate.frozen.AliasBankScreen(z, alpha, odds, -8., state=30,
                                                    grid=64, distance=F(2, 5), _coefficients=old)
        replay = candidate.AliasBankScreen(z, alpha, odds, -8., state=30,
                                            grid=64, _coefficients=control)
        for weight in (1, 2, 3):
            np.testing.assert_array_equal(replay.regional(weight), previous.regional(weight))
        self.assertEqual(replay.constant, previous.constant)
        self.assertEqual(replay.witness()['construction'], previous.witness()['construction'])
        first = candidate.selected_screen(replay, 2)
        second = candidate.frozen.selected_screen(previous, 2, support_method='exact-profiles')
        self.assertEqual(first['log2_upper_screen'], second['log2_upper_screen'])
        self.assertEqual(first['largest_profile'], second['largest_profile'])

    def test_state31_layout_metadata_and_tiny_selected_api(self):
        coefficients = candidate.AliasCoefficients(.5, .125)
        self.assertEqual(coefficients.composed.shape, (4, 65, 65))
        self.assertTrue(np.isneginf(coefficients.single).all())
        bank = candidate.AliasBankScreen(.5, .125, ((-2., -3.),), -8., grid=64,
                                          _coefficients=coefficients)
        row = candidate.selected_screen(bank, 1)
        witness = bank.witness()
        for result in (row, witness):
            self.assertEqual(result['state_dimension'], 31)
            self.assertEqual(result['emission_minimum'], 34)
            self.assertEqual(result['refresh_diagonal_fullfield_products_per_packet'], 62)
            self.assertEqual(result['state_denominator'], str(candidate.FIELD_ORDER**31-1))
            self.assertEqual(result['refresh_density'],
                             str(F(candidate.FIELD_ORDER**31-1, (candidate.FIELD_ORDER-1)**31)))
            self.assertFalse(result['explicit_state30_control'])
            self.assertFalse(result['directed'])
            self.assertFalse(result['whole_code_certificate'])
            self.assertTrue(result['all_pair_before_alpha'])
            self.assertIn('hd_field128_gamma31_all_pair_screen.py', result['source_sha256'])
        self.assertFalse(row['common_potential_used'])
        self.assertEqual(row['profile_pairs_evaluated'], 1)
        self.assertNotIn(31, candidate.gamma.STATES)
        with self.assertRaisesRegex(RuntimeError, 'Pinned'):
            with patch.object(candidate, 'NOTE', candidate.NOTE.with_name('HD_FIELD128_GAMMA_RETURN.md')):
                bank.check_pins()

    def test_guards_and_no_frozen_whitelist_mutation(self):
        for state in (16, 24, 28, 32, 31.0, True):
            with self.assertRaises(ValueError):
                candidate.construction(state)
        for kwargs in (dict(q=256), dict(z=0), dict(alpha=0), dict(alpha=2)):
            args = dict(z=.5, alpha=.125)
            args.update(kwargs)
            with self.assertRaises(ValueError):
                candidate.AliasCoefficients(**args)
        coefficients = candidate.AliasCoefficients(.5, .125, state=30)
        with self.assertRaises(ValueError):
            candidate.AliasBankScreen(.5, .125, ((-2., -3.),), -8., _coefficients=coefficients)
        with self.assertRaises(TypeError):
            candidate.selected_screen(object(), 1)
        self.assertEqual(candidate.gamma.STATES, (16, 24, 28, 30))


if __name__ == '__main__':
    unittest.main()
