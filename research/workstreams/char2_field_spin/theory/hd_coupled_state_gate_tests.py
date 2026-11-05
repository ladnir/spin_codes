from fractions import Fraction as F
from itertools import product
from unittest import TestCase, main

import numpy as np

import hd_coupled_state_gate as gate


def bank_with(active, maximum=2):
    size = len(active)
    with np.errstate(divide='ignore'):
        identity = np.log2(np.eye(size))
    return {(m, f): identity if m == 0 else np.asarray(active)
            for m in range(maximum+1) for f in range(min(1, m//2)+1)}


class CoupledStateGateTests(TestCase):
    def test_log_lp_handles_tiny_weights_and_large_point_potential(self):
        matrix = np.full((3, 3), -np.inf)
        matrix[0, 1], matrix[1, 2], matrix[2, 0] = 0., -128., -8.
        answer = gate.propose_potential(bank_with(matrix), [0., -60., -60.])
        np.testing.assert_allclose(answer['log2_potential'], [0., -136/3, 112/3], atol=1e-10)
        self.assertAlmostEqual(answer['log2_lambda'], 44/3, places=10)
        self.assertAlmostEqual(answer['endpoint_log2'], 136/3, places=10)
        self.assertTrue(answer['actual_ratios_recomputed'])
        self.assertFalse(answer['raw_tiny_potentials_used_in_lp'])

    def test_edge_row_budget_suffices_but_full_ratio_is_recomputed(self):
        with np.errstate(divide='ignore'):
            matrix = np.log2([[.7, .1], [.01, .8]])
        bank, moments = bank_with(matrix), np.array([0., -3., -4.])
        objective, rows, rhs, bounds = gate.edge_budget_program(bank, moments)
        self.assertEqual(rows.shape[1], 4)
        self.assertEqual(bounds[0], (0., 0.))
        self.assertEqual(objective[-1], 1/256)
        answer = gate.propose_potential(bank, moments)
        self.assertLessEqual(answer['log2_lambda'], answer['lp_t']+1e-10)
        independent = gate.recompute_ratios(bank, moments, answer['log2_potential'])
        self.assertEqual(answer['log2_lambda'], independent['log2_lambda'])
        # Deliberately supply an inferior potential. No asserted LP tolerance
        # can hide its actual row violation in the independently recomputed cap.
        inferior = gate.recompute_ratios(bank, moments, [0., -100.])
        self.assertGreater(inferior['log2_lambda'], 90.)

    def test_actual_ratio_bound_preserves_noncommuting_chronology(self):
        first = [[F(1, 2), F(1, 3)], [F(1, 64), F(1, 4)]]
        second = [[F(1, 3), F(1, 5)], [F(1, 32), F(1, 3)]]
        bank = bank_with(np.log2(np.array(first, dtype=float)))
        bank[2, 1] = np.log2(np.array(second, dtype=float))
        answer = gate.recompute_ratios(bank, [0., 0., 0.], [0., -2.])
        potential = [F(1), F(1, 4)]
        exact_lambda = max(F(1), *(sum(row[j]*potential[j] for j in range(2))/potential[i]
            for matrix in (first, second) for i, row in enumerate(matrix)))
        self.assertAlmostEqual(2**answer['log2_lambda'], float(exact_lambda))
        for word in product((first, second), repeat=4):
            row = [F(1), F(0)]
            for matrix in word:
                row = [sum(row[i]*matrix[i][j] for i in range(2)) for j in range(2)]
            self.assertLessEqual(sum(row), F(4)*exact_lambda**4)

    def test_complete_histograms_and_nulls_preserve_exact_zero(self):
        bank = bank_with(np.array([[0., -np.inf], [-2., -3.]]))
        records = gate.serialize_bank(bank)
        self.assertEqual(len(records), 4)
        self.assertEqual(records[-1]['log2_matrix'][0], [0., None])
        restored = {(row['active_weight'], row['full_pairs']):
            np.array([[-np.inf if value is None else value for value in r]
                      for r in row['log2_matrix']]) for row in records}
        for key in bank:
            np.testing.assert_array_equal(bank[key], restored[key])

    def test_frozen_scalar_shape_and_zero_order(self):
        moments = gate.shape(18)
        self.assertEqual(moments[0], 0.)
        self.assertEqual(moments[18], gate.LOG_PROBABILITY+18*gate.LOG_ATOM)
        self.assertEqual(len({(m, f) for m in range(19) for f in range(min(1, m//2)+1)}), 36)

    def test_fail_closed_guards(self):
        with self.assertRaises(ValueError):
            gate.checked_bank({}, 18)
        with self.assertRaises(ValueError):
            gate.edge_budget_program(bank_with(np.zeros((2, 2))), [0., 1., 2.], log_bound=0.)
        with self.assertRaises(ValueError):
            gate.recompute_ratios(bank_with(np.zeros((2, 2))), [0., 1., 2.], [1., 0.])
        with self.assertRaises(ValueError):
            gate.run(max_seconds=31)


if __name__ == '__main__':
    main()
