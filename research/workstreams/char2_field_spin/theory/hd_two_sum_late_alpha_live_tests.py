from itertools import product
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_two_sum_late_alpha_live as gate


class QuietLiveTests(TestCase):
    def test_log_matrix_product_preserves_order_and_zeros(self):
        left = np.array([[1., 2.], [3., 0.]])
        right = np.array([[2., 0.], [1., 4.]])
        with np.errstate(divide='ignore'):
            logs = gate.log_matrix_product(np.log2(left), np.log2(right))
        np.testing.assert_allclose(2**logs, left@right, rtol=1e-14)
        self.assertFalse(np.array_equal(left@right, right@left))
        identity = gate.log_matrix_power(logs, 0)
        np.testing.assert_array_equal(2**identity, np.eye(2))

    def test_all_noncommuting_words_are_dominated(self):
        matrices = [np.array([[1., 0.], [0., .2]]),
                    np.array([[.2, .3], [.1, .4]]),
                    np.array([[.1, .04], [.3, .07]])]
        moments = np.array([1., .5, .3])
        with np.errstate(divide='ignore'):
            curves = np.log2(np.array([matrix.ravel() for matrix in matrices]).T)
        log_bound, _ = gate.shared_matrix(curves, np.log2(moments))
        bound = 2**log_bound
        for word in product(range(3), repeat=4):
            actual, scalar = np.eye(2), 1.
            for m in word:
                actual = actual@matrices[m]
                scalar *= moments[m]
            upper = scalar*np.linalg.matrix_power(bound, 4)
            self.assertTrue(np.all(actual <= upper+1e-14))
            self.assertLessEqual(actual[0].sum(), upper[0].sum()+1e-14)

    def test_ratio_uses_m0_and_every_later_occupancy(self):
        curves = np.array([[0., -4., -6., -7.], [-np.inf, -8., -5., -8.],
                           [-np.inf, -2., -9., -9.], [-3., -9., -9., -5.]])
        moments = np.array([0., -2., -4., -6.])
        bound, contacts = gate.shared_matrix(curves, moments)
        np.testing.assert_array_equal(bound, [[0., -1.], [0., 1.]])
        self.assertEqual(contacts, [[0], [2], [1], [3]])

    def test_terminal_vector_is_full_reference_mass(self):
        matrix = np.array([[1., 2.], [0., 1.]])
        with np.errstate(divide='ignore'):
            logs = gate.log_matrix_power(np.log2(matrix), 5)
        self.assertAlmostEqual(float(np.logaddexp2.reduce(logs[0])), np.log2(11))
        # Normalizing rows or retaining only the final-zero entry changes this.
        self.assertEqual(2**logs[0, 0], 1)

    def test_guards(self):
        with self.assertRaises(ValueError):
            gate.shared_matrix(np.zeros((4, 3)), np.ones(3))
        with self.assertRaises(ValueError):
            gate.log_matrix_power(np.zeros((2, 2)), -1)
        with self.assertRaises(ValueError):
            gate.run(max_seconds=31)

    def test_one_mock_bank_saved_all_four_curves_and_partial_scope(self):
        curve = -np.arange(19, dtype=float)
        bank = {(m, f): np.array([[curve[m], -100.], [-101., curve[m]-2.]])
                for m in range(19) for f in range(m//2+1)}
        source = dict(rowset_count='12870^18', aggregate_input_weight_floor=1296)
        word = dict(comparison=dict(direct_first_moment_alpha1_log2=-432.))
        with patch.object(gate, 'load_inputs', return_value=(source, word, curve, 7)):
            with patch.object(gate, 'forced_nonzero_regions', return_value=bank) as provider:
                answer = gate.run()
        self.assertEqual(provider.call_count, 1)
        self.assertEqual(np.asarray(answer['regional_log2_curves']).shape, (4, 19))
        self.assertEqual(answer['matrix_entry_order'], ['00', '01', '10', '11'])
        self.assertEqual(answer['regional_matrix_count'], 256)
        self.assertTrue(answer['all_regional_state_paths_included'])
        self.assertFalse(answer['nonquiet_outer_moments_included'])
        self.assertFalse(answer['whole_code_certificate'])


if __name__ == '__main__':
    main()
