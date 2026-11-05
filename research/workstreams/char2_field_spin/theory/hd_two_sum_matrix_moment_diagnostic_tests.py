from itertools import product
from unittest import TestCase, main

import numpy as np

import hd_two_sum_matrix_moment as frozen
import hd_two_sum_matrix_moment_diagnostic as gate


class MatrixDiagnosticTests(TestCase):
    def test_maxplus_backtracking_matches_every_small_state_degree_word(self):
        matrix = np.log2(np.array([[[1., 2., 5.], [7., 3., 1.]],
                                  [[3., 9., 2.], [2., 4., 8.]]]))
        for target in range(7):
            best = -np.inf
            for states in product(range(2), repeat=3):
                path = (0, *states)
                for degrees in product(range(3), repeat=3):
                    if sum(degrees) == target:
                        value = sum(matrix[path[i], path[i+1], degrees[i]] for i in range(3))
                        best = max(best, value)
            actual = gate.largest_column_word(matrix, 3, target)
            self.assertAlmostEqual(actual['log2_term'], best)
            self.assertEqual(sum(row['degree'] for row in actual['columns']), target)
            self.assertEqual(actual['columns'][0]['initial'], 0)
            self.assertEqual(actual['columns'][-1]['final'], actual['final_state'])
            self.assertAlmostEqual(sum(row['coefficient_log2'] for row in actual['columns']), best)

    def test_regional_path_terms_sum_back_to_column(self):
        weights = np.log2([1., 2., 3., 4.])
        measure = dict(log2_probability=-2., log2_atom00=-1., log2_atom_cross=-.5,
                       log2_scale01=-1.5, log2_scale10=-4., log2_scale11=-3.)
        matrix, _ = frozen.column_matrix(weights, 2, measure, arity=2)
        for i, j, d in product(range(2), range(2), range(5)):
            row = gate.path_breakdown(weights, 2, measure, 2, i, j, d)
            self.assertAlmostEqual(row['coefficient_log2'], matrix[i, j, d])
            self.assertLessEqual(row['top_terms'][0]['log2_fraction_of_coefficient'], 1e-12)

    def test_guards(self):
        with self.assertRaises(ValueError):
            gate.run(max_seconds=31)
        with self.assertRaises(ValueError):
            gate.largest_column_word(np.zeros((2, 2, 2)), 2, 3)
        with self.assertRaises(ValueError):
            gate.largest_column_word(np.full((2, 2, 2), -np.inf), 2, 1)


if __name__ == '__main__':
    main()
