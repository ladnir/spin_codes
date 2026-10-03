from itertools import product
import unittest

import numpy as np

import diagnose24
import occupancy24


class PathTests(unittest.TestCase):
    def test_boundary_path_statistics(self):
        matrix = np.array([[.3, .2, .1], [.4, .5, .2], [.1, .7, .3]])
        mass, best = 0., 0.
        transitions, terminal = np.zeros((3, 3)), np.zeros(3)
        for suffix in product(range(3), repeat=4):
            path = (0,)+suffix
            weight = np.prod([matrix[a, b] for a, b in zip(path, path[1:])])
            mass += weight
            best = max(best, weight)
            terminal[path[-1]] += weight
            for a, b in zip(path, path[1:]):
                transitions[a, b] += weight
        actual = diagnose24.path_statistics(np.log(matrix), 4)
        self.assertAlmostEqual(np.exp(actual['log_moment']), mass)
        self.assertAlmostEqual(np.exp(actual['maximum_path_log_moment']), best)
        np.testing.assert_allclose(actual['terminal_probabilities'], terminal/mass)
        np.testing.assert_allclose(actual['expected_boundary_transitions'], transitions/mass)
        costs = np.array([[0., 1., 2.], [3., 0., 4.], [0., 5., 1.]])
        marked = np.full((3, 3), -np.inf)
        marked[costs > 0] = np.log((matrix*costs)[costs > 0])
        derivative = occupancy24.marked_matrix_expectation(np.log(matrix), marked, regions=4)
        self.assertAlmostEqual(derivative, float((costs*transitions).sum()/mass))


if __name__ == '__main__':
    unittest.main()
