"""Independent finite-field and positive-arithmetic checks of the new transfer."""
import itertools
import math
import unittest

import numpy as np

import occupation_refresh_v1 as refresh
from test_activation_q1 import gf4_mul


class RefreshTests(unittest.TestCase):
    def setUp(self):
        self.model = refresh.Epochs(4, 2, {2: 2, 4: 1}, [1, 0, 2, 0, 1])

    @staticmethod
    def actual(j, z):
        maps = [0, 15, 10, 5]
        matrix = np.zeros((4, 4))
        for positions in itertools.combinations(range(4), j):
            u = sum(1 << p for p in positions)
            syndrome = ((u & 15).bit_count() % 2) | (((u & 10).bit_count() % 2) << 1)
            for x in range(4):
                for alpha in (1, 2, 3):
                    matrix[x, gf4_mul(alpha, x) ^ syndrome] += z**((u ^ maps[x]).bit_count())/(3*math.comb(4, j))
        return matrix

    def test_every_occupation_and_class_against_exact_field(self):
        # The extreme distributions of L put equal mass on M-1 states.
        classes = [[np.array([1., 0., 0., 0.])],
                   [np.eye(4)[i] for i in (1, 2, 3)],
                   [np.array([0., 1/3, 1/3, 1/3])],
                   [np.array([0.]+[0. if j == i else .5 for j in (1, 2, 3)]) for i in (1, 2, 3)]]
        for z in (.15, .6, .99):
            envelope = np.exp(self.model.at(-math.log(z)))
            for j in range(5):
                actual = self.actual(j, z)
                for values in itertools.product((0., .3, 1.), repeat=4):
                    future = np.array(values)
                    potentials = np.array([max(d@future for d in group) for group in classes])
                    for state, group in enumerate(classes):
                        exact = max(d@actual@future for d in group)
                        self.assertLessEqual(exact, envelope[j, state]@potentials+2e-14)

    def test_region_against_direct_positive_polynomial(self):
        epoch = self.model.at(.7, 4)
        raw = np.exp(epoch)*np.array([math.comb(4, j) for j in range(5)])[:, None, None]
        direct = np.zeros((5, 4, 4))
        for a in range(5):
            for b in range(5-a):
                direct[a+b] += raw[a]@raw[b]
        direct /= np.array([math.comb(8, j) for j in range(5)])[:, None, None]
        np.testing.assert_allclose(np.exp(refresh.region_logs(epoch, 4, 8, 4)), direct, rtol=2e-14, atol=1e-16)

    def test_repeated_exact_epochs_with_zero_syndromes(self):
        envelope = np.exp(self.model.at(.5))
        for js in itertools.product(range(5), repeat=3):
            exact = np.eye(4)
            upper = np.eye(4)
            for j in js:
                exact = exact@self.actual(j, math.exp(-.5))
                upper = upper@envelope[j]
            self.assertLessEqual(sum(exact[0]), sum(upper[0])+2e-14)

    def test_mixture_and_terminal_against_direct_products(self):
        epoch = self.model.at(.6)
        theta = .27
        direct = sum(math.comb(4, j)*theta**j*(1-theta)**(4-j)*np.exp(epoch[j]) for j in range(5))
        matrix = refresh.epoch_mixture(epoch, 4, [theta])
        np.testing.assert_allclose(np.exp(matrix[0]), direct, rtol=2e-14)
        self.assertAlmostEqual(math.exp(refresh.terminal_logs(matrix, 7)[0]), sum(np.linalg.matrix_power(direct, 7)[0]), places=13)

    def test_dense_selected_boxes_cover_all_integer_types(self):
        model = refresh.DenseModel({2: 2, 4: 1}, 4, 4, 2, {2: 2, 4: 1},
                                   [1, 0, 2, 0, 1], 8, [-1., 0.], bands=[[2], [4]])
        result = model.search(minimum=3, maximum_nodes=15, target_bits=20)
        for zero in range(6):
            for low in range(9-zero):
                counts = np.array([zero, low, 8-zero-low])
                covers = [box for box in result['selected_boxes']
                          if np.all(counts >= box['lower']) and np.all(counts <= box['upper'])]
                self.assertEqual(len(covers), 1)
        combined = float(np.logaddexp.reduce([b['own_log_bound'] for b in result['selected_boxes']]))
        self.assertAlmostEqual(combined, result['log_union_upper'], places=11)
        # Re-evaluate each reported fixed witness directly; the proposal
        # interpolation table is deliberately not used in this replay.
        for box in result['selected_boxes']:
            witness = box['witness']
            epoch = self.model.at(math.exp(witness['log_surprisal']))
            theta = np.array(witness['proposal'])@witness['probabilities']
            moment = refresh.terminal_logs(refresh.epoch_mixture(epoch, 4, [theta]), 8)[0]
            corners = refresh.typed.vertices(box['lower'], box['upper'], 8)
            values = refresh.typed.point_logs(corners, 8, 4, witness['log_density_costs'],
                                               witness['proposal'], moment, 3,
                                               math.exp(witness['log_surprisal']))
            expected = float(max(values))+refresh.typed.lattice_log_count(box['lower'], box['upper'])
            self.assertAlmostEqual(expected, box['own_log_bound'], places=11)


if __name__ == '__main__':
    unittest.main()
