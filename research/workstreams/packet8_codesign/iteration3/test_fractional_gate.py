"""Independent path-power and hidden-position clipping checks."""
from itertools import combinations, product
from math import comb, exp
import unittest

import numpy as np

import fractional_gate as fractional


def evaluate(matrices):
    row = np.array([1., 0.])
    for matrix in matrices:
        row = row@matrix
    return row.sum()


class FractionalGateTests(unittest.TestCase):
    def test_positive_path_subadditivity(self):
        random = np.random.default_rng(719)
        for length in range(1, 6):
            matrices = [random.random((2, 2))*.3 for _ in range(length)]
            for alpha in (.25, .4, .7, 1.):
                self.assertLessEqual(evaluate(matrices)**alpha,
                    evaluate([m**alpha for m in matrices])+1e-14)

    def test_power_is_after_label_mixture(self):
        local = np.array([[[1., 0.], [0., .2]], [[.02, .04], [.08, .16]],
                          [[.1, .3], [.02, .25]]])
        expected = fractional.marked.potential_operators(local, 0., packet_bits=2)**.4
        actual = fractional.fractional_operators(local, .4, packet_bits=2)
        self.assertTrue(np.array_equal(actual, expected))
        wrong = fractional.marked.potential_operators(local**.4, 0., packet_bits=2)
        self.assertGreater(np.max(np.abs(actual-wrong)), .001)
        self.assertTrue(np.array_equal(fractional.fractional_operators(local, 1., packet_bits=2),
            fractional.marked.potential_operators(local, 0., packet_bits=2)))

    def test_whole_conditional_clipping_with_hidden_positions(self):
        # There are two two-slot physical steps. The j=1 matrix genuinely
        # depends on which slot is occupied; averaging is conditional on J.
        matrices = {
            (): np.array([[.8, .0], [.0, .3]]),
            (0,): np.array([[.15, .05], [.02, .08]]),
            (1,): np.array([[.04, .25], [.18, .07]]),
            (0, 1): np.array([[.08, .02], [.09, .04]])}
        local = [matrices[()], (matrices[(0,)]+matrices[(1,)])/2, matrices[(0, 1)]]
        placements = list(combinations(range(4), 2))
        records = []
        for placement in placements:
            subsets = [tuple(x-2*i for x in placement if x//2 == i) for i in range(2)]
            records.append((tuple(map(len, subsets)), evaluate([matrices[s] for s in subsets])))
        subset_union = comb(4, 2)
        for message_factor in (2., 10., 40.):
            direct_clipped = sum(min(1., message_factor*v) for _, v in records)/len(records)
            conditional_clipped = sum(min(1., message_factor*evaluate([local[j] for j in counts]))
                for counts, _ in records)/len(records)
            self.assertLessEqual(direct_clipped, conditional_clipped+1e-14)
            for alpha in (.25, .4, .7, 1.):
                powered = sum(evaluate([local[j]**alpha for j in counts])
                    for counts, _ in records)/len(records)
                # C(4,2) is outside alpha; only the conditional message
                # factor and its conditional moment receive that power.
                self.assertLessEqual(subset_union*direct_clipped,
                    subset_union*message_factor**alpha*powered+1e-13)

    def test_good_event_marker_and_fractional_majorant(self):
        local = [np.array([[.7, 0.], [0., .4]]), np.array([[.1, .2], [.05, .3]]),
                 np.array([[.05, .1], [.08, .12]])]
        alpha, nu, h, message_factor = .55, 1.2, 3, 7.
        for counts in product(range(3), repeat=3):
            cost = sum(min(j, 2) for j in counts)
            direct = min(1., message_factor*evaluate([local[j] for j in counts])) if cost >= h else 0.
            majorant = message_factor**alpha*exp(nu*(cost-h))*evaluate([local[j]**alpha for j in counts])
            self.assertLessEqual(direct, majorant+1e-13)

    def test_invalid_powers(self):
        local = np.ones((3, 2, 2))*.2
        for alpha in (0., -1., 1.01, float('nan')):
            with self.assertRaises(ValueError):
                fractional.fractional_operators(local, alpha)


if __name__ == '__main__':
    unittest.main()
