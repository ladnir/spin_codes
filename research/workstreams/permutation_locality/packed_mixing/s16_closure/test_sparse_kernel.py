"""Small exhaustive checks of the uniform-GL fixed-occupancy envelope."""
from fractions import Fraction as Q
from itertools import combinations, product
import unittest

from flint import arb, ctx
import sparse_kernel as kernel
import refresh_kernel
import birth_classes_probe
import occupancy_kernel


def prepare(rows, columns):
    images = [0]*(1 << len(rows))
    for q in range(1, len(images)):
        bit = q & -q
        images[q] = images[q ^ bit] ^ rows[bit.bit_length()-1]
    data = refresh_kernel.prepare(images, columns, len(rows), 4)
    data = occupancy_kernel.prepare(data)
    data = birth_classes_probe.attach(data, images)
    del data['updates']
    return dict(data, distribution='uniform_gl'), images


def fixed_inputs(windows, occupancy):
    for positions in combinations(range(windows), occupancy):
        for labels in product(range(1, 16), repeat=occupancy):
            yield sum(label << (4*p) for p, label in zip(positions, labels))


def actual(data, images, states, occupancy, z):
    inputs = list(fixed_inputs(data['windows'], occupancy))
    S = len(images)
    result = [Q(0)]*S
    for q, prior in states:
        for x in inputs:
            feedback = 0
            for j, column in enumerate(data['columns']):
                if (x >> j) & 1:
                    feedback ^= column
            mass = prior*z**((x ^ images[q]).bit_count())/len(inputs)
            if q == 0:
                result[feedback] += mass
            else:
                for fresh in range(1, S):
                    result[fresh ^ feedback] += mass/(S-1)
    return result


class SparseKernelTests(unittest.TestCase):
    def setUp(self):
        ctx.prec = 192

    def test_all_small_states_and_occupancies(self):
        for rows, columns in (([1, 6, 12], [1, 2, 4, 3]),
                              ([0x13, 0x36, 0xC4], [1, 2, 4, 3, 5, 6, 7, 1])):
            data, images = prepare(rows, columns)
            z = Q(3, 5)
            matrices = kernel.outward_at_z(data, kernel.aq(z))
            L = len(images)-1
            for j, matrix in enumerate(matrices):
                zero = actual(data, images, [(0, Q(1))], j, z)
                self.assertTrue(matrix[0, 0] >= kernel.aq(zero[0]))
                for index, level in enumerate(data['birth_class_levels'], 3):
                    mass = sum(p for q, p in enumerate(zero) if q and images[q].bit_count() == level)
                    self.assertTrue(matrix[0, index] >= kernel.aq(mass))
                for q in range(1, len(images)):
                    exact = actual(data, images, [(q, Q(1))], j, z)
                    level = images[q].bit_count()
                    index = 3+list(data['birth_class_levels']).index(level)
                    for row in (1, index):
                        self.assertTrue(matrix[row, 0] >= kernel.aq(exact[0]))
                        self.assertTrue(all(matrix[row, 2]/L >= kernel.aq(p) for p in exact[1:]))
                exact = actual(data, images, [(q, Q(1, L)) for q in range(1, L+1)], j, z)
                self.assertTrue(matrix[2, 0] >= kernel.aq(exact[0]))
                self.assertTrue(all(matrix[2, 2]/L >= kernel.aq(p) for p in exact[1:]))

    def test_quiet_epoch_has_no_return_or_lazy_transition(self):
        data, _ = prepare([1, 6, 12], [1, 2, 4, 3])
        matrix = kernel.outward_at_z(data, arb(1))[0]
        self.assertEqual(matrix[0, 0], 1)
        for i in range(1, matrix.nrows()):
            self.assertEqual(matrix[i, 0], 0)
            self.assertEqual(matrix[i, 2], 1)
            self.assertTrue(all(matrix[i, j] == 0 for j in range(matrix.ncols()) if j != 2))

    def test_nonzero_feedback_cancellation_floor(self):
        data, _ = prepare([1, 2, 4, 8], [1, 2, 4, 8])
        matrix = kernel.outward_at_z(data, arb(1))[1]
        for i in range(1, matrix.nrows()):
            expected = kernel.aq(Q(1, 15))
            self.assertTrue(matrix[i, 0] >= expected)
            self.assertTrue(matrix[i, 0]-expected < arb(2)**-180)

    def test_explicit_distribution_and_input_validation(self):
        data, _ = prepare([1, 6, 12], [1, 2, 4, 3])
        for changed in (dict(data, distribution='transvections'), dict(data, distribution=None),
                        dict(data, updates=32), dict(data, bits=0), dict(data, windows=0)):
            with self.assertRaises(ValueError):
                kernel.outward_at_z(changed, arb(1))
        for z in (arb(0), arb(2), arb('nan')):
            with self.assertRaises(ValueError):
                kernel.outward_at_z(data, z)
        with self.assertRaises(ValueError):
            kernel.outward(data, 0)


if __name__ == '__main__':
    unittest.main()
