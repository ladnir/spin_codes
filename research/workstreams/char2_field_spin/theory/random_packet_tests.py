"""Independent exact local-law and normalization checks; no timings."""

from fractions import Fraction
import itertools
import math
import unittest

import numpy as np

import random_packet_screen as packet
from field_rate_tradeoff import source_moment


def binary_dot(a, b):
    return (a & b).bit_count() % 2


def binary_toeplitz(rows, columns, diagonals, vector):
    return sum((sum(((diagonals >> (i-j+columns-1)) & 1)*((vector >> j) & 1)
                    for j in range(columns)) % 2) << i for i in range(rows))


class RandomPacketTests(unittest.TestCase):
    def test_rectangular_toeplitz_projection_uniform(self):
        for rows, columns in itertools.product(range(1, 4), repeat=2):
            for vector in range(1, 1 << columns):
                counts = [0]*(1 << rows)
                for diagonals in range(1 << (rows+columns-1)):
                    counts[binary_toeplitz(rows, columns, diagonals, vector)] += 1
                self.assertEqual(counts, [1 << (columns-1)]*(1 << rows))

    def test_binary_packet_exhaustive(self):
        z = Fraction(2, 3)
        # B ranges over GL(2,2), or is fixed identity. A and C range over
        # every 2x1 and 1x2 matrix; these are also rectangular Toeplitz.
        invertibles = [(a, b) for a in range(1, 4) for b in range(1, 4) if a != b]
        for identity in (False, True):
            bases = [(1, 2)] if identity else invertibles
            for x in range(4):
                for incoming in (0, 1):
                    weights = [Fraction(0), Fraction(0)]
                    for basis, a, c in itertools.product(bases, range(4), range(4)):
                        bx = binary_dot(basis[0], x) | (binary_dot(basis[1], x) << 1)
                        y = bx ^ (a if incoming else 0)
                        updated = incoming ^ binary_dot(c, x)
                        weights[updated != 0] += z**y.bit_count()
                    weights = [value/(len(bases)*16) for value in weights]
                    empty, active = packet.local_matrices(z, 2, 1, q=2,
                                                         input_weight=x.bit_count() if identity else None)
                    self.assertEqual(weights, (active if x else empty)[incoming])

    def test_empty_packets_preserve_nonzero_state(self):
        for width, state in itertools.product((1, 2, 4, 16), (1, 2, 4)):
            empty, active = packet.local_matrices(Fraction(1, 2), width, state)
            self.assertEqual(empty[0], [1, 0])
            self.assertEqual(empty[1][0], 0)
            self.assertEqual(sum(empty[1]), Fraction(17, 32)**width)
            self.assertEqual(sum(active[1]), sum(empty[1]))

    def test_two_dimensional_toeplitz_state(self):
        z = Fraction(2, 3)
        for x, incoming in itertools.product(range(4), repeat=2):
            weights = [Fraction(0), Fraction(0)]
            for a, c in itertools.product(range(8), repeat=2):
                y = x ^ binary_toeplitz(2, 2, a, incoming)
                updated = incoming ^ binary_toeplitz(2, 2, c, x)
                weights[updated != 0] += z**y.bit_count()/64
            empty, active = packet.local_matrices(z, 2, 2, q=2, input_weight=x.bit_count())
            self.assertEqual(weights, (active if x else empty)[incoming != 0])

    def test_positive_aggregate_equals_weight_sum(self):
        for width, state, identity in itertools.product((2, 4, 16), (1, 3), (False, True)):
            x, z = .021, .731
            expected = sum(math.comb(width, w)*(15*x)**w*(1+x)**(width-w)*matrix
                           for w, matrix in enumerate(packet.local_by_weight(z, width, state, identity)))
            np.testing.assert_allclose(packet.aggregate(x, z, width, state, identity), expected, rtol=2e-14)

    def test_all_words_at_z_one(self):
        for active, identity in itertools.product((1, 2, 3), (False, True)):
            regions = packet.log_regions(1., active, rows=4, width=2, state=1, identity=identity)
            actual = source_moment(regions, active, 4, 2, exact_systematic=True)
            self.assertAlmostEqual(actual, active*math.log(16**2-1), places=12)

    def test_small_source_measure_exact(self):
        z = Fraction(2, 3)
        for active, identity in itertools.product((1, 2), (False, True)):
            matrices = [packet.local_matrices(z, 2, 1)[0]]
            matrices += [packet.local_matrices(z, 2, 1, input_weight=w if identity else None)[1] for w in (1, 2)]
            parity = [[sum(Fraction(math.comb(active, w)*15**w, 16**active)*matrices[w][i][j]
                           for w in range(active+1)) for j in range(2)] for i in range(2)]
            expected = 15**active*sum(matrices[active][0][i]*parity[i][j] for i in range(2) for j in range(2))
            regions = packet.log_regions(float(z), active, rows=2, width=2, state=1, identity=identity)
            actual = source_moment(regions, active, 2, 1, exact_systematic=True)
            self.assertAlmostEqual(actual, math.log(float(expected)), places=12)

    def test_positive_pascal_matches_original_source_count(self):
        for active, identity in itertools.product((1, 2, 4, 8), (False, True)):
            regions = packet.log_regions(.82, active, rows=16, width=4, state=2, identity=identity)
            old = source_moment(regions, active, 6, 4, exact_systematic=True)
            new = packet.positive_source_moment(regions, active, 6, 4)
            self.assertAlmostEqual(new, old, places=11)

    def test_recorded_floating_midweight_witness(self):
        margin = packet.sparse_bound(128, .9703716207340908, length=64, dimension=40,
                                     width=16, state=5, identity=True)
        self.assertAlmostEqual(margin, 208.4352424136817, places=6)
        # A regression of a floating calculation is not a directed bound.


if __name__ == '__main__':
    unittest.main()
