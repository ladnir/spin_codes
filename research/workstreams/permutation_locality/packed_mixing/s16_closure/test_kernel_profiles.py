"""Exhaustive checks of the independent expansion-profile partition."""
from fractions import Fraction as Q
from math import comb
from types import SimpleNamespace
import unittest

import numpy as np
from flint import arb, ctx
import kernel_maps
import kernel_profiles as kernel
import adapters_regional
import scalar_cover as sc


def endpoint(value):
    m, e = value.upper().man_exp()
    return Q(int(m))*Q(2)**int(e)


class ProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.images = tuple(((3 if s&1 else 0) ^ (17 if s&2 else 0) ^ (68 if s&4 else 0)) for s in range(8))
        cls.columns = [1, 2, 4, 3, 5, 7, 6, 1]
        cls.base = kernel_maps.prepare_maps(cls.images, cls.columns, 3)
        cls.data = kernel.attach(cls.base)

    def setUp(self):
        ctx.prec = 192

    def exact(self, occupancy, z):
        result = [[Q(0)]*8 for _ in range(8)]
        for state in range(8):
            for word in range(256):
                if bool(word&15)+bool(word>>4) != occupancy:
                    continue
                syndrome = 0
                for bit, column in enumerate(self.columns):
                    if word >> bit & 1:
                        syndrome ^= column
                mass = z**(word^self.images[state]).bit_count()/(comb(2, occupancy)*15**occupancy)
                if state:
                    for target in range(1, 8):
                        result[state][target^syndrome] += mass/7
                else:
                    result[state][syndrome] += mass
        return result

    def test_partition_is_separate_and_complete(self):
        np.testing.assert_array_equal(self.base['birth_class_levels'], self.data['birth_class_levels'])
        self.assertGreater(self.data['profile_count'], len(self.data['birth_class_levels']))
        self.assertEqual(self.data['map_sha256'], self.base['map_sha256'])
        self.assertEqual(self.data['profile_state_indices'][0], -1)

    def test_every_state_and_fixed_occupancy_composition(self):
        data = self.data; local = kernel.outward_local_at_z(data, arb(3)/4, Q(1, 4))
        n = local[0].nrows(); exact = [self.exact(j, Q(3, 4)) for j in range(3)]
        for occupancy, matrix in enumerate(local):
            self.assertGreaterEqual(endpoint(matrix[0, 0]), exact[occupancy][0][0])
            for profile in range(data['profile_count']):
                residual = sum(max(Q(0), exact[occupancy][0][s]-endpoint(matrix[0, 2])/7)
                    for s in range(1, 8) if data['profile_state_indices'][s] == profile)
                self.assertGreaterEqual(endpoint(matrix[0, 3+profile]), residual)
            for state in range(1, 8):
                source = 3+data['profile_state_indices'][state]
                self.assertGreaterEqual(endpoint(matrix[source, 0]), exact[occupancy][state][0])
                for target in range(1, 8):
                    self.assertGreaterEqual(endpoint(matrix[source, 2])/7, exact[occupancy][state][target])
        values = [Q(1)]*8; bounded = arb(1)*local[0]
        sequence = (0, 1, 2, 1, 0, 2)
        bounded = None
        for occupancy in sequence:
            values = [sum(a*b for a, b in zip(row, values)) for row in exact[occupancy]]
            bounded = local[occupancy] if bounded is None else local[occupancy]*bounded
            for state in range(8):
                source = 0 if not state else 3+data['profile_state_indices'][state]
                self.assertGreaterEqual(endpoint(sum((bounded[source, j] for j in range(n)), arb(0))), values[state])

    def test_iid_float_and_outward_and_adapter(self):
        probabilities = sc.probabilities(Q(1, 5))
        proposed = kernel.floating(self.data, list(map(float, probabilities)), .25)
        exact = kernel.outward(self.data, probabilities, Q(1, 4))
        n = exact.nrows()
        np.testing.assert_allclose(proposed, [[float(exact[i, j]) for j in range(n)] for i in range(n)],
            rtol=2e-10, atol=2e-12)
        local = adapters_regional.local_operators(SimpleNamespace(data=self.data),
            dict(parameters=['1/4'], s16_birth_density_activity='1/5'))
        self.assertEqual(local[0].nrows(), 3+self.data['profile_count'])


if __name__ == '__main__':
    unittest.main()
