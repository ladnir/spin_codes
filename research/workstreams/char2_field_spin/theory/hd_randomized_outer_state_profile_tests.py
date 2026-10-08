"""Small exact checks of the flat-density transfer and its log implementation."""

from fractions import Fraction as F
from itertools import combinations, product
import unittest

from hd_blockwise_shared_transfer import multiply
from hd_randomized_outer_state_profile import (
    local_profiles, punctured_state_moment, one_row_receipt, route_bound_one_row,
)


def feedback(vector):
    first = second = 0
    for node, value in enumerate(vector):
        first ^= value
        second ^= multiply(node, value, 4)
    return first, second


def emission(state):
    return tuple(state[0] ^ multiply(node, state[1], 4) for node in range(4))


class WeightedStateProfileTests(unittest.TestCase):
    def test_nonuniform_column_caps_and_zero_atom_removal(self):
        # Toy maps: F=I, A=(1,node), C=A^T over GF4. CA=0 and A is MDS.
        tilt = F(1, 3)
        single = (1+3*tilt)/4
        states = tuple(product(range(4), repeat=2))
        for width in range(5):
            for support in combinations(range(4), width):
                vectors = []
                for values in product(range(4), repeat=width):
                    vector = [0]*4
                    for index, value in zip(support, values):
                        vector[index] = value
                    vectors.append(tuple(vector))
                for target in states:
                    for scalar in range(1, 4):
                        weighted = F(0)
                        for vector in vectors:
                            syndrome = feedback(vector)
                            state = tuple(multiply(scalar, x ^ y, 4)
                                          for x, y in zip(target, syndrome))
                            if any(state):
                                output = tuple(x ^ y for x, y in zip(vector, emission(state)))
                                weighted += tilt**sum(x != 0 for x in output)
                        weighted /= 4**width
                        cap = single**width
                        if target == (0, 0):
                            cap -= F(1, 4**width)
                        self.assertLessEqual(weighted, cap)

    def test_perturbed_feedthrough_is_invertible(self):
        for scalar in range(1, 4):
            images = set()
            for vector in product(range(4), repeat=4):
                state = tuple(multiply(scalar, x, 4) for x in feedback(vector))
                images.add(tuple(x ^ y for x, y in zip(vector, emission(state))))
            self.assertEqual(len(images), 4**4)

    def test_uniform_state_enumerator(self):
        from math import log2
        for length in range(5):
            exact = sum(F(1, 3)**sum(x != 0 for x in emission(state)[:length])
                        for state in product(range(4), repeat=2) if any(state))/15
            actual = punctured_state_moment(length, log2(1/3), 4, 2)
            self.assertAlmostEqual(actual, log2(exact), places=12)

    def test_empty_packets_preserve_zero_and_contract_residual(self):
        from math import isinf, log2
        log_z = -.25
        packet = local_profiles(log_z)[0]
        self.assertEqual(packet[0][0], 0)
        self.assertTrue(isinf(packet[0][1]) and packet[0][1] < 0)
        self.assertTrue(isinf(packet[0][2]) and packet[0][2] < 0)
        self.assertTrue(isinf(packet[1][0]) and packet[1][0] < 0)
        self.assertTrue(isinf(packet[2][0]) and packet[2][0] < 0)
        self.assertAlmostEqual(packet[2][2], -16+61*log_z)
        self.assertAlmostEqual(packet[2][1], log2(1-1/65536)+61*log_z)
        self.assertTrue(isinf(local_profiles(log_z, 'uniform')[0][2][2]))

    def test_exact_one_row_region_and_directed_receipt(self):
        from flint import arb
        receipt = one_row_receipt()
        self.assertTrue(arb(receipt['failure_log2_upper']) < -391)
        one = route_bound_one_row(-1/512, 1, steps=1)
        two = route_bound_one_row(-1/512, 1, steps=2)
        self.assertAlmostEqual(one, two, places=7)
        self.assertAlmostEqual(one, float(arb(receipt['failure_log2_upper'])), places=7)


if __name__ == '__main__':
    unittest.main()
