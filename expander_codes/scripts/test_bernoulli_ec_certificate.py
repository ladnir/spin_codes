import math
import unittest

from flint import arb

from bernoulli_ec_certificate import (
    bernoulli_wrapping_matrix_arb,
    dense_bound_arb,
    q_from_density_arb,
)
from ea_certificate import hamming_ball_bound_arb
from exact_row_ec_certificate import wrapping_enumerator_matrix_arb


class BernoulliECCertificateTests(unittest.TestCase):
    def test_hamming_ball_bound_covers_exact_small_balls(self) -> None:
        for length in (8, 15, 24):
            for cutoff in range((length - 1) // 2 + 1):
                exact = sum(math.comb(length, weight) for weight in range(cutoff + 1))
                self.assertGreaterEqual(hamming_ball_bound_arb(length, cutoff), exact)

    def test_scaled_bivariate_matrix_equals_direct_matrix(self) -> None:
        q = arb("0.137")
        z = arb("0.73")
        memory = 4
        direct = bernoulli_wrapping_matrix_arb(q, z, memory)
        bivariate = wrapping_enumerator_matrix_arb(q / (1 - q), z, memory) * (1 - q)
        for row in range(memory + 1):
            for column in range(memory + 1):
                self.assertTrue(direct[row, column].overlaps(bivariate[row, column]))

    def test_activation_probability_is_monotone(self) -> None:
        density = arb("0.0001")
        values = [q_from_density_arb(density, weight) for weight in range(1, 8)]
        self.assertTrue(all(left < right for left, right in zip(values, values[1:])))

    def test_uniform_input_rows_have_common_sum(self) -> None:
        z = arb("0.4")
        matrix = bernoulli_wrapping_matrix_arb(arb(1) / 2, z, 5)
        expected = (1 + z) / 2
        for row in range(6):
            row_sum = sum((matrix[row, column] for column in range(6)), arb(0))
            self.assertTrue(row_sum.overlaps(expected))

    def test_dense_partition_covers_weightwise_bound(self) -> None:
        k = 32
        n = 160
        cutoff = 20
        density = arb("0.1")
        start = 8
        bound, blocks = dense_bound_arb(
            k=k,
            n=n,
            cutoff=cutoff,
            density=density,
            start=start,
        )
        output_ball = sum(math.comb(n, weight) for weight in range(cutoff + 1))
        weightwise = sum(
            arb(math.comb(k, message_weight))
            * output_ball
            * (1 - q_from_density_arb(density, message_weight)) ** n
            for message_weight in range(start, k + 1)
        )
        self.assertGreaterEqual(bound, weightwise)
        self.assertGreaterEqual(blocks, 1)

    def test_dense_partition_adds_early_stops(self) -> None:
        bound, blocks = dense_bound_arb(
            k=1024,
            n=2048,
            cutoff=100,
            density=arb("0.01"),
            start=100,
        )
        self.assertGreater(bound, 0)
        self.assertGreaterEqual(blocks, 10)


if __name__ == "__main__":
    unittest.main()
