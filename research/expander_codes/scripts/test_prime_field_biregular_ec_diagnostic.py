import math
import unittest
from itertools import product

import numpy as np

from prime_field_biregular_ec_diagnostic import (
    balanced_constraint_slot_trace_matrix,
    balanced_singleton_constraint_slot_trace_matrix,
    balanced_slot_trace_matrix,
    balanced_point_mass_logterm,
    balanced_ec_block_logbound,
    candidate_parameters,
    exact_balanced_ec_prefix_logbound,
    constraint_balanced_ec_block_logbound,
    constraint_exact_balanced_ec_band_logbound,
    constraint_full_support_ec_logbound,
    constraint_uniform_occupancy_trace_transfers,
    nonwrapping_constraint_trace_matrices,
    nonwrapping_singleton_constraint_trace_matrices,
    nonwrapping_trace_matrices,
    balanced_occupancy_prefix,
    tail_subtracted_region_matrix,
    singleton_constraint_uniform_support_trace_transfers,
    singleton_constraint_uniform_support_trace_transfers_scaled,
    singleton_constraint_exact_balanced_ec_band_at_markers,
    singleton_constraint_exact_balanced_ec_band_at_markers_scaled,
    singleton_constraint_exact_point_at_markers_scaled,
    singleton_constraint_low_equation_exact_point_logbound,
    singleton_constraint_low_equation_log_point_logbound,
    _truncated_bivariate_matrix_product,
    uniform_occupancy_trace_transfers,
)


class PrimeFieldBiregularECDiagnosticTests(unittest.TestCase):
    def test_constraint_trace_counts_active_zero_equation(self) -> None:
        empty, occupied = nonwrapping_constraint_trace_matrices(
            memory=2, output_marker=0.4, equation_marker=0.2
        )
        np.testing.assert_allclose(empty[0], [0.4, 0.2, 0.0])
        np.testing.assert_allclose(empty[1], [0.4, 0.0, 0.2])
        np.testing.assert_allclose(empty[:2], occupied[:2])
        np.testing.assert_allclose(empty[2], [0.0, 0.0, 1.0])
        np.testing.assert_allclose(occupied[2], [0.4, 0.0, 0.2])

    def test_constraint_group_saddle_factor(self) -> None:
        empty, occupied = nonwrapping_constraint_trace_matrices(
            memory=1, output_marker=0.6, equation_marker=0.2
        )
        x = 0.15
        matrix = balanced_constraint_slot_trace_matrix(
            memory=1,
            group_size=4,
            slot_marker=x,
            output_marker=0.6,
            equation_marker=0.2,
        )
        np.testing.assert_allclose(
            matrix, empty + (math.pow(1.0 + x, 4) - 1.0) * occupied
        )

    def test_singleton_group_cannot_cancel_from_zero_state(self) -> None:
        empty, singleton, collision = (
            nonwrapping_singleton_constraint_trace_matrices(
                memory=2, output_marker=0.4, equation_marker=0.2
            )
        )
        np.testing.assert_allclose(empty[2], [0.0, 0.0, 1.0])
        np.testing.assert_allclose(singleton[2], [0.4, 0.0, 0.0])
        np.testing.assert_allclose(collision[2], [0.4, 0.0, 0.2])

    def test_singleton_saddle_separates_linear_group_term(self) -> None:
        empty, singleton, collision = (
            nonwrapping_singleton_constraint_trace_matrices(
                memory=1, output_marker=0.6, equation_marker=0.2
            )
        )
        x = 0.15
        matrix = balanced_singleton_constraint_slot_trace_matrix(
            memory=1,
            group_size=4,
            slot_marker=x,
            output_marker=0.6,
            equation_marker=0.2,
        )
        np.testing.assert_allclose(
            matrix,
            empty
            + 4 * x * singleton
            + (math.pow(1.0 + x, 4) - 1.0 - 4 * x) * collision,
        )

    def test_singleton_exact_support_slices_match_two_groups(self) -> None:
        empty, singleton, collision = (
            nonwrapping_singleton_constraint_trace_matrices(
                memory=1, output_marker=0.6, equation_marker=0.2
            )
        )
        slices = singleton_constraint_uniform_support_trace_transfers(
            region_length=2,
            group_size=2,
            max_weight=2,
            memory=1,
            output_marker=0.6,
            equation_marker=0.2,
        )
        expected_weight_two = (
            empty @ collision
            + 4 * (singleton @ singleton)
            + collision @ empty
        ) / math.comb(4, 2)
        np.testing.assert_allclose(slices[2], expected_weight_two)

    def test_scaled_singleton_slices_match_unscaled_calculation(self) -> None:
        ordinary = singleton_constraint_uniform_support_trace_transfers(
            region_length=5,
            group_size=3,
            max_weight=6,
            memory=2,
            output_marker=0.6,
            equation_marker=0.2,
        )
        scaled, log_scales = (
            singleton_constraint_uniform_support_trace_transfers_scaled(
                region_length=5,
                group_size=3,
                max_weight=6,
                memory=2,
                output_marker=0.6,
                equation_marker=0.2,
            )
        )
        for expected, matrix, log_scale in zip(ordinary, scaled, log_scales):
            np.testing.assert_allclose(
                math.exp(float(log_scale)) * matrix,
                expected,
                rtol=2e-12,
                atol=2e-14,
            )

    def test_scaled_singleton_band_matches_unscaled_calculation(self) -> None:
        parameters = {
            "prime": 101,
            "k": 6,
            "n": 12,
            "cutoff": 4,
            "region_count": 2,
            "memory": 2,
            "support_start": 1,
            "support_limit": 3,
            "structural_markers": (0.6, 0.2),
            "field_markers": (0.6, 0.2),
        }
        ordinary = singleton_constraint_exact_balanced_ec_band_at_markers(
            **parameters
        )
        scaled = singleton_constraint_exact_balanced_ec_band_at_markers_scaled(
            **parameters
        )
        self.assertAlmostEqual(scaled.structural_log2, ordinary.structural_log2)
        self.assertAlmostEqual(scaled.field_log2, ordinary.field_log2)
        self.assertAlmostEqual(scaled.total_log2, ordinary.total_log2)

    def test_truncated_bivariate_product_matches_direct_convolution(self) -> None:
        rng = np.random.default_rng(17)
        left = rng.random((3, 2, 3, 3))
        right = rng.random((2, 3, 3, 3))
        actual = _truncated_bivariate_matrix_product(
            left, right, max_x_degree=2, max_equation_degree=2
        )
        expected = np.zeros_like(actual)
        for x_left in range(left.shape[0]):
            for a_left in range(left.shape[1]):
                for x_right in range(right.shape[0]):
                    for a_right in range(right.shape[1]):
                        x = x_left + x_right
                        a = a_left + a_right
                        if x <= 2 and a <= 2:
                            expected[x, a] += (
                                left[x_left, a_left]
                                @ right[x_right, a_right]
                            )
        np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-12)

    def test_low_equation_exact_point_is_below_structural_chernoff(self) -> None:
        parameters = {
            "prime": 127,
            "k": 12,
            "n": 24,
            "cutoff": 4,
            "region_count": 4,
            "memory": 2,
            "message_weight": 3,
        }
        z = 0.19
        exact = singleton_constraint_low_equation_exact_point_logbound(
            **parameters, output_marker=z
        )
        chernoff = singleton_constraint_exact_point_at_markers_scaled(
            **parameters, output_marker=z, equation_marker=0.31
        )
        self.assertLessEqual(exact, chernoff + 2e-9)
        log_exact = singleton_constraint_low_equation_log_point_logbound(
            **parameters, output_marker=z
        )
        self.assertAlmostEqual(log_exact, exact, places=10)

    def test_binary_active_state_matches_uniform_output(self) -> None:
        empty, occupied = nonwrapping_trace_matrices(
            prime=2, memory=3, output_marker=0.4, equation_marker=0.2
        )
        for state in range(3):
            self.assertAlmostEqual(empty[state, 0], 0.2)
            self.assertAlmostEqual(empty[state, state + 1], 0.5)
            np.testing.assert_allclose(empty[state], occupied[state])

    def test_zero_state_distinguishes_empty_and_occupied(self) -> None:
        empty, occupied = nonwrapping_trace_matrices(
            prime=5, memory=2, output_marker=0.3, equation_marker=0.1
        )
        np.testing.assert_allclose(empty[2], [0.0, 0.0, 1.0])
        np.testing.assert_allclose(occupied[2], [0.3, 0.0, 0.1])

    def test_group_saddle_factor(self) -> None:
        empty, occupied = nonwrapping_trace_matrices(
            prime=7, memory=1, output_marker=0.6, equation_marker=0.2
        )
        x = 0.15
        matrix = balanced_slot_trace_matrix(
            prime=7,
            memory=1,
            group_size=4,
            slot_marker=x,
            output_marker=0.6,
            equation_marker=0.2,
        )
        np.testing.assert_allclose(
            matrix, empty + (math.pow(1.0 + x, 4) - 1.0) * occupied
        )

    def test_rate_half_candidate_dimensions(self) -> None:
        k, n, ell = candidate_parameters(2**21, 30)
        self.assertEqual((k, n, ell), (1_048_575, 2_097_150, 69_905))
        self.assertEqual(k // ell, 15)

    def test_uniform_slice_transfer_at_length_one(self) -> None:
        empty, occupied = nonwrapping_trace_matrices(
            prime=5, memory=1, output_marker=0.3, equation_marker=0.1
        )
        slices = uniform_occupancy_trace_transfers(
            length=1,
            max_occupied=1,
            prime=5,
            memory=1,
            output_marker=0.3,
            equation_marker=0.1,
        )
        np.testing.assert_allclose(slices[0], empty)
        np.testing.assert_allclose(slices[1], occupied)

    def test_constraint_uniform_slice_transfer_at_length_one(self) -> None:
        empty, occupied = nonwrapping_constraint_trace_matrices(
            memory=1, output_marker=0.3, equation_marker=0.1
        )
        slices = constraint_uniform_occupancy_trace_transfers(
            length=1,
            max_occupied=1,
            memory=1,
            output_marker=0.3,
            equation_marker=0.1,
        )
        np.testing.assert_allclose(slices[0], empty)
        np.testing.assert_allclose(slices[1], occupied)

    def test_constraint_bounds_are_finite(self) -> None:
        exact = constraint_exact_balanced_ec_band_logbound(
            prime=17,
            k=30,
            n=60,
            cutoff=17,
            region_count=6,
            memory=1,
            support_start=1,
            support_limit=3,
        )
        block = constraint_balanced_ec_block_logbound(
            prime=17,
            k=30,
            n=60,
            cutoff=17,
            region_count=6,
            memory=1,
            support_start=4,
            support_limit=7,
        )
        full = constraint_full_support_ec_logbound(
            prime=17,
            n=60,
            cutoff=17,
            memory=1,
            message_weight=30,
        )
        for bound in (exact, block, full):
            self.assertTrue(math.isfinite(bound.total_log2))
            self.assertLessEqual(max(bound.structural_markers[:2]), 1.0)
            self.assertLessEqual(max(bound.field_markers[:2]), 1.0)

    def test_constraint_trace_dominates_tiny_exhaustive_encoder(self) -> None:
        prime = 3
        multiplicities = (0, 2, 1, 2)
        occupied_flags = tuple(value > 0 for value in multiplicities)
        outcomes = 0
        low_weight = [0] * (len(multiplicities) + 1)
        label_count = sum(multiplicities)
        for labels in product(range(1, prime), repeat=label_count):
            inputs = []
            cursor = 0
            for count in multiplicities:
                inputs.append(sum(labels[cursor:cursor + count]) % prime)
                cursor += count
            for feedback in product(range(prime), repeat=len(inputs)):
                previous = 0
                weight = 0
                for value, alpha in zip(inputs, feedback):
                    output = (value + alpha * previous) % prime
                    weight += output != 0
                    previous = output
                outcomes += 1
                for cutoff in range(weight, len(low_weight)):
                    low_weight[cutoff] += 1

        z = 0.43
        empty, occupied = nonwrapping_constraint_trace_matrices(
            memory=1,
            output_marker=z,
            equation_marker=1.0 / (prime - 1),
        )
        transfer = np.eye(2)
        for is_occupied in occupied_flags:
            transfer = transfer @ (occupied if is_occupied else empty)
        trace = float(np.sum(transfer[1]))
        for cutoff, count in enumerate(low_weight):
            probability = count / outcomes
            self.assertLessEqual(probability, trace * z ** (-cutoff) + 1e-14)

    def test_point_mass_bound_is_finite(self) -> None:
        value = balanced_point_mass_logterm(
            prime=17,
            k=30,
            n=60,
            cutoff=17,
            region_count=6,
            message_weight=12,
        )
        self.assertTrue(math.isfinite(value))

    def test_joint_exact_prefix_bound_is_finite(self) -> None:
        bound = exact_balanced_ec_prefix_logbound(
            prime=17,
            k=30,
            n=60,
            cutoff=17,
            region_count=6,
            memory=1,
            support_limit=3,
        )
        self.assertEqual(bound.support_start, 1)
        self.assertEqual(bound.support_limit, 3)
        self.assertGreater(bound.output_marker, 0.0)
        self.assertLessEqual(bound.output_marker, 1.0)
        self.assertGreater(bound.equation_marker, 0.0)
        self.assertLessEqual(bound.equation_marker, 1.0)
        self.assertTrue(math.isfinite(bound.total_log2))

    def test_joint_exact_support_band_is_finite(self) -> None:
        bound = exact_balanced_ec_prefix_logbound(
            prime=17,
            k=30,
            n=60,
            cutoff=17,
            region_count=6,
            memory=1,
            support_start=2,
            support_limit=4,
        )
        self.assertEqual((bound.support_start, bound.support_limit), (2, 4))
        self.assertTrue(math.isfinite(bound.total_log2))

    def test_saddle_support_block_is_finite(self) -> None:
        bound, markers = balanced_ec_block_logbound(
            prime=17,
            k=30,
            n=60,
            cutoff=17,
            region_count=6,
            memory=1,
            support_start=4,
            support_limit=7,
        )
        self.assertEqual((bound.support_start, bound.support_limit), (4, 7))
        self.assertEqual(set(markers), {"structural", "field"})
        self.assertTrue(math.isfinite(bound.total_log2))

    def test_occupancy_prefix_matches_single_weight_helper(self) -> None:
        prefix = balanced_occupancy_prefix(
            region_length=5, group_size=2, limit=4
        )
        from prime_field_regular_ea_trace_diagnostic import (
            balanced_occupancy_size_distribution,
        )
        for draws, law in enumerate(prefix):
            np.testing.assert_allclose(
                law,
                balanced_occupancy_size_distribution(5, 2, draws),
            )

    def test_complete_tail_subtraction_is_zero(self) -> None:
        tail = tail_subtracted_region_matrix(
            prime=5,
            k=6,
            region_length=3,
            group_size=2,
            memory=1,
            slot_marker=0.4,
            output_marker=0.7,
            equation_marker=0.2,
            prefix_limit=7,
        )
        np.testing.assert_allclose(tail, np.zeros((2, 2)), atol=2e-14)


if __name__ == "__main__":
    unittest.main()
