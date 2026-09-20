#!/usr/bin/env python3

from __future__ import annotations

import itertools
import math
import unittest
from collections import Counter

from binary_biregular_diagnostic import (
    biregular_parity_shell_counts,
    biregular_parity_shell_distribution,
    central_fourier_l1_log2,
    degree_five_conditional_variances,
    degree_five_dense_scan,
    degree_three_conditional_variances,
    degree_three_dense_scan,
    degree_three_parity_shell_counts,
    exact_biregular_fixed_marker_block,
    geometric_blocks,
    parity_probability,
    saddle_biregular_logterm,
    saddle_biregular_block_logterm,
)


class BinaryBiregularDiagnosticTests(unittest.TestCase):
    def test_direct_degree_three_shell_counts_match_generic_counts(self) -> None:
        for support_size in range(16):
            self.assertEqual(
                degree_three_parity_shell_counts(
                    left_vertices=15, support_size=support_size
                ),
                biregular_parity_shell_counts(
                    left_vertices=15,
                    right_degree=3,
                    support_size=support_size,
                ),
            )

    def test_shell_distribution_matches_exhaustive_subsets(self) -> None:
        left_vertices = 6
        right_degree = 2
        support_size = 3
        counts: Counter[int] = Counter()
        for subset in itertools.combinations(range(left_vertices), support_size):
            parities = [0] * (left_vertices // right_degree)
            for slot in subset:
                parities[slot // right_degree] ^= 1
            counts[sum(parities)] += 1
        expected = {
            weight: count / math.comb(left_vertices, support_size)
            for weight, count in counts.items()
        }
        self.assertEqual(
            set(expected),
            set(biregular_parity_shell_distribution(
                left_vertices=left_vertices,
                right_degree=right_degree,
                support_size=support_size,
            )),
        )
        actual = biregular_parity_shell_distribution(
            left_vertices=left_vertices,
            right_degree=right_degree,
            support_size=support_size,
        )
        for weight, probability in expected.items():
            self.assertAlmostEqual(actual[weight], probability)

    def test_full_support_exposes_right_degree_parity(self) -> None:
        odd = biregular_parity_shell_distribution(
            left_vertices=12, right_degree=3, support_size=12
        )
        even = biregular_parity_shell_distribution(
            left_vertices=12, right_degree=2, support_size=12
        )
        self.assertEqual(odd, {4: 1.0})
        self.assertEqual(even, {0: 1.0})

    def test_marked_slot_parity_probability(self) -> None:
        x = 0.3
        q = parity_probability(3, x)
        direct = sum(
            math.comb(3, j) * x**j for j in (1, 3)
        ) / (1.0 + x) ** 3
        self.assertAlmostEqual(q, direct)

    def test_central_fourier_sum_matches_krawtchouk_definition(self) -> None:
        k = 15
        degree = 3
        region_length = k // degree
        r = (k - 1) // 2
        total = 0.0
        for groups in range(region_length + 1):
            slots = degree * groups
            krawtchouk = sum(
                (-1) ** selected
                * math.comb(slots, selected)
                * math.comb(k - slots, r - selected)
                for selected in range(max(0, r - (k - slots)), min(slots, r) + 1)
            )
            total += (
                math.comb(region_length, groups)
                * abs(krawtchouk)
                / math.comb(k, r)
            )
        self.assertAlmostEqual(
            central_fourier_l1_log2(left_vertices=k, right_degree=degree),
            math.log2(total),
        )

    def test_degree_five_conditional_variances_at_half(self) -> None:
        import numpy as np

        even, odd = degree_five_conditional_variances(np.array([0.5]))
        self.assertAlmostEqual(float(even[0]), 5.0 / 16.0)
        self.assertAlmostEqual(float(odd[0]), 5.0 / 16.0)

    def test_degree_five_dense_scan_is_finite(self) -> None:
        result = degree_five_dense_scan(
            k=25,
            left_degree=10,
            cutoff=4,
            support_start=5,
            support_limit=20,
        )
        self.assertTrue(math.isfinite(result.worst_log2_term))
        self.assertGreaterEqual(result.worst_support, 5)
        self.assertLessEqual(result.worst_support, 20)

    def test_degree_three_conditional_variances_at_half(self) -> None:
        import numpy as np

        even, odd = degree_three_conditional_variances(np.array([0.5]))
        self.assertAlmostEqual(float(even[0]), 3.0 / 16.0)
        self.assertAlmostEqual(float(odd[0]), 3.0 / 16.0)

    def test_degree_three_dense_scan_is_finite(self) -> None:
        result = degree_three_dense_scan(
            k=15,
            left_degree=6,
            cutoff=4,
            support_start=6,
            support_limit=9,
        )
        self.assertTrue(math.isfinite(result.worst_log2_term))
        self.assertGreaterEqual(result.worst_support, 6)
        self.assertLessEqual(result.worst_support, 9)

    def test_saddle_is_finite(self) -> None:
        term = saddle_biregular_logterm(
            code="ec",
            k=18,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            message_weight=5,
        )
        self.assertTrue(math.isfinite(term.log2_bound))

    def test_exact_fixed_marker_block_matches_single_support(self) -> None:
        block = exact_biregular_fixed_marker_block(
            code="ec",
            k=15,
            left_degree=10,
            right_degree=5,
            cutoff=4,
            memory=3,
            support_start=2,
            support_limit=2,
            output_marker=0.7,
        )
        self.assertEqual(block.worst_support, 2)
        self.assertAlmostEqual(block.worst_log2_term, block.summed_log2_bound)
        self.assertTrue(math.isfinite(block.log2_terms[0]))

    def test_block_saddle_is_finite(self) -> None:
        term = saddle_biregular_block_logterm(
            code="ec",
            k=18,
            left_degree=6,
            right_degree=3,
            cutoff=4,
            memory=3,
            support_start=4,
            support_limit=6,
        )
        self.assertTrue(math.isfinite(term.log2_bound))

    def test_geometric_blocks_are_contiguous(self) -> None:
        blocks = geometric_blocks(7, 100, 0.1)
        self.assertEqual(blocks[0][0], 7)
        self.assertEqual(blocks[-1][1], 100)
        for left, right in zip(blocks, blocks[1:]):
            self.assertEqual(left[1] + 1, right[0])


if __name__ == "__main__":
    unittest.main()
