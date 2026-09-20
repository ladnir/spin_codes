#!/usr/bin/env python3

import unittest

import numpy as np

from streaming_ec_heuristic_audit import (
    audit,
    constant_slot_witness,
    modular_rank,
    striped_expander,
)


class StreamingEcHeuristicAuditTest(unittest.TestCase):
    def test_unsigned_stripes_have_constant_slot_kernel(self) -> None:
        matrix = striped_expander(
            left_degree=10,
            right_degree=5,
            region_size=7,
            seed=9,
            signed=False,
        )
        witness = constant_slot_witness(
            right_degree=5,
            region_size=7,
            prime=127,
        )
        self.assertTrue(np.all((witness @ matrix) % 127 == 0))
        self.assertLess(modular_rank(matrix, 127), matrix.shape[0])

    def test_signed_and_permutation_models_pass_reduced_rank_gate(self) -> None:
        rows = audit(
            left_degree=10,
            right_degree=5,
            region_sizes=[3, 7],
            prime=127,
            seeds=8,
        )
        for row in rows:
            self.assertEqual(row.signed_min_rank, row.message_size)
            self.assertEqual(row.signed_rank_failures, 0)
            self.assertEqual(row.permutation_min_rank, row.message_size)
            self.assertEqual(row.permutation_rank_failures, 0)


if __name__ == "__main__":
    unittest.main()
