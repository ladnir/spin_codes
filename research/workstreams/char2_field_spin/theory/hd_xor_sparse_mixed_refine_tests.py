"""Scope, normalization and bounded-search tests for the sparse mixed gate."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from hd_xor_sparse_mixed_refine import (MAX_EVALUATIONS, MU_BINS, evaluate_at,
    fit_cell, geometry_count, inherited_literal, metadata, source_counts)


class SparseMixedRefineTests(unittest.TestCase):
    def test_geometry_partition_counts_each_fixed_assignment_once(self):
        self.assertEqual([mu for lo, hi in MU_BINS for mu in range(lo, hi+1)], list(range(1, 17)))
        self.assertEqual(sum(geometry_count(*bounds) for bounds in MU_BINS), 28**16)
        self.assertEqual(geometry_count(16, 16), 28)
        self.assertEqual(geometry_count(15, 15), 12096)
        with self.assertRaises(ValueError):
            geometry_count(1, 16)

    def test_source_scope_and_inherited_witnesses(self):
        self.assertEqual(len(source_counts()), 81)
        for source in source_counts():
            for bounds in MU_BINS:
                z, p = inherited_literal(source, bounds)
                self.assertTrue(0 < z < 1 and p <= 0)

    def test_floor_and_projective_normalizations(self):
        spec = metadata()
        self.assertEqual(spec['sparse_support_floor'], 17*33)
        self.assertEqual(spec['sparse_support_normalizer_exponent'], 16*33)
        self.assertEqual(spec['projective_normalizations'], 1)
        self.assertEqual(spec['endpoint_potential_factors'], 1)
        self.assertFalse(spec['directed'])
        self.assertFalse(spec['whole_code_certificate'])

    def test_bounded_search_never_exceeds_budget(self):
        def fake(source, bounds, z, p):
            return dict(mark=str(z), log2_potential=str(p),
                failure_log2_screen=((float(z)-.994)/.001)**2+(float(p)+65)**2)
        with patch('hd_xor_sparse_mixed_refine.evaluate_at', side_effect=fake):
            result = fit_cell()
        self.assertLessEqual(result['evaluations'], MAX_EVALUATIONS)
        self.assertLessEqual(result['best']['failure_log2_screen'], result['proposals'][0]['failure_log2_screen'])

    def test_original_corner_reproduced(self):
        z, p = inherited_literal((8, 0), (9, 14))
        row = evaluate_at((8, 0), (9, 14), z, p)
        self.assertAlmostEqual(row['failure_log2_screen'], -7.200543368064245, places=8)
        self.assertEqual(row['cohort_colorings'], 17)
        self.assertEqual(row['one_row_choices'], 8)
        self.assertEqual(row['two_row_geometries'], 3439311631744464)


if __name__ == '__main__':
    unittest.main()
