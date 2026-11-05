from unittest import TestCase, main
import numpy as np
import hd_post_refresh_diagnostic as diagnostic


class DiagnosticTests(TestCase):
    def test_distinct_cross_decomposition_matches_complete_tiny_columns(self):
        weights = np.array([0., 2., 3., 2.])
        measure = dict(log2_probability=-2., log2_atom00=-1.,
            log2_atom01=-3., log2_atom10=-2., log2_scale01=.5,
            log2_scale10=-4., log2_scale11=-3.)
        matrix, _ = diagnostic.engine.column_matrix(weights, 2, measure, arity=2)
        for initial in range(2):
            for final in range(2):
                for degree in range(5):
                    row = diagnostic.path_breakdown(weights, 2, measure, 2, initial, final, degree)
                    self.assertAlmostEqual(row['coefficient_log2'], matrix[initial, final, degree])


if __name__ == '__main__':
    main()
