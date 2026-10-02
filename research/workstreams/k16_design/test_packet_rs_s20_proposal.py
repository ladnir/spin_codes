"""Small floating-objective tests; no full state preparation or certificates."""
from fractions import Fraction as Q
from math import comb, log, log1p
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb_mat
import packet_rs_s20_proposal as proposal


class ProposalTests(unittest.TestCase):
    def test_scaled_power_matches_direct_small_powers(self):
        matrix = np.array([[.5, .2], [.1, .3]])
        for exponent in (0, 1, 2, 7, 16, 63):
            expected = log(np.linalg.matrix_power(matrix, exponent)[0].sum())
            self.assertAlmostEqual(proposal.log_power(matrix, exponent), expected, places=12)

    def test_scaled_power_avoids_initial_and_long_power_underflow(self):
        for base, exponent in ((1e-250, 2), (1e-200, 16384), (.999, 100000000)):
            self.assertAlmostEqual(proposal.log_power(np.array([[base]]), exponent),
                                   exponent * log(base), places=6)

    def test_invalid_matrix_is_not_an_optimistic_infinite_score(self):
        for matrix in (np.array([[np.nan]]), np.array([[-1.0]]), np.array([[np.inf]]),
                       np.zeros((2, 3))):
            with self.assertRaises(ValueError):
                proposal.log_power(matrix, 2)
        with self.assertRaises(FloatingPointError):
            proposal.log_power(np.zeros((2, 2)), 2)

    def test_iid_macro_matches_scalar_binomial_identity(self):
        local = np.array([[[.5**j]] for j in range(33)])
        for probability in (.001, .1, .5, 1.0):
            self.assertAlmostEqual(proposal.iid_macro(local, probability)[0, 0],
                                   (1 - 15 * probability / 32)**32, places=13)

    def _model(self):
        with patch.object(proposal.sparse, 'validated', return_value=20), \
                patch.object(proposal.sparse, 'source_snapshot', return_value={'frozen': 'source'}):
            return proposal.ProposalModel({'map_sha256': 'toy'}, {'map': 'toy'})

    def test_score_matches_full_conditioning_formula(self):
        model = self._model()
        tilt, probability, q = Q(1, 100), Q(1, 16), 32
        local = np.array([[[.5**j]] for j in range(33)])
        with patch.object(model, 'local', return_value=local):
            score = model.score(q, tilt, probability)
        lc = log(comb(8192, q))
        expected = (lc + q * model.log_beta + .01 * (2**21 // 10)
                    - 64 * (lc + q * log(1/16) + (8192-q) * log1p(-1/16))
                    + 64 * 8192 * log1p(-15/32/16))
        self.assertAlmostEqual(score, expected, places=8)

    def test_marker_optimization_recovers_analytic_saddlepoint(self):
        model = self._model()
        local = np.array([[[.5**j]] for j in range(33)])
        q = 64
        ratio, a = q / 8192, 15/32
        expected = ratio / (1-a+a*ratio)
        with patch.object(model, 'local', return_value=local):
            score, probability = model._marker_fit(q, Q(1, 100), (), 80)
            baseline = model.score(q, Q(1, 100), Q(q, 8192))
        self.assertAlmostEqual(float(probability), expected, places=7)
        self.assertLess(score, baseline)

    def test_proposals_are_replayable_rationals_not_endpoints(self):
        model = self._model()
        local = [arb_mat([[1]]) for _ in range(33)]
        with patch.object(proposal.q1.kernel_t64, 'local_operators', return_value=local) as operators, \
                patch.object(proposal.sparse, 'source_snapshot', return_value={'frozen': 'source'}):
            record = model.propose([32, 8192], ['.002', '.001'])
        self.assertEqual(operators.call_count, 2)
        self.assertTrue(record['proposal_only'])
        self.assertFalse(record['whole_code_certificate'])
        self.assertFalse(record['has_numerical_upper_endpoints'])
        self.assertEqual(record['tilts'], ['1/1000'])
        self.assertTrue(all(0 < Q(p) <= 1 for p in record['marker_probabilities']))
        self.assertEqual(record['witnesses']['8192']['marker_probability'], '1')
        self.assertTrue(all(w['requires_outward_replay'] for w in record['witnesses'].values()))
        self.assertNotIn('union_upper', record)

    def test_endpoint_markers_and_unbounded_settings_are_rejected(self):
        model = self._model()
        for q, probability in ((32, 1), (32, 0), (32, -1), (8193, .5)):
            with self.assertRaises(ValueError):
                model.score(q, Q(1, 100), probability)
        for occupancy in ([], [2, 1], [1, 1], [True]):
            with self.assertRaises(ValueError):
                model.propose(occupancy, ['.01'])
        with self.assertRaises(ValueError):
            model.propose([32], ['.01'], tilt_refinements=9)


if __name__ == '__main__':
    unittest.main()
