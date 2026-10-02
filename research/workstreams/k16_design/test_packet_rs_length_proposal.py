"""Small exact comparisons for floating proposals; no map census or certificate."""
from fractions import Fraction as Q
from math import comb, log
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, arb_mat, fmpq_mat, ctx
import packet_rs_length_proposal as proposal


class LengthProposalTests(unittest.TestCase):
    def _model(self, bits=17):
        with patch.object(proposal.sparse, 'validated', return_value=bits), \
                patch.object(proposal.sparse, 'source_snapshot', return_value={'frozen': 'source'}):
            return proposal.ProposalModel({'map_sha256': 'toy'}, {'map': 'toy'})

    def test_placement_matches_noncommuting_rational_products(self):
        exact_ops = [fmpq_mat([[1, 1], [0, 1]]), fmpq_mat([[1, 0], [1, 1]]),
                     fmpq_mat([[2, 1], [0, 1]])]
        floats = np.array([[[float(m[i, j]) for j in range(2)] for i in range(2)]
                           for m in exact_ops])
        for epochs in (1, 2, 3, 7):
            degree = min(2 * epochs, 5)
            exact = proposal.q1.placement(exact_ops, epochs, 2, fmpq_mat,
                                         lambda x: x, maximum_groups=degree)
            actual = proposal.scaled_placement(floats, degree, epochs=epochs, windows=2)
            reconstructed = actual.matrices * np.exp(actual.log_scales)[:, None, None]
            expected = np.array([[[float(m[i, j]) for j in range(2)] for i in range(2)]
                                 for m in exact])
            np.testing.assert_allclose(reconstructed, expected, rtol=2e-13, atol=2e-13)

    def test_exact_placement_binomial_scalar_identity(self):
        # T_j=c*a^j gives R_j=c^epochs*a^j for every placement of j packets.
        a, c, epochs, degree = .75, .875, 5, 19
        local = np.array([[[c * a**j]] for j in range(33)])
        regional = proposal.scaled_placement(local, degree, epochs=epochs)
        for q in (0, 3, 11, degree):
            mixed, scale = proposal.regional_uniform(regional, q)
            actual = proposal.log_power(mixed, 64) + 64 * scale
            expected = 64 * (epochs * log(c) + q * log((1 + 15*a)/16))
            self.assertAlmostEqual(actual, expected, places=10)

    def test_long_placement_rescaling_prevents_global_underflow(self):
        local = np.full((3, 1, 1), 1e-100)
        regional = proposal.scaled_placement(local, 3, epochs=12, windows=2)
        np.testing.assert_allclose(regional.log_scales, 12 * log(1e-100), atol=2e-11)
        mixed, scale = proposal.regional_uniform(regional, 3)
        self.assertAlmostEqual(proposal.log_power(mixed, 64) + 64*scale,
                               64*12*log(1e-100), places=8)

    def test_binomial_mixture_preserves_rare_dominant_occupancy(self):
        # P[J=0]=16^-400 underflows as a float. The corresponding R_0 dominates.
        q = 400
        regional = proposal.ScaledPlacement(np.ones((q+1, 1, 1)),
                                            np.array([0.] + [-10000.] * q))
        mixed, scale = proposal.regional_uniform(regional, q)
        self.assertAlmostEqual(proposal.log_power(mixed, 64) + 64*scale,
                               -64*q*log(16), places=8)

    def test_16_and_17_state_score_geometry_matches_closed_form(self):
        a, c, tilt, q = .875, .9375, Q(1, 100), 7
        local = np.array([[[c*a**j]] for j in range(33)])
        for bits in (16, 17):
            model = self._model(bits)
            for K in (4096, 8192, 12288):
                with patch.object(model, 'local', return_value=local):
                    value = model.score(K, q, tilt)
                expected = (log(comb(K//128, q)) + q*model.log_beta
                    + float(tilt)*(2*K//10)
                    + 64*((K//4096)*log(c) + q*log((1+15*a)/16)))
                self.assertAlmostEqual(value, expected, places=9)

    def test_multiple_lengths_reuse_locals_and_largest_q_prefix(self):
        model = self._model()
        family = [arb_mat([[1]]) for _ in range(33)]
        with patch.object(proposal.q1.kernel_t64, 'local_operators', return_value=family) as local, \
                patch.object(proposal.sparse, 'source_snapshot', return_value={'frozen': 'source'}):
            record = model.propose([4096, 8192], {4096: [3, 8], 8192: [5, 11]}, ['.02', '.01'])
            self.assertEqual(local.call_count, 2)
            self.assertEqual(model.placement_evaluations, 4)
            model.propose([4096, 8192], [3, 5], ['.01'])
            self.assertEqual(model.placement_evaluations, 4)
        self.assertTrue(record['proposal_only'])
        self.assertFalse(record['whole_code_certificate'])
        self.assertFalse(record['has_numerical_upper_endpoints'])
        self.assertEqual(record['state_bits'], 17)
        self.assertEqual(record['lengths']['4096']['tilts'], ['1/100'])
        self.assertEqual(record['lengths']['8192']['occupancy_values'], [5, 11])
        self.assertNotIn('union_upper', record)
        self.assertTrue(all(w['requires_outward_replay']
            for r in record['lengths'].values() for w in r['witnesses'].values()))

    def test_local_uses_upper_endpoint_and_restores_precision(self):
        model = self._model()
        matrix = arb_mat([[arb(1, '0.0625')]])
        previous = ctx.prec
        with patch.object(proposal.q1.kernel_t64, 'local_operators', return_value=[matrix]*33):
            result = model.local('.01')
        self.assertGreater(result[0, 0, 0], 1)
        self.assertEqual(ctx.prec, previous)

    def test_invalid_geometry_and_grids_are_rejected(self):
        model = self._model()
        for K in (0, -4096, 4097, True, 4096.0):
            with self.assertRaises(ValueError):
                proposal.geometry(K)
        for grid in ([], [3, 3], [5, 3], [2], [33], [True]):
            with self.assertRaises(ValueError):
                model.propose([4096], grid, ['.01'])
        for tilts in ([], ['0'], ['-.1'], ['.1', '1/10'], ['nan'], ['inf'], '.1'):
            with self.assertRaises(ValueError):
                model.propose([4096], [3], tilts)
        with self.assertRaises(ValueError):
            model.propose([4096], {8192: [3]}, ['.01'])
        with self.assertRaises(ValueError):
            self._model(23)

    def test_nonfinite_and_underflow_are_errors_not_infinite_margins(self):
        for local in (np.array([[[np.nan]]]), np.array([[[-1.]]]), np.zeros((1, 1, 1))):
            with self.assertRaises((ValueError, FloatingPointError)):
                proposal.scaled_placement(local, 0, epochs=1, windows=1)
        model = self._model()
        tiny = arb_mat([[arb('1e-1000')]])
        with patch.object(proposal.q1.kernel_t64, 'local_operators', return_value=[tiny]*33):
            with self.assertRaises(FloatingPointError):
                model.local('.01')
        regional = proposal.ScaledPlacement(np.ones((2, 2, 2)), np.array([0., 0.]))
        bad_matrices = regional.matrices.copy()
        bad_matrices[:, 0, 1] = np.nextafter(0., 1.)
        # Scaling a positive subnormal down must not silently erase this entry.
        bad_matrices[:, 0, 0] = 2.
        with self.assertRaises(FloatingPointError):
            proposal.regional_uniform(proposal.ScaledPlacement(bad_matrices, regional.log_scales), 1)

    def test_changed_source_refuses_proposal_receipt(self):
        model = self._model()
        with patch.object(model, 'local', return_value=np.ones((33, 1, 1))), \
                patch.object(proposal.sparse, 'source_snapshot', return_value={'changed': 'source'}):
            with self.assertRaises(RuntimeError):
                model.propose([4096], [3], ['.01'])


if __name__ == '__main__':
    unittest.main()
