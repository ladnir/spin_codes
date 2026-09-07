"""Compare positive contractions against the independent log evaluator."""
import unittest

import numpy as np

import balanced_occupation as balanced
import composition_boxes as original
import occupation_composition_positive_v1 as positive


class PositiveContractions(unittest.TestCase):
    def test_shifted_ordinary_rare_and_deterministic_laws(self):
        rng = np.random.default_rng(573)
        for scale in (1., 80., 900.):
            regions = -rng.random((47, 4, 4))*scale
            regions[:, 3, 1] = -np.inf
            for p in (0., .03, .5, .97, 1.):
                law = original.binomial_logs(19, p)
                expected = original.shifted_mixture_logs(regions, law, 27)
                actual = positive.shifted_mixture(regions, law, 27)
                np.testing.assert_array_equal(np.isfinite(actual), np.isfinite(expected))
                mask = np.isfinite(expected)
                np.testing.assert_allclose(actual[mask], expected[mask], rtol=0, atol=2e-11)

    def test_envelope_folding_and_support_fallback(self):
        rng = np.random.default_rng(671)
        p = np.array([.1, .5, .9])
        envelope = balanced.Envelope(np.array([.3, .1, .35]), np.log(p), np.log1p(-p))
        for n, scale in ((1, 1.), (19, 1.), (129, 30.), (513, 30.), (65, 900.)):
            logs = -rng.random((n, 4, 4))*scale
            logs[:, 2, 1] = -np.inf
            logs[1:, 0, 1] = -np.inf
            expected = positive.log_fold(logs, envelope)
            actual = positive.envelope_fold(logs, envelope)
            np.testing.assert_array_equal(np.isfinite(actual), np.isfinite(expected))
            mask = np.isfinite(expected)
            np.testing.assert_allclose(actual[mask], expected[mask], rtol=0, atol=2e-10)

    def test_complete_composition_matrix(self):
        counts = {2: 10, 4: 20, 6: 10, 8: 1}
        bands = [[2], [4, 6], [8]]
        probabilities = [.2, .5, .9]
        baseline = positive.lazy.CompositionBoxes(counts, 8, bands, probabilities, 33)
        fast = positive.CompositionBoxes(counts, 8, bands, probabilities, 33)
        rng = np.random.default_rng(922)
        for scale in (5., 900.):
            regions = -rng.random((34, 4, 4))*scale
            regions[:, 2, 1] = -np.inf
            for lower in ([0, 0, 0], [5, 7, 1], [15, 8, 10]):
                expected = baseline.matrix(regions, lower, 33)
                actual = fast.matrix(regions, lower, 33)
                np.testing.assert_array_equal(np.isfinite(actual), np.isfinite(expected))
                mask = np.isfinite(expected)
                np.testing.assert_allclose(actual[mask], expected[mask], rtol=0, atol=2e-10)


if __name__ == '__main__':
    unittest.main()
