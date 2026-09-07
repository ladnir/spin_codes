"""Endpoint, orientation, and finite-field regressions for candidate screening."""
import math
import unittest
import numpy as np
import inner_candidate_screen as screen


class CandidateTests(unittest.TestCase):
    def test_binomial_endpoints_and_mass(self):
        law = np.exp(screen.binomial_logs(4, [0., .3, 1.]))
        np.testing.assert_array_equal(law[0], [1,0,0,0,0])
        np.testing.assert_array_equal(law[2], [0,0,0,0,1])
        np.testing.assert_allclose(law.sum(axis=1), 1., rtol=1e-14)
        with self.assertRaises(ValueError): screen.binomial_logs(4, [1.01])

    def test_mixture_and_epoch_count(self):
        prepared = screen.refresh.Epochs(4,2,{2:2,4:1},[1,0,2,0,1])
        epoch = prepared.at(.3)
        probabilities = np.exp(screen.binomial_logs(4,[.3]))[0]
        expected = sum(p*np.exp(e) for p,e in zip(probabilities,epoch))
        actual = screen.mixture(epoch, screen.binomial_logs(4,[.3]))
        np.testing.assert_allclose(np.exp(actual[0]),expected,rtol=1e-14,atol=1e-16)
        for count in (1,4,8):
            self.assertAlmostEqual(float(screen.refresh.terminal_logs(actual,count)[0]),
                math.log(np.linalg.matrix_power(expected,count)[0].sum()),places=12)

    def test_all_one_input_from_zero_is_exact(self):
        prepared = screen.refresh.Epochs(4,2,{2:2,4:1},[1,0,2,0,1])
        actual = screen.mixture(prepared.at(.7),screen.binomial_logs(4,[1.]))
        self.assertAlmostEqual(float(screen.refresh.terminal_logs(actual,16)[0]),-.7*64,places=12)


if __name__ == '__main__': unittest.main()
