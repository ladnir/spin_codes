from fractions import Fraction as F
import math
import unittest

import numpy as np
from flint import arb, ctx

import certify_density_types as cert


class TypeArbTests(unittest.TestCase):
    def test_exact_proposal_normalization(self):
        pi, ps = cert.probabilities(dict(proposal=[.1, .2, .700000000000001], probabilities=[0., .4, 1.]))
        self.assertEqual(sum(pi), 1)
        self.assertEqual(ps[1], F(.4))
        with self.assertRaises(ValueError):
            cert.probabilities(dict(proposal=[0., 1.], probabilities=[0., .5]))

    def test_toy_moment_against_independent_log_transfer(self):
        ctx.prec = 256
        checker = cert.Checker.__new__(cert.Checker)
        checker.spec = dict(rows=8, cutoff=204)
        checker.t, checker.s = 4, 2
        checker.ac, checker.kernel = {2: 2, 4: 1}, {0: 1, 2: 2, 4: 1}
        log_tilt, theta = F(-1), F(1, 3)
        actual = checker.moment(log_tilt, theta)
        model = cert.transfer.reference.Epochs(4, 2, checker.ac, [1, 0, 2, 0, 1])
        matrix = cert.transfer.reference.base.epoch_mixture(model.at(math.exp(-1)), 4, [float(theta)])
        expected = cert.transfer.reference.base.terminal_logs(matrix, 512)[0]+204*math.exp(-1)
        self.assertAlmostEqual(float(actual), float(expected), places=10)


if __name__ == '__main__':
    unittest.main()
