"""Positive finite-sum checks of the source envelope, not certificates."""

import itertools
import math
import unittest

import numpy as np

import random_packet_screen as packet
import source_envelope as envelope


class SourceEnvelopeTests(unittest.TestCase):
    def test_entrywise_domination(self):
        regions = packet.log_regions(.72, 8, rows=16, width=4, state=2, identity=True)
        for logv in (-5., -1., 0., 1.):
            matrix = envelope.dominate(regions, logv)
            self.assertTrue(np.all(regions <= matrix+np.arange(9)[:, None, None]*logv+1e-14))

    def test_actual_source_moment_is_bounded(self):
        for active, z, logv in itertools.product((1, 2, 4, 8), (.4, .9, 1.), (-5., -1., 0.)):
            regions = packet.log_regions(z, active, rows=16, width=4, state=2, identity=True)
            actual = packet.positive_source_moment(regions, active, 6, 4)
            upper = envelope.source_envelope(regions, active, 6, 4, logv)
            self.assertLessEqual(actual, upper+1e-10)

    def test_arbitrary_positive_matrices(self):
        rng = np.random.default_rng(3927)
        for active in (1, 2, 4):
            regions = np.log(rng.uniform(.01, 1., size=(active+1, 2, 2)))
            exact = packet.positive_source_moment(regions, active, 5, 3)
            for lv, pv in itertools.product((-2., 0., .8), repeat=2):
                upper = envelope.source_envelope(regions, active, 5, 3, lv, pv)
                self.assertLessEqual(exact, upper+1e-10)

    def test_nonzero_enumerator(self):
        for k, lv in itertools.product((1, 4, 32), (-20., -2., 0.)):
            direct = sum(math.comb(k, w)*math.exp(w*(math.log(15)+lv)) for w in range(1, k+1))
            self.assertAlmostEqual(envelope.log_nonzero_binomial(k, lv), math.log(direct), places=11)


if __name__ == '__main__':
    unittest.main()
