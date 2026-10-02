import unittest
from fractions import Fraction as Q
from math import log

import numpy as np
from flint import arb, arb_mat, ctx
from packet_outer_cost_search import Model, estimate, replay, geometry
from rs_uniform_envelope import UniformInputEnvelope


class CostSearchTests(unittest.TestCase):
    def test_scalar_formula_including_sparse(self):
        env = UniformInputEnvelope(2, 1, 4, 1)
        local = [arb_mat([[arb(9)/10*(arb(4)/5)**j]]) for j in range(3)]
        got = estimate(local, K=32, envelope=env, occupancies=[1, 2, 4],
                       tilt=Q(1, 20), windows=2)
        for q in (1, 2, 4):
            expected = 2*(4*log(.9)+q*log((1+15*.8)/16))
            self.assertAlmostEqual(got['witnesses'][str(q)]['log_moment'], expected, places=11)
        self.assertTrue(got['proposal_only'])
        self.assertFalse(got['whole_code_certificate'])

    def test_continuous_state(self):
        numeric = np.array([[.8, .1], [.02, .6]])
        matrix = arb_mat([[arb(4)/5, arb(1)/10], [arb(1)/50, arb(3)/5]])
        env = UniformInputEnvelope(2, 1, 4, 1)
        got = estimate([matrix]*3, K=16, envelope=env, occupancies=[1],
                       tilt=Q(1, 20), windows=2)
        expected = log(np.linalg.matrix_power(numeric, 4)[0].sum())
        reset = 2*log(np.linalg.matrix_power(numeric, 2)[0].sum())
        self.assertAlmostEqual(got['witnesses']['1']['log_moment'], expected, places=12)
        self.assertGreater(abs(expected-reset), .001)

    def test_arb_replay_agrees_with_proposal(self):
        old = ctx.prec
        try:
            ctx.prec = 192
            env = UniformInputEnvelope(2, 1, 4, 1)
            local = [arb_mat([[arb(9)/10*(arb(4)/5)**j]]) for j in range(3)]
            options = dict(K=32, envelope=env, occupancies=[1, 2, 4], tilt=Q(1, 20), windows=2)
            screen, exact = estimate(local, **options), replay(local, **options)
            for q in ('1', '2', '4'):
                m, e = exact['values'][q]['upper']
                margin = -(log(m)/log(2)+e)
                self.assertAlmostEqual(margin, screen['witnesses'][q]['estimated_margin_bits'], places=10)
            self.assertFalse(exact['whole_code_certificate'])
        finally:
            ctx.prec = old

    def test_no_underflow_of_local_mass(self):
        tiny = arb(2)**-10000
        env = UniformInputEnvelope(2, 1, 4, 1)
        result = estimate([arb_mat([[tiny]])]*3, K=16, envelope=env,
                          occupancies=[1], tilt=Q(1, 20), windows=2)
        self.assertAlmostEqual(result['witnesses']['1']['log_moment'], -40000*log(2), places=8)
        fallback = estimate([arb_mat([[tiny]])]*3, K=16, envelope=env,
                            occupancies=[1], tilt=Q(1, 20), windows=2, backend='scaled-or-log')
        self.assertEqual(fallback['backend'], 'log')
        self.assertEqual(fallback['witnesses'], result['witnesses'])

    def test_fast_scaled_agrees(self):
        env = UniformInputEnvelope(2, 1, 4, 1)
        local = [arb_mat([[arb(9)/10*(arb(4)/5)**j]]) for j in range(3)]
        options = dict(K=32, envelope=env, occupancies=[1, 2, 4], tilt=Q(1, 20), windows=2)
        slow = estimate(local, **options)
        fast = estimate(local, backend='scaled-or-log', **options)
        self.assertEqual(fast['backend'], 'scaled')
        for q in ('1', '2', '4'):
            self.assertAlmostEqual(slow['witnesses'][q]['estimated_margin_bits'],
                                   fast['witnesses'][q]['estimated_margin_bits'], places=10)

    def test_geometry_and_rejection(self):
        env = UniformInputEnvelope(16, 8, 4, 8)
        result = geometry(1 << 20, env)
        self.assertEqual((result['groups'], result['regions'], result['epochs_per_region']), (4096, 128, 128))
        with self.assertRaises(ValueError):
            geometry((1 << 20)+128, env)
        with self.assertRaises(ValueError):
            geometry(1024, UniformInputEnvelope(3, 1, 4, 1))

    def test_authenticated_model_cannot_change_macro_width(self):
        # The rejection precedes access to prepared data or any census.
        model = Model.__new__(Model)
        for function in (model.estimate, model.replay):
            with self.assertRaises(ValueError):
                function(tilt=Q(1, 10), windows=16)


if __name__ == '__main__':
    unittest.main()
