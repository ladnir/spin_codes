import unittest

import numpy as np

import conditioned_diagnostics as cd


class DiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.local = np.array([[[.8, 0., 0.], [0., .7, 0.], [0., .6, 0.]],
            [[.1, 0., .4], [.2, .5, 0.], [.15, .4, 0.]],
            [[.2, 0., .3], [.1, .6, 0.], [.3, .2, 0.]]])

    def test_batched_ordered_placement(self):
        families = np.stack([self.local, self.local*1.13, self.local*np.array([1., .7, 1.2])[:, None, None]])
        for q in (0, 2, 5):
            moments, terminals = cd.batch_moments(families, q, epochs=3, windows=2, regions=2)
            for i, family in enumerate(families):
                expected = cd.logarithmic_moment(family, q, epochs=3, windows=2, regions=2)
                self.assertAlmostEqual(moments[i], expected, places=12)
                self.assertAlmostEqual(terminals[i].sum(), 1., places=13)

    def test_active_and_potential_derivatives(self):
        epsilon, nu, q = 1e-5, .7, 2
        families, labels = cd.make_perturbations(self.local, nu, epsilon=epsilon, packet_bits=2)
        moments, terminals = cd.batch_moments(families, q, epochs=3, windows=2, regions=2)
        counts = dict(zip(labels, (moments[2::2]-moments[1::2])/(2*epsilon)))
        self.assertAlmostEqual(sum(counts[f'potential_{j}'] for j in range(3)), 6., places=8)
        self.assertAlmostEqual(sum(j*counts[f'potential_{j}'] for j in range(3)), 4., places=8)
        self.assertAlmostEqual(sum(counts[f'active_{k}'] for k in range(3)), 6., places=8)
        self.assertAlmostEqual(counts['zero_source']+counts['nonzero_source'], 6., places=8)
        self.assertAlmostEqual(counts['birth_transition']-counts['return_transition'],
                               1-terminals[0, 0], places=8)
        # A block-upper-triangular transfer propagates an exact first
        # derivative independently of finite differences.
        for name, domain, selector in (
            ('potential_1', 'potential', (1, slice(None), slice(None))),
            ('active_1', 'active', (1, slice(None), slice(None))),
            ('zero_source', 'active', (slice(None), 0, slice(None))),
            ('active_zero_feedback_zero_to_zero', 'active', (slice(1, None), 0, 0))):
            base = cd.cap_gate.marked(self.local, nu, packet_bits=2)
            derivative = np.zeros_like(self.local)
            derivative[selector] = (base if domain == 'potential' else self.local)[selector]
            if domain == 'active':
                derivative = cd.cap_gate.marked(derivative, nu, packet_bits=2)
            dual = np.zeros((3, 6, 6))
            dual[:, :3, :3] = dual[:, 3:, 3:] = base
            dual[:, :3, 3:] = derivative
            regional, _ = cd.cap_gate.gate.prior.placement(dual, q, epochs=3, windows=2, force_log=True)
            total = np.linalg.matrix_power(np.exp(regional[q]), 2)
            expected = total[0, 3:].sum()/total[0, :3].sum()
            self.assertAlmostEqual(counts[name], expected, places=8)


if __name__ == '__main__':
    unittest.main()
