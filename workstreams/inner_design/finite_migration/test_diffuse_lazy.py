"""Exact small-support checks of the diagnostic lazy-branch envelope."""
from collections import Counter
from fractions import Fraction as F
import math
import unittest

import numpy as np
import parameter_diffuse_lazy as subject


class DiffuseTests(unittest.TestCase):
    def test_pointwise_measure_dominance(self):
        # No relation between A and C is required. Include rank-deficient C.
        for columns_a, columns_c, s in (
            ([1, 2, 4, 3, 5, 6], [1, 3, 5, 2, 7, 6], 3),
            ([1, 2, 3, 1], [1, 1, 1, 1], 2),
        ):
            t, m = len(columns_a), (1 << s) - 1
            images = [sum(((q & c).bit_count() % 2) << i
                          for i, c in enumerate(columns_a)) for q in range(m + 1)]
            shells = {}
            for q in range(1, m + 1):
                shells.setdefault(images[q].bit_count(), []).append(q)
            levels = sorted(shells)
            spectrum = {v: len(shells[v]) for v in levels}
            supports = [[] for _ in range(t + 1)]
            for x in range(1 << t):
                syndrome = 0
                for i, c in enumerate(columns_c):
                    if x >> i & 1:
                        syndrome ^= c
                supports[x.bit_count()].append((x, syndrome))
            fibers = [Counter(c for _, c in row) for row in supports]
            kernel = [row[0] for row in fibers]
            caps = [{'cap': max((n for c, n in row.items() if c), default=0)}
                    for row in fibers]
            for z in (F(1, 4), F(1, 2), F(3, 4), F(1)):
                epochs = np.full((t + 1, len(levels) + 2, len(levels) + 2), -np.inf)
                for j, inputs in enumerate(supports):
                    for i, v in enumerate(levels):
                        marginal = sum(z ** (images[shells[v][0]] ^ x).bit_count()
                                       for x, _ in inputs) / len(inputs)
                        epochs[j, i + 2, 1] = math.log(float(marginal / 2))
                result = subject.diffuse(epochs, spectrum, kernel, caps, t, -math.log(float(z)))
                for j, inputs in enumerate(supports):
                    if kernel[j] == len(inputs):
                        np.testing.assert_array_equal(result[j], epochs[j])
                        continue
                    for i, v in enumerate(levels):
                        a = len(shells[v])
                        exact = [F() for _ in range(m + 1)]
                        lazy_nonzero = [F() for _ in range(m + 1)]
                        for q in shells[v]:
                            for x, c in inputs:
                                weight = z ** (images[q] ^ x).bit_count() / (a * len(inputs))
                                exact[q ^ c] += weight / 2
                                if c:
                                    lazy_nonzero[q ^ c] += weight
                                for u in range(1, m + 1):
                                    exact[u ^ c] += weight / (2 * m)
                        nu = F(len(inputs) - kernel[j], len(inputs))
                        h = F(caps[j]['cap'], len(inputs))
                        point_cap = min(nu, a * h) / a * z ** abs(v - j)
                        for y in range(m + 1):
                            self.assertLessEqual(lazy_nonzero[y], point_cap)
                        self.assertEqual(result[j, i + 2, 1], -np.inf)
                        for k, w in enumerate(levels):
                            envelope = math.exp(result[j, i + 2, k + 2]) / len(shells[w])
                            for y in shells[w]:
                                self.assertLessEqual(float(exact[y]), envelope + 2e-14)

    def test_zero_feedback_leaves_existing_rows(self):
        epochs = np.arange(27, dtype=float).reshape(3, 3, 3)
        result = subject.diffuse(epochs, {1: 3}, [1, 2, 1],
                                 [{'cap': 0}] * 3, 2, .5)
        np.testing.assert_array_equal(result, epochs)
        self.assertIsNot(result, epochs)


if __name__ == '__main__':
    unittest.main()
