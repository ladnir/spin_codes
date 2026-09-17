"""Exact small-state checks for the new shell-supported transfer rows."""
from collections import Counter
from fractions import Fraction as F
import unittest

import q1_injection_shells as shells
import q1_slack


class InjectionShellTests(unittest.TestCase):
    def test_coefficient_minimum_requires_identical_coverage(self):
        self.assertEqual(q1_slack.minimum({1: F(1, 2)}, {1: F(1, 3)}), {1: F(1, 3)})
        with self.assertRaises(ValueError):
            q1_slack.minimum({1: F(1, 2)}, {2: F(1, 3)})

    def test_new_rows_dominate_exact_small_state_laws(self):
        a = [1, 2, 3, 4, 5]
        b = [3, 5, 6, 7, 1]
        t, m = len(a), 7
        image = lambda q: sum(((q & col).bit_count() & 1) << i for i, col in enumerate(a))
        weights = {q: image(q).bit_count() for q in range(1, m+1)}
        spectrum = Counter(weights.values())
        levels = sorted(spectrum)
        old_n = len(levels) + 2
        iw = [weights[q] for q in b]
        cw = [(image(q) ^ (1 << i)).bit_count() for i, q in enumerate(b)]
        for z in (F(1, 2), F(9, 10), F(1)):
            base0 = [F()] * (old_n**2)
            base1 = base0.copy()
            base0[0], base1[1] = F(1), z
            rz, ra, n = shells.extend(base0, base1, spectrum, iw, cw, z, t)
            for j, transfer in enumerate((rz, ra)):
                for q in range(m+1):
                    if q == 0 and j == 0:
                        continue
                    index = 0 if q == 0 else old_n + levels.index(weights[q])
                    bound = transfer[index*n:(index+1)*n]
                    exact = [F()] * (m+1)
                    inputs = [(0, 0)] if j == 0 else [(1 << i, b[i]) for i in range(t)]
                    for x, feedback in inputs:
                        mass = z**((x ^ image(q)).bit_count()) / len(inputs)
                        if q == 0:
                            exact[feedback] += mass
                        else:
                            exact[q ^ feedback] += mass/2
                            for u in range(1, m+1):
                                exact[u ^ feedback] += mass/(2*m)
                    self.assertLessEqual(exact[0], bound[0])
                    residual = {u: max(F(), exact[u] - bound[2+levels.index(weights[u])]/spectrum[weights[u]])
                                for u in range(1, m+1)}
                    uncovered = sum((max(F(), sum(residual[u] for u in residual if weights[u] == v)
                                               - bound[old_n+i]) for i, v in enumerate(levels)), F())
                    self.assertLessEqual(uncovered, bound[1], (z, j, q))

    def test_zero_activation_preserves_exact_shell_masses(self):
        spectrum = {1: 2, 2: 1}
        n = len(spectrum) + 2
        zero, one = [F()]*(n*n), [F()]*(n*n)
        zero[0], one[1] = F(1), F(3, 4)
        rz, ra, width = shells.extend(zero, one, spectrum, [1, 1, 2], [0, 2, 1], F(3, 4), 3)
        self.assertEqual(ra[1], 0)
        self.assertEqual(ra[n:n+2], (F(1, 2), F(1, 4)))
        self.assertEqual(sum(ra[:width]), F(3, 4))


if __name__ == '__main__':
    unittest.main()
