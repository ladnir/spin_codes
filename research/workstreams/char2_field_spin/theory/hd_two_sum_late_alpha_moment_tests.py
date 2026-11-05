from fractions import Fraction as F
from itertools import combinations, product
from math import comb, log2
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_two_sum_late_alpha_moment as gate
from hd_two_sum_closed_kernel import zero_two_moments


def multiply(left, right):
    answer = [F(0)]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i+j] += a*b
    return answer


def product_polynomials(rows):
    answer = [F(1)]
    for row in rows:
        answer = multiply(answer, row)
    return answer


def polynomial(n, moments, q=4):
    return [comb(n, m)*(q-1)**m*moments[m] for m in range(n+1)]


class LateAlphaMomentTests(TestCase):
    def test_quiet_positive_majorization(self):
        for q in (4, 8, 256):
            eta = F(1, (q-1)**2)
            for w in range(1, min(q, 16)+1):
                self.assertLessEqual(zero_two_moments(w, q), eta*(q-1)**w)
            self.assertEqual(1-eta+eta, 1)

    def test_exact_common_A_mixture_identity(self):
        q, a, r = 4, 3, 2
        eta, moments = F(1, 9), [F(1), F(1, 4), F(1, 8)]
        direct = [F(0)]*(a*r+1)
        for masks in product(range(1 << a), repeat=r):
            weights = [mask.bit_count() for mask in masks]
            mass = F(1)
            for weight in weights:
                mass *= F(1) if weight == 0 else eta*(q-1)**weight
            for row in range(a):
                mass *= moments[sum(mask >> row & 1 for mask in masks)]
            direct[sum(weights)] += mass
        mixed = [F(0)]*(a*r+1)
        for h in range(r+1):
            term = product_polynomials([polynomial(h, moments)]*a)
            weight = comb(r, h)*(1-eta)**(r-h)*eta**h
            for k, value in enumerate(term):
                mixed[k] += weight*value
        self.assertEqual(direct, mixed)
        actual = gate.quiet_mixture_coefficients(r, a, np.log2(list(map(float, moments))), q=q)
        np.testing.assert_allclose(2**actual, list(map(float, direct)), rtol=1e-12)

    def test_strong_coefficientwise_moment_convexity(self):
        atoms = ((F(0), F(2, 5)), (F(1), F(1, 5)), (F(3), F(2, 5)))
        moments = [sum((probability*value**m for value, probability in atoms), F(0))
                   for m in range(7)]
        # At m0, including0**0=1 represents a probability moment correctly.
        self.assertEqual(moments[0], 1)
        for a in range(1, 6):
            for b in range(1, a+1):
                large = multiply(polynomial(a+1, moments), polynomial(b-1, moments))
                small = multiply(polynomial(a, moments), polynomial(b, moments))
                self.assertTrue(all(x >= y for x, y in zip(large, small)))

    def test_arbitrary_fixed_rowsets_endpoint_bound(self):
        # Enumerate EVERY labelled pair of a2 row sets on four regions.
        q, group, a, r = 4, 4, 2, 2
        moments = [F(1), F(1, 5)*2, F(1, 5)*4]
        sets = tuple(combinations(range(group), a))
        for rowsets in product(sets, repeat=r):
            for selected in product((0, 1), repeat=r):
                h = sum(selected)
                incidence = [sum(active and region in rowset
                                 for active, rowset in zip(selected, rowsets)) for region in range(group)]
                self.assertEqual(sum(incidence), a*h)
                actual = product_polynomials([polynomial(n, moments) for n in incidence])
                bound = product_polynomials([polynomial(h, moments)]*a)
                self.assertTrue(all(x <= y for x, y in zip(actual, bound)))

    def test_two_atom_curve_and_guards(self):
        curve = np.array([0., -10., -12., -14.])
        moments, metadata = gate.two_atom_envelope(curve, 1, 3)
        self.assertEqual(metadata['log2_nonzero_atom'], -2.)
        self.assertEqual(metadata['log2_probability'], -8.)
        np.testing.assert_array_equal(moments, curve)
        with self.assertRaises(ValueError):
            gate.two_atom_envelope(curve, 0, 3)
        with self.assertRaises(ArithmeticError):
            gate.two_atom_envelope(np.array([0., 2., 2.]), 1, 2)
        with self.assertRaises(ValueError):
            gate.run(max_seconds=31)

    def test_two_scopes_do_not_share_rowset_factor(self):
        # Test the production accounting without evaluating a production bank.
        curve = np.arange(gate.GROUPS+1, dtype=float)*-1.
        source = dict(regional_envelope={'log2_slope': -1.},
                      comparison={'late_column_log2': -10.})
        with patch.object(gate, 'read_source', return_value=(source, curve)):
            with patch.object(gate, 'column_coefficients', return_value=(np.array([-10.]), None)):
                word = np.full(gate.MINIMUM*gate.ARITY*gate.GROUPS+1, -np.inf)
                word[-1] = -10.
                with patch.object(gate, 'log_power', return_value=word):
                    with patch.object(gate, 'quiet_mixture_coefficients', return_value=np.zeros(145)):
                        with patch.object(gate, 'word_tail', return_value=-20.):
                            result = gate.run()
        common, arbitrary = result['controls']
        difference = (gate.GROUPS-1)*log2(comb(16, gate.ARITY))
        self.assertAlmostEqual(arbitrary['log2_prefix']-common['log2_prefix'], difference)
        self.assertFalse(common['arbitrary_fixed_rowset_vectors'])
        self.assertTrue(arbitrary['arbitrary_fixed_rowset_vectors'])
        self.assertEqual(result['raw_regional_banks_built'], 0)
        self.assertFalse(result['whole_code_certificate'])


if __name__ == '__main__':
    main()
