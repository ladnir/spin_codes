"""Tiny exact/Arb oracles; no provider or production contraction."""
from fractions import Fraction as F
from itertools import combinations
from math import nextafter, prod
import unittest

import numpy as np
from flint import arb, ctx

import hd_bit_syndrome_atoms_directed as atoms
import hd_positive_upward as up


def exact(value):
    return atoms.upper_fraction(value)


class DirectedAtomTests(unittest.TestCase):
    def test_pinned_frozen_literals_only(self):
        self.assertEqual(atoms.frozen_parameters(), atoms.PARAMETERS)

    def test_dyadic_rounding_and_checked_roots(self):
        with ctx.workprec(256):
            for value in (F(0), F(1, 3), F(2), F(1, 1 << 1500), F(1 << 1600)):
                upper = atoms.dyadic_up(value, 80)
                self.assertGreaterEqual(upper, value)
                if value:
                    self.assertLessEqual(upper/value, 1+F(1, 1 << 78))
            for base, power in ((F(2), F(3, 7)), (F(0), F(2, 3)),
                                (F(1, 17), F(5, 2)), (F(7, 3), F(0))):
                upper = exact(atoms.power_upper(base, power))
                self.assertGreaterEqual(upper**power.denominator, base**power.numerator)

    def test_fresh_curves_dominated_and_probability_complement(self):
        with ctx.workprec(256):
            caps = np.empty((4, 19), dtype=object)
            for m in range(19):
                for e, value in enumerate((F(1) if not m else F(1, 1 << (25+7*m)),
                                           F(1, 1 << (1+2*m)), F(1, 1 << (125+2*m)), F(1, 1 << 68))):
                    caps[e, m] = atoms.rational(value).upper()
            result = atoms.build_atoms(caps)
            a, w = result['atoms'], result['masses']
            self.assertEqual(w[0]+w[1], 1)
            for m in range(19):
                upper = (F(1) if not m else w[1]*a[1]**m, w[2]*a[2]**m,
                         w[3]*a[3]**m, w[4])
                for e in range(4):
                    self.assertGreaterEqual(upper[e], exact(caps[e, m]))
            self.assertTrue(result['all_76_moments_verified'])
            self.assertFalse(result['saved_float_bounds_used'])

    def test_failed_zero_moment_and_probability_guard(self):
        with ctx.workprec(256):
            curves = np.array([[arb(0) for _ in range(19)] for _ in range(4)], dtype=object)
            curves[0, 0] = arb(2)
            with self.assertRaises(ValueError):
                atoms.build_atoms(curves)
            curves[0, 0] = arb(1)
            curves[0, 1] = arb(1)
            with self.assertRaises(ValueError):
                atoms.build_atoms(curves)

    def test_cover_weights_against_independent_arb_oracle(self):
        measure = dict(atoms=(F(0), F(1, 3), F(2), F(3, 5), F(1)),
            masses=(F(2, 3), F(1, 3), F(2, 5), F(1, 11), F(3, 7)), precision=256)
        covers = (F(0), F(1, 3), F(2, 3), F(1))
        result = atoms.primitive_weights(measure, covers, ((1, 2, 3),), ((0, 1, 2),), groups=2)
        with ctx.workprec(384):
            for ci, c in enumerate(covers):
                oracle = [[arb(0) for _ in range(2)] for _ in range(2)]
                for t, mass, (i, j) in zip(measure['atoms'], measure['masses'], atoms.ENDPOINTS):
                    term = ((1+255*atoms.rational(t)**atoms.rational(1/c))**atoms.rational(2*c)
                            if c else atoms.rational(max(1, t)**2))
                    oracle[i][j] += atoms.rational(mass)*term
                for i in range(2):
                    for j in range(2):
                        self.assertGreaterEqual(exact(result['cover_matrices'][ci, i, j]), exact(oracle[i][j]))
            self.assertEqual(exact(result['empty_matrix'][0, 0]), 1)

    def test_terminal_exact_labelled_masks_and_no_edge_mass(self):
        measure = dict(atoms=(F(0), F(1, 3), F(2), F(3, 5), F(1)),
            masses=(F(2, 3), F(1, 3), F(2, 5), F(1, 11), F(3, 7)), precision=256)
        quiet, outside = ((1, 2, 3), (2, 0, 7)), ((0, 1, 2), (0, 3, 1))
        result = atoms.primitive_weights(measure, (F(0),), quiet, outside, groups=3)
        for ci, count in enumerate(result['count_order']):
            markers = tuple(t for t, number in zip(measure['atoms'], count) for _ in range(number))
            elementary = [sum((prod(row, start=F(1)) for row in combinations(markers, degree)), F(0))
                          for degree in range(3)]
            for name, table in (('terminal_quiet', quiet), ('terminal_nonquiet', outside)):
                for cls, coefficients in enumerate(table):
                    expected = sum(n*e for n, e in zip(coefficients, elementary))**3
                    self.assertGreaterEqual(exact(result[name][cls, ci]), expected)
                    self.assertLessEqual(exact(result[name][cls, ci])-expected,
                                         max(F(1), expected)*F(1, 1 << 220))
        self.assertFalse(result['edge_masses_in_terminals'])

    def test_scaled_export_normal_subnormal_and_zero(self):
        with ctx.workprec(256):
            values = [F(0), F(1, 3), F(2)**1024, F(2)**-1500]
            raw = np.array([atoms.rational(x).upper() for x in values], dtype=object)
            result = atoms.export_scaled(raw)
            self.assertEqual(result['mantissas'][0], 0.)
            self.assertGreater(result['mantissas'][-1], 0.)
            for x, y in zip(raw, result['mantissas']):
                self.assertGreaterEqual(F(float(y))*F(2)**result['exponent'], exact(x))
            tiny = atoms.export_scaled([atoms.rational(F(2)**-1500)], exponent=0)
            self.assertEqual(tiny['mantissas'][0], nextafter(0., 1.))
            exact_float = nextafter(1., 2.)
            self.assertEqual(atoms.export_scaled([arb(exact_float)], exponent=0)['mantissas'][0], exact_float)
            with self.assertRaises(OverflowError):
                atoms.export_scaled([atoms.rational(F(2)**2000)], exponent=0)

    def test_cross_layer_export_after_precision_context(self):
        with ctx.workprec(384):
            endpoints = np.array([arb(2).sqrt().upper(),
                                  (arb(2)**-1500).upper(), arb(0)], dtype=object)
            original = tuple(exact(value) for value in endpoints)
        # Neither reading an exact endpoint nor converting the exported
        # dyadic floats may depend on the caller's ambient Arb precision.
        with ctx.workprec(32):
            self.assertEqual(tuple(exact(value) for value in endpoints), original)
            exported = atoms.export_scaled(endpoints)
            imported = up.from_export(exported)
            for wanted, value in zip(original, imported.mantissas):
                self.assertGreaterEqual(F(float(value))*F(2)**imported.exponent, wanted)
            self.assertEqual(imported.mantissas[-1], 0.)

    def test_guards(self):
        with self.assertRaises(ValueError):
            atoms.upper_fraction(arb(-1))
        with self.assertRaises(ValueError):
            atoms.export_scaled([])
        with self.assertRaises(ValueError):
            atoms.build_atoms(np.array([[arb(0)]]))
        with self.assertRaises(ValueError):
            atoms.power_upper(F(-1), F(1, 2))


if __name__ == '__main__':
    unittest.main()
