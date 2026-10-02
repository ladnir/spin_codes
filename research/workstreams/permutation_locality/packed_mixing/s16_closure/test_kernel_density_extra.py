"""Exact representation, relative-accuracy and fail-closed regressions.

These tests enumerate tiny input spaces, not production proof domains.
The overflow cases deliberately supply adversarial synthetic Fourier data
to exercise arithmetic guards that physical probability laws should avoid.
"""
from fractions import Fraction as Q
from math import comb
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, ctx
import kernel_maps as kernel
import kernel_birth_density as density
import scalar_cover as sc
import sparse_kernel


def endpoint(value):
    mantissa, exponent = value.upper().man_exp()
    return Q(int(mantissa))*Q(2)**int(exponent)


def images_from_rows(rows):
    images = [0]*(1 << len(rows))
    for state in range(len(images)):
        for j, row in enumerate(rows):
            if state >> j & 1:
                images[state] ^= row
    return tuple(images)


def exact_conditional(columns, bits, z):
    """Direct all-input enumeration, independent of Fourier/profile code."""
    windows = len(columns)//4
    laws = [[Q(0)]*(1 << bits) for _ in range(windows+1)]
    for word in range(1 << len(columns)):
        occupancy = sum(bool((word >> offset) & 15)
                        for offset in range(0, len(columns), 4))
        state = 0
        for j, column in enumerate(columns):
            if word >> j & 1:
                state ^= column
        laws[occupancy][state] += z**word.bit_count()/(comb(windows, occupancy)*15**occupancy)
    return laws


def exact_iid(columns, bits, z, activity):
    """Independent product-law enumeration, without conditioning on J."""
    law = [Q(0)]*(1 << bits)
    for word in range(1 << len(columns)):
        probability = Q(1)
        for offset in range(0, len(columns), 4):
            probability *= activity/15 if (word >> offset) & 15 else 1-activity
        state = 0
        for j, column in enumerate(columns):
            if word >> j & 1:
                state ^= column
        law[state] += probability*z**word.bit_count()
    return law


class DensityExtraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.columns = (1, 2, 4, 3, 5, 7, 6, 1)
        cls.cases = []
        # The second A has several states in each expansion-weight class,
        # so the residual checks exercise class aggregation, not singletons.
        for rows in ((1, 6, 120), (0x11, 0x22, 0x44)):
            images = images_from_rows(rows)
            classes = kernel.prepare_maps(images, cls.columns, 3, birth_density='classes')
            capped = kernel.prepare_maps(images, cls.columns, 3, birth_density='capped')
            cls.cases.append((images, classes, capped))

    def setUp(self):
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)
        ctx.prec = 256

    def assert_representation(self, data, row, exact):
        """Zero mass plus U/L and per-class residual budgets cover m_s."""
        self.assertEqual(endpoint(row[1]), 0)
        self.assertGreaterEqual(endpoint(row[0]), exact[0])
        cap = endpoint(row[2])/((1 << data['bits'])-1)
        self.assertGreaterEqual(cap, 0)
        for target, level in enumerate(data['birth_class_levels'], 3):
            needed = sum((max(Q(0), exact[state]-cap)
                          for state in range(1, len(exact))
                          if data['map_images'][state].bit_count() == level), Q(0))
            self.assertGreaterEqual(endpoint(row[target]), needed)

    def test_capped_iid_rows_cover_exact_conditional_mixture(self):
        for _, classes, capped in self.cases:
            for z in (Q(1), Q(3, 4), Q(1, 64)):
                conditional = exact_conditional(self.columns, 3, z)
                za = kernel.aq(z)
                local = sparse_kernel.outward_at_z(capped, za)
                for activity in (Q(0), Q(1, 5), Q(3, 4), Q(1)):
                    with self.subTest(z=z, activity=activity,
                                      expansion=classes['map_sha256']):
                        masses = [Q(comb(2, j))*activity**j*(1-activity)**(2-j)
                                  for j in range(3)]
                        mixture = [sum((m*conditional[j][state] for j, m in enumerate(masses)), Q(0))
                                   for state in range(8)]
                        self.assertEqual(mixture, exact_iid(self.columns, 3, z, activity))
                        probabilities = sc.probabilities(activity)
                        baseline = kernel.outward_at_z(classes, probabilities, za)
                        chosen = kernel.outward_at_z(capped, probabilities, za)
                        self.assert_representation(capped, [chosen[0, j] for j in range(chosen.ncols())], mixture)

                        # Audit every offered IID row, not only the selected
                        # row. The oracle above never calls Fourier inversion.
                        values, _ = kernel.birth_classes.base.fourier_parts(capped, activity, za, 0)
                        nums, denominators = density.upper_laws(capped, [[v] for v in values])
                        original = [baseline[0, j] for j in range(baseline.ncols())]
                        for row in density.outward_rows(capped, original, nums[:, 0], denominators[0]):
                            self.assert_representation(capped, row, mixture)

                        # Independently chosen conditional rows need not
                        # match IID entries; their mixture must still cover
                        # exactly the same state measure.
                        refined = density.refine_local(capped, local, za, activity)
                        combined = [sum((kernel.aq(m)*refined[j][0, k]
                                         for j, m in enumerate(masses)), arb(0))
                                    for k in range(chosen.ncols())]
                        self.assert_representation(capped, combined, mixture)

    def test_per_column_scaling_preserves_tiny_mass_relative_accuracy(self):
        data = self.cases[0][2]
        ctx.prec = 1024
        for exponent in (80, 600):
            z = Q(1, 1 << exponent)
            exact = exact_conditional(self.columns, 3, z)
            numerators, denominators = density.conditional_laws(data, kernel.aq(z))
            # J=0,1,2 have radically different total tilted masses. A single
            # absolute 44-bit scale cannot satisfy this relative-error test.
            self.assertGreater(denominators[1], denominators[0]*(1 << (exponent-2)))
            self.assertGreater(denominators[2], denominators[1]*(1 << (exponent-2)))
            for j in range(3):
                total = sum(exact[j], Q(0))
                self.assertGreater(total, 0)
                for state in range(8):
                    upper = Q(int(numerators[state, j]), denominators[j])
                    self.assertGreaterEqual(upper, exact[j][state])
                    self.assertLessEqual(upper-exact[j][state], total/Q(1 << 38))

    def test_insufficient_precision_fails_before_walsh_then_fresh_retry_works(self):
        data = self.cases[0][2]
        ctx.prec = 32
        with patch.object(density, 'walsh', side_effect=AssertionError('unsafe transform reached')):
            with self.assertRaisesRegex(ArithmeticError, 'insufficient precision'):
                density.conditional_laws(data, arb(3)/5)
        self.assertEqual(ctx.prec, 32)
        # Explicit caller retry: no stale low-precision coefficient is reused.
        ctx.prec = 256
        numerators, denominators = density.conditional_laws(data, arb(3)/5)
        exact = exact_conditional(self.columns, 3, Q(3, 5))
        for j in range(3):
            for state in range(8):
                self.assertGreaterEqual(Q(int(numerators[state, j]), denominators[j]), exact[j][state])

    def test_signed_error_budget_overflow_is_rejected_before_walsh(self):
        # Deliberately corrupt multiplicities. A rounding-width sum that
        # cannot fit the audited budget must never reach the int64 transform.
        data = dict(bits=1, records=[(), ()], birth_character_indices=np.array([0, 1]),
                    multiplicities=[1, 1 << 61])
        with patch.object(density, 'walsh', side_effect=AssertionError('unsafe transform reached')):
            with self.assertRaisesRegex(OverflowError, 'signed Walsh accumulation'):
                density.upper_laws(data, [[arb(1)], [arb(1)/3]])

    def test_residual_sum_overflow_is_rejected_after_bounded_signed_transform(self):
        # A 14-bit bent sign function has a Walsh transform of magnitude
        # 128. Coefficients/intermediates fit int64, but the sum of positive
        # dyadic numerators exceeds the separate class-reduction budget.
        bits = 14
        indices = np.array([((x & (x >> 1) & 0x1555).bit_count() & 1)
                            for x in range(1 << bits)], dtype=np.int64)
        data = dict(bits=bits, records=[(), ()], birth_character_indices=indices,
                    multiplicities=np.bincount(indices, minlength=2))
        with self.assertRaisesRegex(OverflowError, 'class residual accumulation'):
            density.upper_laws(data, [[arb(1)], [arb(-1)]])


if __name__ == '__main__':
    unittest.main()
