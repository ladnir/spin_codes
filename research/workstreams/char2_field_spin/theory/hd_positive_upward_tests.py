"""Exact-rational controls for positive directed binary64 arithmetic."""
from fractions import Fraction as F
from random import Random
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_positive_upward as up


def scalar(value, exponent=0):
    answer = F(float(value))
    return answer * (1 << exponent) if exponent >= 0 else answer / (1 << -exponent)


def exact(value):
    return np.array([scalar(x, value.exponent) for x in value.mantissas.flat], dtype=object).reshape(value.mantissas.shape)


class UpwardTests(TestCase):
    def assert_upper(self, got, wanted):
        values = exact(got) if isinstance(got, up.Scaled) else np.array([F(float(x)) for x in np.asarray(got).flat], dtype=object).reshape(np.asarray(got).shape)
        wanted = np.asarray(wanted, dtype=object)
        self.assertEqual(values.shape, wanted.shape)
        self.assertTrue(np.all(values >= wanted))

    def test_platform_and_subnormal_primitives(self):
        self.assertTrue(up.ieee_checks()['gradual_underflow'])
        tiny = np.nextafter(0., 1.)
        self.assertEqual(up.upper_multiply(tiny, tiny), tiny)
        self.assertEqual(up.upper_ldexp(tiny, -3000), tiny)
        self.assertEqual(up.upper_multiply(0., tiny), 0)
        self.assertEqual(up.upper_add(0., 0.), 0)
        self.assertEqual(up.upper_ldexp(0., -3000), 0)

    def test_random_exact_positive_primitives(self):
        rng = Random(90210)
        for _ in range(256):
            a = np.ldexp(0.5 + rng.getrandbits(52) * 2.**-53, rng.randrange(-1000, 800))
            b = np.ldexp(0.5 + rng.getrandbits(52) * 2.**-53, rng.randrange(-1000, 200))
            self.assert_upper(up.upper_add(a, b), F(float(a)) + F(float(b)))
            self.assert_upper(up.upper_multiply(a, b), F(float(a)) * F(float(b)))
            shift = rng.randrange(-1000, 100)
            self.assert_upper(up.upper_ldexp(a, shift), scalar(a, shift))

    def test_normalization_and_large_exact_scales(self):
        a = up.normalize([0., 0.125, 0.75], 5000)
        b = up.normalize([0., 0.625, 0.25], -5000)
        self.assert_upper(up.add(a, b), exact(a) + exact(b))
        self.assert_upper(up.multiply(a, b), exact(a) * exact(b))
        self.assert_upper(up.concatenate((a, b)), np.concatenate((exact(a), exact(b))))
        self.assert_upper(up.add(up.zeros((3,)), b), exact(b))
        self.assertEqual(up.normalize([0.], 5000).exponent, 0)

    def test_broadcast_products_and_noncommuting_matrix_power(self):
        a = up.normalize(np.array([[[.5, .125], [.25, .75]], [[.25, .5], [0., .125]]]))
        b = up.normalize(np.array([[.375, .125], [.5, .25]]))
        want = np.stack([matrix @ exact(b) for matrix in exact(a)])
        self.assert_upper(up.matmul(a, b), want)
        self.assertFalse(np.array_equal(exact(a)[0] @ exact(b), exact(b) @ exact(a)[0]))
        result = up.matrix_power(a, 16)
        wanted = []
        for matrix in exact(a):
            product = np.eye(2, dtype=object)
            for _ in range(16):
                product = product @ matrix
            wanted.append(product)
        self.assert_upper(result, np.array(wanted, dtype=object))

    def test_explicit_pairwise_sum(self):
        rng = Random(44)
        data = np.array([[rng.random(), rng.random()] for _ in range(495)])
        value = up.normalize(data, -200)
        wanted = np.array([sum(exact(value)[:, i], F(0)) for i in range(2)], dtype=object)
        self.assert_upper(up.sum_axis0(value), wanted)

    def test_injective_scatter_and_input_preservation(self):
        dest = up.normalize([[.5, .25], [.125, 0.], [0., .75]], 100)
        values = up.normalize([[.25, .75], [.125, .5]], -50)
        before = dest.mantissas.copy()
        wanted = exact(dest).copy()
        wanted[[2, 0]] += exact(values)
        self.assert_upper(up.injective_add(dest, [2, 0], values), wanted)
        np.testing.assert_array_equal(dest.mantissas, before)
        with self.assertRaises(ValueError):
            up.injective_add(dest, [0, 0], values)

    def test_invalid_and_overflow_fail_closed(self):
        for value in (2**53+1, F(1, 3), np.array([1], dtype=np.float32)):
            with self.assertRaises(ValueError):
                up.normalize(value)
        for value in (-1., np.nan, np.inf):
            with self.assertRaises(ValueError):
                up.normalize(value)
        with self.assertRaises(ValueError):
            up.normalize(1., 0.5)
        with self.assertRaises((ArithmeticError, FloatingPointError)):
            up.upper_multiply(np.finfo(float).max, 2.)
        with self.assertRaises(ValueError):
            up.matrix_power(up.identity(), -1)

    def test_flush_to_zero_probe_fails_closed(self):
        real_multiply = np.multiply
        def flushing_multiply(a, b):
            result = real_multiply(a, b)
            return np.where((result > 0) & (result < np.finfo(float).tiny), 0., result)
        with patch.object(np, 'multiply', side_effect=flushing_multiply):
            with self.assertRaises(ArithmeticError):
                up.ieee_checks()


if __name__ == '__main__':
    main()
