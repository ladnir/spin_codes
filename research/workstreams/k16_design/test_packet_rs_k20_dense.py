"""Lightweight exact algebra checks; no selected-map census or dense proof run."""
from fractions import Fraction as Q
from itertools import product
from math import comb
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat
import packet_rs_k20_dense as dense


def rational(value):
    return fmpq(value.numerator, value.denominator)


def endpoint_fraction(value):
    mantissa, exponent = dense.q1.endpoint(value)
    return Q(mantissa) * Q(2)**exponent


def exact_iid(local, probability):
    return dense.iid_macro(local, probability, matrix=fmpq_mat,
                           rational=rational, rounding=lambda x: x)


def entrywise_le(left, right):
    return all(left[i, j] <= right[i, j]
               for i in range(left.nrows()) for j in range(left.ncols()))


class DenseConditioningTests(unittest.TestCase):
    def test_k20_geometry_and_rejections(self):
        g = dense.geometry()
        self.assertEqual((g.K, g.N, g.group_count, g.regions), (2**20, 2**21, 8192, 64))
        self.assertEqual(g.macros_per_region, 256)
        self.assertEqual(g.regions * g.macros_per_region, 16384)
        self.assertEqual(g.N // 10, 209715)
        self.assertEqual(dense.geometry(512).macros_per_region, 16)
        for invalid in (0, -32, 31, 33, True, 32.0):
            with self.assertRaises(ValueError):
                dense.geometry(invalid)

    def test_exact_conditioning_probability_and_symmetry(self):
        for groups in range(1, 17):
            for q in range(groups + 1):
                p = Q(q, groups)
                probability = comb(groups, q) * p**q * (1 - p)**(groups - q)
                loss = dense.conditioning_loss_exact(groups, q)
                self.assertEqual(loss * probability, 1)
                self.assertGreaterEqual(loss, 1)
                self.assertEqual(loss, dense.conditioning_loss_exact(groups, groups - q))
        for groups, q in ((0, 0), (4, -1), (4, 5), (True, 0), (4, True)):
            with self.assertRaises(ValueError):
                dense.conditioning_loss_exact(groups, q)

    def test_arb_conditioning_loss_encloses_exact_loss(self):
        ctx.prec = 128
        for groups in (1, 2, 7, 16, 32):
            for q in range(groups + 1):
                upper = endpoint_fraction(dense.log_conditioning_loss(groups, q).exp())
                self.assertGreaterEqual(upper, dense.conditioning_loss_exact(groups, q))
        self.assertEqual(dense.log_conditioning_loss(8192, 8192), arb(0))

    def test_binomial_mixture_and_endpoint_weights(self):
        local = [fmpq_mat([[1, j], [j + 1, 2]]) for j in range(4)]
        for p in (Q(0), Q(1, 3), Q(15, 16), Q(1)):
            weights = dense.binomial_weights(3, p)
            self.assertEqual(sum(weights), 1)
            expected = fmpq_mat(2, 2)
            for pattern in product((0, 1), repeat=3):
                j = sum(pattern)
                expected += local[j] * rational(p**j * (1 - p)**(3 - j))
            self.assertEqual(exact_iid(local, p), expected)
        self.assertEqual(dense.binomial_weights(0, Q(1, 2)), (Q(1),))

    def test_exact_noncommuting_marker_conditioning(self):
        # Three local occupancies, two slots per macro, two macros per region.
        local = [fmpq_mat([[1, 1], [0, 1]]), fmpq_mat([[1, 0], [1, 1]]),
                 fmpq_mat([[2, 1], [0, 1]])]
        self.assertNotEqual(local[0] * local[1], local[1] * local[0])
        groups, windows, epochs, activity = 4, 2, 2, Q(2, 3)
        regional = dense.q1.placement(local, epochs=epochs, windows=windows,
            matrix=fmpq_mat, rounding=lambda x: x, maximum_groups=groups)
        identity = fmpq_mat([[1, 0], [0, 1]])
        for q in range(groups + 1):
            p = Q(q, groups)
            iid = exact_iid(local, activity * p)**epochs
            fixed = fmpq_mat(2, 2)
            for j, weight in enumerate(dense.binomial_weights(q, activity)):
                fixed += regional[j] * rational(weight)
            # Enumerate unmarked, marked-zero, marked-nonzero states explicitly.
            event_moment = fmpq_mat(2, 2)
            unconditioned = fmpq_mat(2, 2)
            event_probability = Q(0)
            for pattern in product((0, 1, 2), repeat=groups):
                weight, matrix = Q(1), identity
                for value in pattern:
                    weight *= (1 - p, p * (1 - activity), p * activity)[value]
                for e in range(epochs):
                    count = sum(value == 2 for value in pattern[e*windows:(e+1)*windows])
                    matrix = matrix * local[count]
                unconditioned += matrix * rational(weight)
                if sum(value != 0 for value in pattern) == q:
                    event_probability += weight
                    event_moment += matrix * rational(weight)
            self.assertEqual(unconditioned, iid)
            self.assertEqual(event_moment / rational(event_probability), fixed)
            loss = rational(dense.conditioning_loss_exact(groups, q))
            self.assertTrue(entrywise_le(fixed, iid * loss))
            # Regional products preserve state: no projection to a scalar at seams.
            self.assertTrue(entrywise_le(fixed**3, (iid * loss)**3))
            self.assertEqual((iid * loss)**3, exact_iid(local, activity*p)**6 * loss**3)
            if q == groups:
                self.assertEqual(fixed, iid)

    def test_continuous_state_power_not_reset(self):
        local = [fmpq_mat([[1, 1], [1, 0]]) for _ in range(33)]
        g = dense.geometry(32)
        matrix = exact_iid(local, Q(15, 32))
        power = matrix**(g.regions * g.macros_per_region)
        continuous = sum(power[0, j] for j in range(2))
        reset = sum(matrix[0, j] for j in range(2))**g.regions
        self.assertNotEqual(continuous, reset)

    def test_outward_bound_matches_direct_scalar_formula(self):
        ctx.prec = 128
        g = dense.geometry(32)
        local = [arb_mat([[31]]) for _ in range(33)]
        # Exact positive scalar local operator, avoiding any selected-map census.
        for q in (1, 16, 32):
            upper = dense.occupancy_upper(local, finite_geometry=g, q=q,
                                          tilt=Q(1, 100), beta=Q(5, 2))
            exact_prefactor = (comb(32, q) * Q(5, 2)**q
                * dense.conditioning_loss_exact(32, q)**64 * 31**64)
            expected = (dense.q1.kernel_t64.aq(exact_prefactor)
                        * (dense.q1.kernel_t64.aq(Q(1, 100)) * (g.N // 10)).exp())
            self.assertGreaterEqual(endpoint_fraction(upper),
                                    endpoint_fraction(expected.lower()))
            # Upper endpoints differ only by outward arithmetic error.
            self.assertTrue((upper / expected - 1).abs_upper() < arb(2)**-100)

    def test_options_reject_ambiguous_scope(self):
        invalid = [dict(q_min=0), dict(q_max=8193), dict(precision=True),
                   dict(tilts=['0']), dict(tilts=['.5', '1/2']), dict(tilts='0.5')]
        base = dict(group_count=8192, q_min=33, q_max=None, precision=128,
                    tilts=['.01'], output=None)
        for update in invalid:
            with self.assertRaises(ValueError):
                dense._options(**(base | update))

    def test_arbitrary_marker_conditioning(self):
        ctx.prec = 128
        for groups in (2, 7, 16):
            for q in range(groups + 1):
                for p in (Q(1, 8), Q(1, 2), Q(7, 8)):
                    exact = dense.conditioning_loss_exact(groups, q, p)
                    self.assertEqual(exact * comb(groups, q) * p**q * (1-p)**(groups-q), 1)
                    upper = endpoint_fraction(dense.log_conditioning_loss(
                        groups, q, marker_probability=p).exp())
                    self.assertGreaterEqual(upper, exact)
        for q in (0, 1, 15):
            with self.assertRaises(ValueError):
                dense.conditioning_loss_exact(16, q, 1)

    def test_physical_iid_square_equals_macro_convolution(self):
        # Sixteen physical packet slots, as in the actual selected-map kernel.
        physical = [fmpq_mat([[1, j], [j+1, 2]]) for j in range(17)]
        macro = []
        for count in range(33):
            value = fmpq_mat(2, 2)
            for left in range(max(0, count-16), min(16, count)+1):
                weight = Q(comb(16, left) * comb(16, count-left), comb(32, count))
                value += physical[left] * physical[count-left] * rational(weight)
            macro.append(value)
        for p in (Q(1, 128), Q(15, 32), Q(15, 16)):
            self.assertEqual(exact_iid(macro, p), exact_iid(physical, p)**2)

    def test_fugacity_scalar_coefficient_formula(self):
        ctx.prec = 128
        g, beta, tilt, moment = dense.geometry(32), Q(5, 2), Q(1, 100), Q(3, 5)**64
        for q in (1, 16, 32):
            for p in (Q(1, 8), Q(1, 2), Q(7, 8)):
                logarithm = dense.coefficient_log_bound(finite_geometry=g, q=q, tilt=tilt,
                    beta=beta, marker_probability=p, log_moment=dense.q1.kernel_t64.aq(moment).log())
                exact = comb(32, q) * beta**q * dense.conditioning_loss_exact(32, q, p)**64 * moment
                expected = dense.q1.kernel_t64.aq(exact) * (dense.q1.kernel_t64.aq(tilt) * (g.N//10)).exp()
                self.assertTrue((logarithm.exp() / expected - 1).abs_upper() < arb(2)**-100)
        endpoint = dense.coefficient_log_bound(finite_geometry=g, q=32, tilt=tilt,
            beta=beta, marker_probability=1, log_moment=dense.q1.kernel_t64.aq(moment).log())
        expected = (dense.q1.kernel_t64.aq(moment * beta**32)
                    * (dense.q1.kernel_t64.aq(tilt) * (g.N//10)).exp())
        self.assertTrue((endpoint.exp() / expected - 1).abs_upper() < arb(2)**-100)

    def test_selector_keeps_complete_rows(self):
        ctx.prec = 128
        local = [arb_mat([[1, 2], [1, 1]]) for _ in range(17)]
        options = ((arb(1), arb(2)), (arb(2), arb(1)))
        candidates = dict(local=local, rows=tuple(options for _ in range(17)), precision=ctx.prec,
            arrays=dense.np.array([[[1., 2.], [1., 1.]]] * 17),
            row_arrays=tuple(dense.np.array([[1., 2.], [2., 1.]]) for _ in range(17)))
        with patch.object(dense.q1.kernel_t64.kernel_birth_density, 'continuation',
                          return_value=dense.np.array([1., 2.])):
            actual, indices = dense.select_physical(candidates, Q(1, 8))
        self.assertEqual(indices, (1,) * 17)
        self.assertEqual(actual[0, 0], arb(2))
        self.assertEqual(actual[0, 1], arb(1))
        self.assertNotEqual(actual[0, 0] + actual[0, 1], arb(2))  # invalid entrywise minimum

    def test_fugacity_run_caches_once_per_tilt_and_records_scope(self):
        from contextlib import redirect_stdout
        from io import StringIO
        data = dict(physical_step_bits=64, bits=16, macro_windows=32,
                    birth_density='capped', map_sha256='toy')
        record = dict(map_sha256='toy')
        with (patch.object(dense.q1.kernel_t64, 'authenticate'),
              patch.object(dense, 'physical_candidates', return_value='cached') as prepare,
              patch.object(dense, 'select_physical', return_value=(arb_mat([[1, 0], [0, 1]]), (0,)*17)) as select,
              redirect_stdout(StringIO())):
            result = dense.run_fugacity(group_count=32, q_min=30, q_max=32, precision=128,
                tilts=('1/100', '1/50'), marker_probabilities=('1/2', '3/4', '1'),
                distance='19/200', data=data, map_record=record)
        self.assertEqual(prepare.call_count, 2)
        self.assertEqual(select.call_count, 6)
        self.assertEqual(result['schema'], dense.FUGACITY_SCHEMA)
        self.assertEqual(result['distance'], '19/200')
        self.assertEqual(result['threshold'], 8192 * 19 // 200)
        self.assertEqual(result['occupancy_covered'], [30, 32])
        self.assertEqual(set(result['occupancy_uppers']), {'30', '31', '32'})
        self.assertFalse(result['whole_code_certificate'])
        self.assertTrue(result['evaluated_every_integer_occupancy'])


if __name__ == '__main__':
    unittest.main()
