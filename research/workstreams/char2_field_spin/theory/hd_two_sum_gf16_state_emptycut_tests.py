"""Exact small labelled path/mask controls; no regional bank or gate."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_state_emptycut as target


def toy_curves(c=.25, maximum=3):
    slopes = np.log2([.8, 1.2, .7, .9])
    intercepts = np.log2([c or 1., .3, .4, .6])
    curves = intercepts[:, None]+slopes[:, None]*np.arange(maximum+1)
    curves[0, 0] = 0.
    if c == 0.:
        curves[0, 1:] = -np.inf
    return curves, slopes


def labelled_column(weights, types, envelope, regions):
    """Independently sum every mask tuple and ordered regional matrix word."""
    result = np.zeros((2, 2))
    b = 2**envelope['slopes']
    a = 2**envelope['intercepts']
    c = 2**envelope['log_c']
    for masks in product(range(1 << regions), repeat=len(types)):
        mass = np.prod([2**weights[arity][mask.bit_count()] for arity, mask in zip(types, masks)])
        word = np.eye(2)
        for region in range(regions):
            occupancy = sum((mask >> region) & 1 for mask in masks)
            matrix = (a*b**occupancy).reshape(2, 2).copy()
            matrix[0, 0] = b[0]**occupancy*(c+(1-c)*(occupancy == 0))
            word = word@matrix
        result += mass*word
    return result


def expanded_column(prepared, types):
    result = np.zeros((2, 2))
    for index, (initial, final, _, multiplicity) in enumerate(prepared['classes']):
        value = prepared['path_a'][index]+sum(prepared['typed'][index, arity-1] for arity in types)
        result[initial, final] += multiplicity*2**value
    return result


class StateEmptyCutTests(unittest.TestCase):
    def test_complete_envelope_and_exact_empty_endpoint(self):
        curves, slopes = toy_curves()
        curves[0, 1] = -np.inf
        envelope = target.empty_envelope(curves, slopes)
        self.assertAlmostEqual(envelope['log_c'], -2.)
        self.assertAlmostEqual(2**envelope['log_one_minus_c'], .75)
        self.assertEqual(envelope['slack'][0, 0], 0.)
        self.assertTrue(np.all(envelope['slack'] >= -1e-14))
        for occupancy in range(4):
            value = 2**(slopes[0]*occupancy)*(.25+.75*(occupancy == 0))
            self.assertGreaterEqual(value+1e-15, 2**curves[0, occupancy])

    def test_small_labelled_masks_paths_and_cut_multiplicities(self):
        regions = 3
        curves, slopes = toy_curves()
        weights = {1: np.log2([1., 2., 3., 4.]), 2: np.log2([2., 1., 5., 3.])}
        prepared = target.prepared_envelope(curves, slopes, weights, regions=regions)
        envelope = target.empty_envelope(curves, slopes)
        for types in ((1,), (2,), (1, 1), (1, 2), (2, 1), (2, 2)):
            np.testing.assert_allclose(expanded_column(prepared, types),
                labelled_column(weights, types, envelope, regions), rtol=5e-14)
        # For every original state path, its cuts sum to probability one.
        original = target.state.path_classes(regions)
        for initial, final, counts, multiplicity in original:
            total = 0.
            for index, row in enumerate(prepared['classes']):
                if row[:3] == (initial, final, counts):
                    allowed = prepared['allowed_counts'][index, 0]
                    total += row[3]*(.25**allowed)*(.75**(counts[0]-allowed))
            self.assertAlmostEqual(total, multiplicity)

    def test_fixed_type_holder_minkowski_with_exact_arity_restriction(self):
        curves, slopes = toy_curves()
        weights = {1: np.log2([1., 2., 3., 4.]), 2: np.log2([2., 1., 5., 3.])}
        prepared = target.prepared_envelope(curves, slopes, weights, regions=3)
        envelope = target.empty_envelope(curves, slopes)
        columns, groups, bound = 3, 2, 3
        exact_union = 0.
        for types in product((1, 2), repeat=groups):
            if sum(types) <= bound:
                column = labelled_column(weights, types, envelope, 3)
                exact_union += np.linalg.matrix_power(column, columns)[0].sum()
        norms = target.exact.restricted_type_norms(prepared['typed'], {groups: (0, bound)}, columns=columns)[groups]
        kernel = 2**target.exact.kernel_from_restricted_norms(prepared['classes'], prepared['path_a'], norms)
        upper = np.linalg.matrix_power(kernel, columns)[0].sum()
        self.assertGreaterEqual(upper, exact_union*(1-1e-13))

    def test_c_endpoints_reduce_without_zero_times_infinity(self):
        weights = {1: np.log2([1., 2., 3., 4.])}
        for c in (0., 1.):
            curves, slopes = toy_curves(c)
            prepared = target.prepared_envelope(curves, slopes, weights, regions=3)
            self.assertEqual(prepared['state_cut_classes'], prepared['original_state_path_classes'])
            self.assertTrue(np.all(np.isfinite(prepared['path_a'])))
            self.assertFalse(np.any(np.isnan(prepared['typed'])))
            np.testing.assert_allclose(expanded_column(prepared, (1, 1)),
                labelled_column(weights, (1, 1), target.empty_envelope(curves, slopes), 3), rtol=5e-14)
            for row, allowed in zip(prepared['classes'], prepared['allowed_counts']):
                self.assertEqual(allowed[0], row[2][0] if c == 1 else 0)

    def test_fractional_cardinality_adapter_matches_pinned_formula(self):
        dense = {a: np.linspace(-2*a, a/4, 17) for a in range(3, 17)}
        markers = np.array([np.linspace(-3., 1., 16), np.linspace(.2, -.7, 16)])
        floors = np.arange(16)/13
        elementary = target.state.elementary_logs(markers)
        for alpha in (F(1), F(1, 5)):
            weights = target.fractional_cardinality_weights(alpha, dense, floors)
            actual = np.stack([target.state.logsum(elementary+weights[a], axis=1) for a in range(1, 17)], axis=1)
            expected = target.fractional.fractional_type_weights(markers, dense, floors, alpha)
            np.testing.assert_allclose(actual, expected, atol=2e-13, rtol=0)

    def test_forbidden_coordinates_have_exact_zero_higher_coefficients(self):
        curves, slopes = toy_curves(0.)
        prepared = target.prepared_envelope(curves, slopes, {1: np.zeros(4)}, regions=3)
        found_empty = False
        for counts, elementary in zip(prepared['allowed_counts'], prepared['elementary']):
            self.assertTrue(np.all(np.isneginf(elementary[sum(counts)+1:])))
            if not sum(counts):
                found_empty = True
                np.testing.assert_array_equal(elementary, [0., -np.inf, -np.inf, -np.inf])
        self.assertTrue(found_empty)

    def test_prefix_scope_and_no_second_fractional_power(self):
        prepared = dict(regions=16, typed=np.zeros((1, 16)), classes=((0, 0, (16, 0, 0, 0), 1),),
                        path_a=np.array([0.]), state_cut_classes=1, log_c=-2., log_one_minus_c=log2(.75))
        identity = np.array([[0., -np.inf], [-np.inf, 0.]])
        with patch.object(target.state, 'kernel_from_types', return_value=identity) as kernel:
            row = target.evaluate_prepared(prepared, 18, F(691, 700), F(1, 5))
        self.assertEqual(kernel.call_args.args[-1], 18)
        expected = .2*(-target.fractional.LIMIT*log2(691/700)-log2(255))+log2(comb(8192, 18))
        self.assertAlmostEqual(row['failure_log2_screen'], expected)
        self.assertFalse(row['cut_coefficient_powered_again'])
        self.assertFalse(row['directed'])
        self.assertFalse(row['continuous_weight_band_covered'])
        with patch.object(target.exact, 'evaluate_prepared', return_value={'directed': False}) as restricted:
            target.evaluate_prepared(prepared, 18, F(691, 700), F(1, 5), cutoff=247)
        self.assertEqual(restricted.call_args.args[1:3], (18, 247))

    def test_invalid_envelope_and_laws_fail_closed(self):
        curves, slopes = toy_curves()
        bad = curves.copy()
        bad[0, 0] = -.1
        with self.assertRaises(ValueError):
            target.empty_envelope(bad, slopes)
        bad = curves.copy()
        bad[0, 2] = 1.
        with self.assertRaises(ValueError):
            target.empty_envelope(bad, slopes)
        for weights in ({}, {2: np.zeros(4)}, {1: np.zeros(3)}, {1: [0., 1., np.nan, 2.]}):
            with self.assertRaises(ValueError):
                target.prepared_envelope(curves, slopes, weights, regions=3)


if __name__ == '__main__':
    unittest.main()
