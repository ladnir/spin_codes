"""Small positive labelled-coverage tests; no saved bank or provider run."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_state_cover as target


def labelled_transitions(weights, required, required_marks, optional):
    first, second = required
    states = [(i, j) for i in range(first+1) for j in range(second+1)]
    result = np.zeros((len(weights), len(states), len(states)))
    marks = [required_marks[0]]*first+[required_marks[1]]*second+list(optional)
    for source, (covered1, covered2) in enumerate(states):
        old = ((1 << covered1)-1) | (((1 << covered2)-1) << first)
        for mask in range(1 << len(marks)):
            union = old | mask
            next1 = (union & ((1 << first)-1)).bit_count()
            next2 = ((union >> first) & ((1 << second)-1)).bit_count()
            destination = next1*(second+1)+next2
            mark = sum(value for index, value in enumerate(marks) if mask & (1 << index))
            result[:, source, destination] += 2**(weights[:, mask.bit_count()]+mark)
    return result


class StateCoverTests(unittest.TestCase):
    def test_transitions_match_every_labelled_mask_and_covered_set(self):
        for required in ((0, 0), (1, 0), (0, 2), (1, 2), (2, 2)):
            optional = np.array([-.7])
            length = sum(required)+len(optional)
            weights = np.log2(np.array([np.arange(1, length+2), np.arange(length+2, 1, -1)]))
            required_marks = [.3, -.4]
            actual = target.type_transition_logs(weights, required, required_marks, optional)
            expected = labelled_transitions(weights, required, required_marks, optional)
            np.testing.assert_allclose(2**actual, expected, rtol=3e-14, atol=0)

    def test_coverage_history_minkowski_preserves_fixed_type_vectors(self):
        weights = np.log2(np.array([[1., 2., 4., 3.], [3., 1., 2., 5.]]))
        transitions = target.type_transition_logs(weights, (1, 1), [-.2, .3], [-.5])
        for groups in (1, 2, 3):
            for columns in (1, 2, 4):
                exact = 0.
                for types in product(range(2), repeat=groups):
                    run = np.eye(4)
                    for a in types:
                        run = run@(2**transitions[a])
                    exact += run[0, -1]**columns
                expected_norm = exact**(1/columns)
                actual_norm = 2**target.cover_norm_from_transitions(transitions, groups, columns=columns)
                self.assertGreaterEqual(actual_norm, expected_norm*(1-1e-13))
                if columns == 1 or groups == 1:
                    self.assertAlmostEqual(actual_norm/expected_norm, 1., places=12)

    def test_all_empty_masks_cannot_cross_an_off_diagonal_step(self):
        weights = np.array([[0., -np.inf, -np.inf, -np.inf],
                            [2., -np.inf, -np.inf, -np.inf]])
        for required in ((1, 0), (0, 1), (1, 1)):
            transitions = target.type_transition_logs(weights, required, [0., 0.], [])
            for groups in (1, 3):
                self.assertTrue(np.isneginf(target.cover_norm_from_transitions(transitions, groups)))

    def test_no_required_coordinates_equals_original_type_norm(self):
        weights = np.log2(np.array([[1., 2., 4., 3.], [3., 1., 2., 5.]]))
        optional = [-.2, .3, -.5]
        transitions = target.type_transition_logs(weights, (0, 0), [0., 0.], optional)
        elementary = target.state.elementary_logs(optional)
        typed = target.state.logsum(weights+elementary, axis=1)
        for groups in (1, 2, 4):
            for columns in (1, 4, 16):
                expected = groups/columns*target.state.logsum(columns*typed)
                self.assertAlmostEqual(target.cover_norm_from_transitions(
                    transitions, groups, columns=columns), expected, places=12)

    def test_positive_tiny_terms_use_log_fallback_without_false_zero(self):
        weights = np.array([[0., -2000.], [-np.inf, -np.inf]])
        transitions = target.type_transition_logs(weights, (1, 0), [0., 0.], [])
        self.assertEqual(transitions[0, 0, 1], -2000.)
        self.assertTrue(np.all(np.isneginf(transitions[1])))
        self.assertEqual(target.cover_norm_from_transitions(transitions, 1), -2000.)

    def test_pathwise_min_and_cache_by_allowed_counts(self):
        prepared = dict(off_diagonal_empty_support_exact=True, regions=16,
            typed=np.array([[2., 3.], [2., 3.], [1., 2.]]),
            allowed_counts=np.array([[0, 1, 1, 0], [0, 1, 1, 0], [2, 0, 0, 0]]),
            coverage_type_weights=np.zeros((2, 17)), slopes=np.zeros(4))
        old = 2/16*target.state.logsum(16*prepared['typed'], axis=1)
        with patch.object(target, 'type_transition_logs', return_value=np.zeros((2, 4, 4))) as build, \
                patch.object(target, 'cover_norm_from_transitions', return_value=old[0]-1.):
            result, diagnostic = target.covered_path_norms(prepared, 2)
        self.assertEqual(build.call_count, 1)
        np.testing.assert_allclose(result, [old[0]-1, old[1]-1, old[2]])
        self.assertEqual(diagnostic['cover_matrices_built'], 1)
        self.assertEqual(diagnostic['unique_allowed_count_classes'], 2)
        with patch.object(target, 'type_transition_logs', return_value=np.zeros((2, 4, 4))), \
                patch.object(target, 'cover_norm_from_transitions', return_value=old[0]+1.):
            result, _ = target.covered_path_norms(prepared, 2)
        np.testing.assert_array_equal(result, old)

    def test_scope_and_prefix_keep_source_count_outside_alpha(self):
        prepared = dict(typed=np.zeros((1, 16)), classes=((0, 0, (16, 0, 0, 0), 1),), path_a=np.array([0.]))
        identity = np.array([[0., -np.inf], [-np.inf, 0.]])
        with patch.object(target, 'covered_path_norms', return_value=(np.array([0.]), {})), \
                patch.object(target.exact, 'kernel_from_restricted_norms', return_value=identity):
            row = target.evaluate_prepared(prepared, 18, F(691, 700), F(1, 5))
        expected = .2*(-target.emptycut.fractional.LIMIT*log2(691/700)-log2(255))+log2(comb(8192, 18))
        self.assertAlmostEqual(row['failure_log2_screen'], expected)
        self.assertFalse(row['inclusion_exclusion_used'])
        self.assertFalse(row['independently_resampled_types'])
        self.assertFalse(row['directed'])
        self.assertFalse(row['continuous_weight_band_covered'])

    def test_guards_and_deadline(self):
        for required in ((-1, 2), (9, 8), (True, 1)):
            with self.assertRaises(ValueError):
                target.coverage_geometry(*required)
        def stop():
            raise TimeoutError('tiny fixture')
        with self.assertRaises(TimeoutError):
            target.cover_norm_from_transitions(np.zeros((2, 4, 4)), 2, check_time=stop)
        curves = np.array([[0., -2.], [-np.inf, -1.], [0., -3.], [-1., -1.]])
        with self.assertRaises(ValueError):
            target.prepared_envelope(curves, np.zeros(4), {a: np.zeros(17) for a in range(1, 17)})


if __name__ == '__main__':
    unittest.main()
