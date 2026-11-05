"""Pure array and mocked optimizer tests; no coefficient banks or file writes."""
from copy import deepcopy
from fractions import Fraction as F
from types import SimpleNamespace
import unittest

import numpy as np

import hd_field128_gamma_portable_residual_search_v28 as search
import hd_field128_gamma_transfer_diagnostic as transfer


class Budget:
    def check(self):
        pass


class Counts:
    def __init__(self, values):
        self.values = np.array(values)
        self.calls = 0
    def log_counts(self, a, b, *, budget):
        self.calls += 1
        with np.errstate(divide='ignore'):
            return np.log2(self.values[a, b])


def witness(markers=None):
    return dict(z='1/8', alpha='1/128', grid=128, state_dimension=30,
                markers=markers or [['1', '2'], ['3', '4']])


def row(score):
    return dict(combined_log2_upper_screen=score, candidate_single_log2_upper_screen=score+3,
        largest_profile=dict(even_full=1, odd_full=0, log2_term=score-1),
        candidate_winning_profile_count=2, profile_pairs_evaluated=4, all_profiles_retained=True)


class Tests(unittest.TestCase):
    def test_marker_union_is_fixed_exact_and_deduplicated(self):
        proposals = [dict(witness=witness()), dict(witness=witness([['3', '4'], ['5', '6']]))]
        before = deepcopy(proposals)
        initial, policy = search.candidate_witness(proposals)
        self.assertEqual(initial['markers'], [['1', '2'], ['3', '4'], ['5', '6']])
        self.assertEqual(policy['policy'], 'deduplicated union')
        self.assertEqual((initial['z'], initial['alpha']), ('1/8', '1/128'))
        initial['markers'][0][0] = '9'
        self.assertEqual(proposals, before)
        big = [dict(witness=witness([[str(i), '1'] for i in range(1, 76)])),
               dict(witness=witness([[str(i), '1'] for i in range(76, 130)]))]
        fallback, policy = search.candidate_witness(big)
        self.assertEqual(fallback, big[0]['witness'])
        self.assertEqual(policy['deduplicated_union_count'], 129)
        self.assertIn('exceeds 128', policy['policy'])
        with self.assertRaises(ValueError):
            search.candidate_witness(proposals[:1])
        proposals[0]['witness']['z'] = '1/3'
        with self.assertRaises(ValueError):
            search.candidate_witness(proposals)

    def test_residual_equals_reviewed_three_witness_screen(self):
        full, budget = np.arange(2), Budget()
        count = Counts([[2, 3], [7, 0]])
        logs = search.count_array(count, full, budget, chunk=1)
        self.assertEqual(count.calls, 2)
        matrices = [np.array([[[1., 2.], [0., 1.]], [[1., 0.], [3., 1.]]]),
                    np.array([[[2., 0.], [1., 1.]], [[1., 1.], [0., 2.]]]),
                    np.array([[[.5, 1.], [2., .5]], [[1., .5], [.25, 1.]]])]
        with np.errstate(divide='ignore'):
            regions = [np.log2(m) for m in matrices]
        prefixes = [-7., 2., -2.]
        terms = [search.class_terms(r, p, logs, transfer, budget, repetitions=8, chunk=1)
                 for r, p in zip(regions, prefixes)]
        baseline = np.minimum(terms[0], terms[1])
        saved = baseline.copy()
        result = search.residual_score(baseline, terms[2], logs, full)
        control = search.profiles.combine_profiles(regions, prefixes, full, count, transfer,
                                                   budget, repetitions=8, chunk=1)
        self.assertAlmostEqual(result['combined_log2_upper_screen'], control['combined_log2_upper_screen'])
        self.assertAlmostEqual(result['candidate_single_log2_upper_screen'], control['single_log2_upper_screens'][2])
        self.assertEqual(result['profile_pairs_evaluated'], 4)
        np.testing.assert_array_equal(saved, baseline)
        # Counts were prepared once; each repeated transfer uses this same array.
        self.assertEqual(count.calls, 4)  # Two additional calls belong to the independent control.

    def test_complete_transfer_and_own_prefix_precede_minimum(self):
        logs, full = np.zeros((1, 1)), np.array([0])
        first = np.log2(np.array([[[1., 8.], [8., 1.]]]))
        second = np.log2(np.array([[[8., 1.], [1., 8.]]]))
        a = search.class_terms(first, 0., logs, transfer, Budget(), repetitions=2)
        b = search.class_terms(second, -2., logs, transfer, Budget(), repetitions=2)
        score = search.residual_score(a, b, logs, full)
        self.assertAlmostEqual(score['combined_log2_upper_screen'], np.log2(81.)-2)
        self.assertEqual(score['candidate_winning_profile_count'], 1)
        bad = search.class_terms(np.minimum(first, second), -2., logs, transfer, Budget(), repetitions=2)
        self.assertLess(search.log_union(bad), score['combined_log2_upper_screen'])
        for prefix in (np.inf, np.nan):
            with self.assertRaises(ValueError):
                search.class_terms(first, prefix, logs, transfer, Budget())
        with self.assertRaises(ArithmeticError):
            search.log_union(np.full((1, 1), -np.inf))

    def test_mock_search_keeps_exact_literals_best_and_duplicate_cache(self):
        observed = []
        def evaluate(w):
            observed.append(deepcopy(w))
            return row(-len(observed))
        def minimize(objective, start, **kwargs):
            self.assertEqual(kwargs['options']['maxfev'], 4)
            self.assertEqual(objective(start), -1)
            self.assertEqual(objective(start), -1)
            self.assertEqual(objective(start+[.5, 0]), -2)
            return SimpleNamespace(success=False, message='mock limit', nfev=3)
        initial = witness()
        result = search.search_candidates(initial, evaluate, Budget(), minimize, max_banks=4)
        self.assertEqual(result['attempted_candidate_banks'], 2)
        self.assertEqual(result['complete_candidate_banks'], 2)
        self.assertEqual(result['best']['combined_log2_upper_screen'], -2)
        self.assertEqual(result['best']['evaluation_index'], 1)
        self.assertEqual(observed[0], initial)
        self.assertEqual(initial, witness())
        for w in observed:
            self.assertEqual(w['markers'], initial['markers'])
            self.assertEqual(F(float(F(w['z']))), F(w['z']))
            self.assertEqual(F(float(F(w['alpha']))), F(w['alpha']))

    def test_interruption_never_promotes_unfinished_candidate(self):
        observed = []
        def evaluate(w):
            observed.append(w)
            if len(observed) == 2:
                raise TimeoutError('toy stop')
            return row(7.)
        def minimize(objective, start, **kwargs):
            objective(start)
            objective(start+[.5, 0])
            self.fail('Timeout must escape the optimizer')
        result = search.search_candidates(witness(), evaluate, Budget(), minimize)
        self.assertEqual(result['attempted_candidate_banks'], 2)
        self.assertEqual(result['complete_candidate_banks'], 1)
        self.assertEqual(result['best']['combined_log2_upper_screen'], 7.)
        self.assertEqual(result['failures'][0]['evaluation_index'], 1)
        self.assertFalse(result['failures'][0]['complete_profile_sum'])
        self.assertEqual(result['optimizer']['resource_stop'], 'TimeoutError')

    def test_bank_cap_and_no_complete_candidate(self):
        def minimize(objective, start, **kwargs):
            objective(start)
            objective(start+[.5, 0])
        result = search.search_candidates(witness(), lambda w: row(3.), Budget(), minimize, max_banks=1)
        self.assertEqual(result['attempted_candidate_banks'], 1)
        self.assertEqual(result['complete_candidate_banks'], 1)
        self.assertEqual(result['optimizer']['resource_stop'], 'TimeoutError')
        def failed(w):
            raise TimeoutError('before first full profile sum')
        result = search.search_candidates(witness(), failed, Budget(), minimize)
        self.assertIsNone(result['best'])
        self.assertEqual(result['history'], [])
        self.assertEqual(result['complete_candidate_banks'], 0)


if __name__ == '__main__':
    unittest.main()
