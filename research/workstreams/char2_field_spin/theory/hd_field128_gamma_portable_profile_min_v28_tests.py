"""Tiny array/validation and mocked orchestration tests; no coefficient jobs."""
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

import hd_field128_gamma_portable_profile_min_v28 as screen
import hd_field128_gamma_transfer_diagnostic as transfer


class Budget:
    def check(self):
        pass


class Counts:
    def __init__(self, values):
        self.values = np.array(values, dtype=float)
        self.calls = 0

    def log_counts(self, first, second, *, budget):
        self.calls += 1
        with np.errstate(divide='ignore'):
            return np.log2(self.values[first, second])


def proposal():
    return dict(active_groups=548, construction=screen.CONSTRUCTION, state_dimension=30,
        target_distance='2/5', field_order=1 << 128, inherited_bound=-99999,
        witness=dict(z='1/8', alpha='1/128', grid=128, state_dimension=30,
                     markers=[['1/2', '3/4'], ['2', '3']], inherited_bound=-99999))


def combine(matrices, prefixes, counts, *, repetitions=4, chunk=1):
    with np.errstate(divide='ignore'):
        logs = [np.log2(m) for m in matrices]
    return screen.combine_profiles(logs, prefixes, np.arange(len(matrices[0])),
        counts, transfer, Budget(), repetitions=repetitions, chunk=chunk)


def oracle(matrices, prefixes, counts, repetitions):
    side = len(matrices[0])
    terms = np.zeros((len(matrices), side, side))
    for j, (regions, prefix) in enumerate(zip(matrices, prefixes)):
        for a in range(side):
            for b in range(side):
                product = np.linalg.matrix_power(regions[a] @ regions[b], repetitions//2)
                terms[j, a, b] = counts[a, b]*2.**prefix*sum(product[0])
    return terms


class Tests(unittest.TestCase):
    def test_independent_noncommuting_transfer_and_disjoint_classes(self):
        matrices = [np.array([[[1., 2.], [0., 1.]], [[1., 0.], [3., 1.]]]),
                    np.array([[[2., 0.], [1., 1.]], [[1., 1.], [0., 2.]]])]
        prefixes = [-3., 2.]
        counts = Counts([[3, 5], [7, 0]])
        actual = combine(matrices, prefixes, counts, repetitions=8)
        terms = oracle(matrices, prefixes, counts.values, 8)
        np.testing.assert_allclose(actual['single_log2_upper_screens'], np.log2(terms.sum(axis=(1, 2))))
        self.assertAlmostEqual(actual['combined_log2_upper_screen'], log2(terms.min(axis=0).sum()), places=12)
        self.assertEqual(actual['profile_pairs_evaluated'], 4)
        self.assertEqual(actual['positive_count_profiles'], 3)
        self.assertEqual(actual['zero_count_profiles'], 1)
        self.assertEqual(sum(actual['witness_winning_profile_counts']), 3)
        self.assertEqual(counts.calls, 2)  # Shared count caps once per chunk, not per witness.
        self.assertTrue(actual['all_profiles_retained'])
        self.assertFalse(actual['entrywise_cross_witness_minimum'])

    def test_distinct_prefixes_must_be_inside_class_minimum(self):
        matrices = [np.array([np.eye(2)*2]), np.array([np.eye(2)])]
        result = combine(matrices, [-4., 0.], Counts([[1]]), repetitions=2)
        self.assertAlmostEqual(result['combined_log2_upper_screen'], -2.)
        self.assertEqual(result['witness_winning_profile_counts'], [1, 0])
        self.assertNotEqual(result['combined_log2_upper_screen'], -4.)

    def test_never_mix_entries_or_even_odd_matrices(self):
        matrices = [np.array([[[1., 8.], [8., 1.]]]),
                    np.array([[[8., 1.], [1., 8.]]])]
        good = combine(matrices, [0., 0.], Counts([[1]]), repetitions=2)
        bad = combine([np.minimum(*matrices)], [0.], Counts([[1]]), repetitions=2)
        self.assertAlmostEqual(good['combined_log2_upper_screen'], log2(81.))
        self.assertAlmostEqual(bad['combined_log2_upper_screen'], 2.)
        self.assertGreater(good['combined_log2_upper_screen'], bad['combined_log2_upper_screen'])
        crossed = [np.array([np.eye(2), np.eye(2)*8]), np.array([np.eye(2)*8, np.eye(2)])]
        complete = combine(crossed, [0., 0.], Counts([[1, 1], [1, 1]]), repetitions=2)
        expected = oracle(crossed, [0., 0.], np.ones((2, 2)), 2).min(axis=0)
        self.assertAlmostEqual(complete['combined_log2_upper_screen'], log2(expected.sum()))

    def test_single_recovery_chunk_invariance_and_first_ties(self):
        regions = np.array([[[1., .5], [.25, 1.]], [[.5, .25], [1., .5]]])
        one = combine([regions], [-7.], Counts([[2, 3], [4, 5]]))
        two = combine([regions, regions], [-7., -7.], Counts([[2, 3], [4, 5]]), chunk=2)
        self.assertAlmostEqual(one['combined_log2_upper_screen'], one['single_log2_upper_screens'][0])
        self.assertAlmostEqual(one['combined_log2_upper_screen'], two['combined_log2_upper_screen'])
        self.assertEqual(two['witness_winning_profile_counts'], [4, 0])

    def test_literal_projection_and_rejections(self):
        saved = proposal()
        projected, key = screen.project_literals(saved, 548)
        self.assertEqual(key, 'markers')
        self.assertNotIn('inherited_bound', projected)
        projected['markers'][0][0] = '7'
        self.assertEqual(saved['witness']['markers'][0][0], '1/2')
        alternate = proposal()
        alternate['witness']['exact_positive_marker_proposals'] = [['7', '8']]
        self.assertEqual(screen.project_literals(alternate, 548)[0]['markers'], [['7', '8']])
        for mutate in (
                lambda p: p.update(active_groups=547), lambda p: p.update(construction='other'),
                lambda p: p.update(target_distance='9/20'), lambda p: p.update(field_order=256),
                lambda p: p['witness'].update(grid=64), lambda p: p['witness'].update(state_dimension=31),
                lambda p: p['witness'].update(z='2/16'), lambda p: p['witness'].update(alpha='0'),
                lambda p: p['witness'].update(z='1'), lambda p: p['witness'].update(z=.125),
                lambda p: p['witness'].update(markers=[['0', '1']]),
                lambda p: p['witness'].update(markers=[['1', '2'], ['1', '2']])):
            bad = proposal(); mutate(bad)
            with self.assertRaises(ValueError):
                screen.project_literals(bad, 548)
        for name in ('../x.json', '/x.json', 'C:/x.json', 'a\\x.json', 'a//x.json', 'a/../x.json'):
            with self.assertRaises(ValueError):
                screen.proposal_name(name)
        for target in (1, 4096, True, 548.):
            with self.assertRaises(ValueError):
                screen.validate_target(target)

    def test_explicit_pin_checks_and_no_ancestry_reads(self):
        with TemporaryDirectory() as folder, patch.object(screen, 'HERE', Path(folder)):
            raw = json.dumps(dict(proposal(), source={'path': 'missing/ancestor.json', 'sha256': 'a'*64})).encode()
            first, second = Path(folder)/'first.json', Path(folder)/'second.json'
            first.write_bytes(raw); second.write_bytes(raw)
            digest = sha256(raw).hexdigest()
            specs = [('first.json', digest), ('second.json', digest)]
            loaded = screen.load_proposals(specs, 548)
            self.assertEqual(len(loaded), 2)
            self.assertNotIn('inherited_bound', loaded[0])
            for bad in (specs[:1], specs*3, [specs[0], specs[0]], [specs[0], ('second.json', 'a'*64)]):
                with self.assertRaises(ValueError):
                    screen.load_proposals(bad, 548)
            first.write_bytes(raw+b' ')
            with self.assertRaises(ValueError):
                screen.load_proposals(specs, 548)

    def test_nonfinite_arrays_and_zero_union_fail_closed(self):
        for matrix in (np.full((1, 2, 2), np.nan), np.full((1, 2, 2), np.inf)):
            with self.assertRaises(ValueError):
                screen.combine_profiles([matrix], [0.], np.array([0]), Counts([[1]]), transfer, Budget())
        with self.assertRaises(ArithmeticError):
            combine([np.zeros((1, 2, 2))], [0.], Counts([[1]]))

    def test_actual_tiny_midpoint_single_witness_control(self):
        import hd_field128_gamma_midpoint_counts as mode
        counts = mode.MidpointCounts(4, 4, max_singletons=4, budget=Budget())
        regions = np.log2(np.array([[[1., .5], [.25, 1.]]]*3))
        actual = screen.combine_profiles([regions], [-3.], np.arange(3), counts,
            transfer, Budget(), repetitions=8)
        old = mode.support_sum(4, 4, regions, counts=counts, repetitions=8, budget=Budget())
        self.assertAlmostEqual(actual['combined_log2_upper_screen'], old['log2_support_sum_screen']-3.)
        self.assertEqual(actual['profile_pairs_evaluated'], old['profile_pairs_evaluated'])


if __name__ == '__main__':
    unittest.main()
