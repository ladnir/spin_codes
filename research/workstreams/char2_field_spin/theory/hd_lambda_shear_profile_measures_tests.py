"""Exact small subset checks and synthetic-bank tests; no real provider run."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import log2
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_shear_profile_measures as target


PINS = {'test-only.py': 'mock'}


def toy_bank(maximum=18):
    values = {}
    for m, f in target.full_keys(maximum):
        values[m, f] = np.array([[0., -np.inf], [-np.inf, -5.]]) if m == 0 else np.array(
            [[-2.-1.5*m, -.5-.25*m], [-10.-.3*m, -5.]])-f/128
    return values


def metadata():
    return dict(mark='9927/10000', maximum=18, reference_cone_unchanged=True,
        nonempty_packet_caps_unchanged=True, exact_cap='1/2',
        exact_cap_sha256=sha256(b'1/2').hexdigest())


class ProfileMeasureTests(unittest.TestCase):
    def test_feasibility_exact_for_every_small_source_subset(self):
        for r in range(1, 8):
            for full in range(r//2+1):
                source = tuple(range(2*full))+tuple(2*full+2*k for k in range(r-2*full))
                for m in range(r+1):
                    actual = set()
                    for subset in combinations(source, m):
                        active = set(subset)
                        actual.add(sum(2*i in active and 2*i+1 in active for i in range(full)))
                    self.assertEqual(actual, set(target.feasible_full_counts(r, m, full)))
        self.assertEqual(target.feasible_full_counts(18, 18, 9), (9,))
        self.assertEqual(target.feasible_full_counts(18, 17, 9), (8,))

    def test_full100_roundtrip_and_exact_null_zeros(self):
        regions = toy_bank()
        data = target.make_bank(regions, metadata(), PINS)
        encoded = json.loads(json.dumps(data, allow_nan=False))
        loaded = target.validate_bank(encoded, PINS)
        self.assertEqual(len(loaded), 100)
        self.assertIsNone(encoded['rows'][0]['log2_matrix'][0][1])
        for key in regions: self.assertTrue(np.array_equal(regions[key], loaded[key]))
        for changed in (dict(complete=False), dict(reviewed_source_sha256={}), dict(maximum=17)):
            bad = deepcopy(data)
            bad.update(changed)
            with self.assertRaises(ValueError): target.validate_bank(bad, PINS)
        bad = deepcopy(data)
        bad['rows'][-1] = bad['rows'][0]
        with self.assertRaises(ValueError): target.validate_bank(bad, PINS)
        bad = deepcopy(data)
        bad['rows'].pop()
        with self.assertRaises(ValueError): target.validate_bank(bad, PINS)

    def test_profile_union_uses_both_lower_bounds_and_retains_ordered_count_interface(self):
        regions = toy_bank(6)
        for key in regions: regions[key] = np.full((2, 2), float(key[1]))
        forward, feasible = target.profile_curves(regions, 6, (1, 3))
        reverse, _ = target.profile_curves(regions, 6, (3, 1))
        self.assertTrue(np.array_equal(forward, reverse))
        self.assertEqual(feasible[6], [1, 3])
        self.assertEqual(feasible[5], [0, 1, 2])
        self.assertEqual(forward[0, 6], 3.)

    def test_five_atom_domination_all_feasible_entries_and_small_r_endpoints(self):
        regions = toy_bank()
        result = target.fit_profile_measure(regions, 18, (3, 8))
        bound = target.moment_curves(result['state_measure'], 18)
        checked = 0
        for m, values in enumerate(result['feasible_full_counts_by_occupancy']):
            for full in values:
                self.assertTrue(np.all(bound[:, m] >= regions[m, full].ravel()))
                checked += 4
        self.assertEqual(result['feasible_matrix_entries_checked'], checked)
        self.assertEqual(bound[0, 0], 0.)
        atoms, _ = target.common_moment_atoms(result['state_measure'])
        self.assertTrue(np.all(atoms <= 0))
        small = target.fit_profile_measure(toy_bank(3), 3, (0, 1))
        self.assertEqual(small['state_measure']['zero_secant_endpoints'], [2, 3])
        self.assertEqual(small['state_measure']['cross_secant_endpoints'], [2, 3])

    def test_above_one_atom_or_missed_zero_degree_rejected(self):
        regions = toy_bank(6)
        for (m, _), matrix in regions.items():
            if m: matrix[0, 0] = -50.+m
        with self.assertRaises(target.MomentFitError): target.fit_profile_measure(regions, 6, (1, 2))
        regions = toy_bank(6)
        regions[0, 0][0, 0] = .001
        with self.assertRaises(target.MomentFitError): target.fit_profile_measure(regions, 6, (1, 2))

    def test_mock_provider_one_call_and_no_profile_fit(self):
        with patch.object(target, 'dependencies', return_value=PINS), \
             patch.object(target.provider, 'forced_nonzero_regions', return_value=(toy_bank(), metadata())) as provider, \
             patch.object(target, 'fit_profile_measure', side_effect=AssertionError('No fit in bank command')):
            result = target.build_bank(max_seconds=1)
        self.assertEqual(provider.call_count, 1)
        self.assertEqual(provider.call_args.args, (F(9927, 10000), 18))
        self.assertTrue(result['complete'])
        with patch.object(target, 'dependencies', side_effect=[PINS, {}]), \
             patch.object(target.provider, 'forced_nonzero_regions', return_value=(toy_bank(), metadata())):
            with self.assertRaises(ValueError): target.build_bank()

    def test_profile_catalog_loader_counts_pins_and_independent_moment_check(self):
        with TemporaryDirectory() as directory, patch.object(target, 'dependencies', return_value=PINS):
            bank_path = Path(directory)/'bank.json'
            bank_path.write_text(json.dumps(target.make_bank(toy_bank(), metadata(), PINS)), encoding='utf-8')
            bank_sha = sha256(bank_path.read_bytes()).hexdigest()
            data = target.profile_inventory(bank_path, bank_sha)
            self.assertTrue(data['all_profiles_have_valid_measures'])
            self.assertEqual(data['profile_entries'], 99)
            self.assertEqual(data['unordered_measure_keys'], 54)
            path = Path(directory)/'profiles.json'
            def store(value):
                path.write_text(json.dumps(value, allow_nan=False), encoding='utf-8')
                return sha256(path.read_bytes()).hexdigest()
            digest = store(data)
            _, profiles = target.load_profile_measures(path, digest)
            self.assertEqual(len(profiles), 99)
            self.assertIs(profiles[1, 3]['state_measure'], profiles[3, 1]['state_measure'])
            self.assertTrue(profiles[1, 3]['source_count_is_upper'])
            count = target.profile_counts(18)[1, 3]
            self.assertEqual(profiles[1, 3]['source_count'], count)
            self.assertAlmostEqual(profiles[1, 3]['common_prefix_log2'],
                log2(count)-target.LIMIT*log2(target.MARK)-log2(255))
            with self.assertRaises(ValueError): target.load_profile_measures(path, '0'*64)
            bad = deepcopy(data)
            bad['measures'][0]['state_measure']['log2_scale11'] = -1000.
            with self.assertRaises(ValueError): target.load_profile_measures(path, store(bad))
            bad = deepcopy(data)
            bad['profiles'][0]['indexed_source_count_upper'] = '1'
            with self.assertRaises(ValueError): target.load_profile_measures(path, store(bad))
            bad = deepcopy(data)
            bad['reviewed_source_sha256'] = {}
            with self.assertRaises(ValueError): target.load_profile_measures(path, store(bad))

    def test_timeout_never_loads_as_complete_bank(self):
        with patch.object(target, 'dependencies', return_value=PINS), \
             patch.object(target.provider, 'forced_nonzero_regions', side_effect=TimeoutError('mock stop')):
            result = target.build_bank()
        self.assertFalse(result['complete'])
        self.assertEqual(result['rows'], [])
        with self.assertRaises(ValueError): target.validate_bank(result, PINS)


if __name__ == '__main__':
    unittest.main()
