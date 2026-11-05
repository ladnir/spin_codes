"""Tiny directed algebra and mock gates; no production tail/provider runs."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_lambda_renyi_marked_control as target


def exact_convolution(left, right):
    result = [F(0)]*(len(left)+len(right)-1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i+j] += x*y
    return result


class MarkedControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = target.DIRECTORY/target.INVENTORY_NAME
        cls.raw_inventory = json.loads(cls.path.read_bytes())

    def test_complete_observed_inventory_and_all_witnesses(self):
        data = target.load_inventory(self.path)
        self.assertEqual(len(data['source_index']), 6884)
        self.assertEqual(len(data['polynomials']), 299)
        self.assertEqual(data['inventory_sha256'], target.INVENTORY_SHA256)
        self.assertEqual(target.dependencies()['hd_lambda_renyi_marked.py'], target.PINS['hd_lambda_renyi_marked.py'])

    def test_inventory_scope_coefficients_and_witness_corruption_rejected(self):
        variants = []
        for key, value in [('complete', False), ('arities', [True, 2, 3, 4, 5]),
                           ('input_floor_applied', True), ('reviewed_source_sha256', {}),
                           ('maps', [])]:
            row = deepcopy(self.raw_inventory)
            row[key] = value
            variants.append(row)
        row = deepcopy(self.raw_inventory)
        row['representatives'][0]['coefficients']['2'][0] = '1'
        variants.append(row)
        row = deepcopy(self.raw_inventory)
        row['representatives'][0]['coefficients']['2'][1] = '-1'
        variants.append(row)
        row = deepcopy(self.raw_inventory)
        row['representatives'][0]['coefficients']['2'].append('0')
        variants.append(row)
        row = deepcopy(self.raw_inventory)
        row['rows'][5]['representative_to_source_permutation'] = list(range(16))
        # Ensure this is an actual nonidentity witness before corrupting it.
        if row['rows'][5] == self.raw_inventory['rows'][5]:
            row['rows'][5]['representative_to_source_permutation'] = list(reversed(range(16)))
        variants.append(row)
        row = deepcopy(self.raw_inventory)
        row['rows'][-1] = row['rows'][0]
        variants.append(row)
        for data in variants:
            with self.assertRaises(ValueError):
                target.checked_inventory(data)
        with patch.object(Path, 'read_bytes', return_value=b'{}'):
            with self.assertRaisesRegex(ValueError, 'hash changed'):
                target.load_inventory(self.path)

    def test_positive_tail_matches_exact_labelled_expansion(self):
        polynomials = ((F(1, 3), F(2), F(0)), (F(7), F(2, 5)), (F(0), F(11, 13)))
        for precision in (256, 384):
            with ctx.workprec(precision):
                for floor in range(7):
                    exact = F(0)
                    for degrees in product(*(range(len(poly)) for poly in polynomials)):
                        if sum(degrees) >= floor:
                            term = F(1)
                            for poly, degree in zip(polynomials, degrees):
                                term *= poly[degree]
                            exact += term
                    upper = target.exact_upper(target.positive_coefficient_tail(polynomials, floor))
                    self.assertGreaterEqual(upper, exact)
                    self.assertLessEqual(upper-exact, max(F(1), exact)*F(2)**-200)

    def test_extreme_coefficient_ranges_keep_small_positive_tail(self):
        # A total-minus-prefix computation at256 bits loses this tail.
        huge, tiny = F(2)**3000, F(1, 255**21)
        for precision in (256, 384):
            with ctx.workprec(precision):
                upper = target.exact_upper(target.positive_coefficient_tail(((huge, tiny),)*4, 4))
                self.assertGreaterEqual(upper, tiny**4)
                self.assertLess(upper, tiny**4*(1+F(2)**-200))
                self.assertEqual(target.exact_upper(target.positive_coefficient_tail(((1, 0),)*2, 1)), 0)

    def test_scheduled_root_multiplicity_and_RS_factor_outside_root(self):
        # Every singleton's schedule visits all16 labels, hence one signature
        # with multiplicity16. This is fake low-degree arithmetic, not a gate.
        checked = dict(source_index={(i,): i for i in range(16)},
            polynomials=tuple({2: (F(1), F(0), F(i+1))} for i in range(16)))
        signatures = target.scheduled_signatures(checked, 1)
        self.assertEqual(signatures, {tuple(range(16)): 16})
        exact = [F(1)]
        for i in range(16):
            exact = exact_convolution(exact, (F(1), F(0), F(i+1)))
        tail = sum(exact[18:])
        for precision in (256, 384):
            with ctx.workprec(precision):
                weights, diagnostics = target.scheduled_weights(checked, 2, arities=(1,))
                upper = target.exact_upper(weights[1])
                self.assertGreaterEqual((upper*255**8/16)**2, tail)
                self.assertLess((upper*255**8/16)**2, tail*(1+F(2)**-190))
                self.assertEqual(diagnostics[0]['labelled_rowsets'], 16)
                self.assertEqual(diagnostics[0]['distinct_scheduled_signatures'], 1)
                self.assertEqual(diagnostics[0]['input_degree_floor'], 18)

    def test_prefix_uses_raw_weights_no_beta_and_single_projective_factor(self):
        with ctx.workprec(256):
            for p in target.ORDERS:
                theta = F(p, 5*(p-1))
                weight = arb(7)/3
                expected = target.A.uniform.rational(theta)*(18*target.A.log2(weight)-target.A.log2(arb(255)))
                expected -= target.A.uniform.rational(F(1, 5))*734003*target.A.log2(target.A.uniform.rational(F(691, 700)))
                self.assertTrue(target.class_prefix(weight, p).overlaps(expected))

    def test_mock_gate_one_bank_all_profiles_and_common_order(self):
        def fake_weights(checked, order, **kwargs):
            return {a: arb(a+order) for a in target.ARITIES}, []
        counts = target.A.profile_counts(18)
        with ctx.workprec(256):
            terms = [(e, o, n, arb(n)*arb(2)**-20000) for (e, o), n in counts.items()]
            uniform = target.A.log2(arb(sum(counts.values())))-20000
            expected = {p: str(target.A.uniform.rational_ceiling(
                target.class_prefix(arb(sum(a+p for a in target.ARITIES)), p)+uniform)) for p in target.ORDERS}
        with patch.object(target, 'scheduled_weights', side_effect=fake_weights), \
                patch.object(target.A.uniform, 'regional_directed', return_value={}) as bank, \
                patch.object(target.A, 'source_profile_terms', return_value=terms):
            result = target.receipt(self.path, 256)
        bank.assert_called_once_with(F(691, 700), F(1, 5), 18, 60)
        self.assertEqual(result['profile_count'], 99)
        self.assertEqual(result['fresh_regional_banks'], 1)
        for control in result['controls']:
            self.assertEqual(control['mixed']['failure_log2_upper'], expected[control['moment_order']])
            self.assertEqual(len(control['homogeneous']), 5)
        self.assertEqual(result['selected_mixed']['failure_log2_upper'], str(min(map(F, expected.values()))))
        self.assertTrue(result['no_extra_density_beta'])
        self.assertTrue(result['homogeneous_classes_not_added_to_mixed'])
        self.assertFalse(result['whole_code_certificate'])

    def test_weights_only_has_no_provider_and_pin_drift_is_rejected(self):
        def fake_weights(checked, order, **kwargs):
            return {a: arb(a+order) for a in target.ARITIES}, []
        with patch.object(target, 'scheduled_weights', side_effect=fake_weights), \
                patch.object(target.A.uniform, 'regional_directed') as bank:
            result = target.receipt(self.path, 384, weights_only=True)
        bank.assert_not_called()
        self.assertEqual(result['fresh_regional_banks'], 0)
        self.assertEqual(result['controls'], [])
        self.assertIsNone(result['uniform_profile_union_log2_upper'])
        with patch.object(target, 'dependencies', side_effect=({'a': 'old'}, {'a': 'new'})), \
                patch.object(target, 'scheduled_weights', side_effect=fake_weights):
            with self.assertRaisesRegex(ValueError, 'source pins changed'):
                target.receipt(self.path, 256, weights_only=True)

    def test_guards_and_incomplete_profiles(self):
        for polynomials, floor in [((), 0), (((-1, 2),), 0), (((1, 2),), -1)]:
            with self.assertRaises(ValueError):
                target.positive_coefficient_tail(polynomials, floor)
        for p in (1, 5, True):
            with self.assertRaises(ValueError):
                target.class_prefix(arb(1), p)
        for precision in (128, True):
            with self.assertRaises(ValueError):
                target.receipt(self.path, precision)
        for budget in (0, 1801, float('nan'), True):
            with self.assertRaises(ValueError):
                target.receipt(self.path, 256, max_seconds=budget)
        with self.assertRaises(ValueError):
            target._profile_rows([])


if __name__ == '__main__':
    unittest.main()
