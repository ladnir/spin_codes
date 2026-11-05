"""Exact root, mixture, floors, scope and replay wiring; no regional replay."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import prod
from pathlib import Path
import unittest
from unittest.mock import patch

from flint import arb, ctx

from hd_shifted_coverage import rational_ceiling
from hd_xor_dense_holder_classes import class_moments
from hd_xor_dense_holder_certificate import (EXPECTED_SCOPE, HELPER_PINS,
    PROPOSAL_SHA256, BASE_PROPOSAL_SHA256, SOURCE_COUNT, checked_partition,
    root_upper, dyadic_power_dominates, holder_branch_caps_arb, receipt,
    checked_rows, rounded_union_512, verify)
from hd_xor_mixed_mask_certificate import bivariate_tail_arb, _rational


def exact_dyadic(value):
    mantissa, exponent = value.man_exp()
    exponent = int(exponent)
    return F(int(mantissa) << exponent) if exponent >= 0 else F(int(mantissa), 1 << -exponent)


def fake_receipt():
    entries = []
    for identifier, label, z, logp, arity, masks in checked_partition():
        nums, denominator = class_moments(masks)
        digest = sha256(json.dumps([list(nums), denominator], separators=(',', ':')).encode()).hexdigest()
        entries.append(dict(type_id=identifier, witness_label=label, mark=str(z), log2_potential=str(logp),
            active_constituents=arity, active_row_masks=list(masks), fixed_row_sets=len(masks),
            exact_moments_sha256=digest, positive_mask_weights=sum(n > 0 for n in nums),
            root_upper_checks_passed=True, failure_log2_upper='-100'))
    return dict(EXPECTED_SCOPE, precision_bits=256, frozen_proposal_sha256=PROPOSAL_SHA256,
        base_proposal_sha256=BASE_PROPOSAL_SHA256, reviewed_helper_sha256=HELPER_PINS,
        source_sha256=sha256(Path(__file__).with_name('hd_xor_dense_holder_certificate.py').read_bytes()).hexdigest(),
        entries=entries)


class HolderCertificateTests(unittest.TestCase):
    def test_exact_root_endpoint_direction_and_small_excess(self):
        with ctx.workprec(128):
            for value in (F(0), F(1), F(2), F(1, 3), F(255)**-128,
                          F((1 << 801)+1, (1 << 279)-3)):
                upper = root_upper(value.numerator, value.denominator)
                dyadic = exact_dyadic(upper)
                self.assertGreaterEqual(dyadic**32, value)
                self.assertLessEqual(dyadic**32, value*(1+F(1, 1 << 100)))
                self.assertTrue(dyadic_power_dominates(upper, value.numerator, value.denominator))
            self.assertFalse(dyadic_power_dominates(arb(1), 2, 1))
            with self.assertRaises(ValueError):
                root_upper(-1, 1)
            with self.assertRaises(ValueError):
                root_upper(1, 0)

    def test_exact_tiny_class_moments_and_singleton_scaling(self):
        with ctx.workprec(128):
            nums, denominator = class_moments((3, 7), 4, 4, 3, 2)
            for n in nums:
                value = root_upper(n, denominator, 3)
                self.assertGreaterEqual(exact_dyadic(value)**3, F(n, denominator))
            # A one-element class has no Holder loss: cube root of
            # N(T)^3/3^|A|, here |A|=3, is N(T)/3 exactly.
            nums, denominator = class_moments((7,), 4, 4, 3, 2)
            from hd_xor_dense_holder import uniform_mask_count
            for mask, n in enumerate(nums):
                upper = exact_dyadic(root_upper(n, denominator, 3))
                exact = F(uniform_mask_count(4, 4, 7, mask), 3)
                self.assertGreaterEqual(upper, exact)
                self.assertLessEqual(upper, exact*(1+F(1, 1 << 100)))

    def test_weighted_mask_mixture_precedes_common_quotient(self):
        with ctx.workprec(128):
            vectors = {(0, 0): {0: (arb(1), arb(2)), 1: (arb(10), arb(1)), 2: (arb(1), arb(10))}}
            weights = {0: arb(3), 1: arb(2), 2: arb(5)}
            empty, nonempty = holder_branch_caps_arb(vectors, weights, arb(1))
            self.assertEqual(empty[0, 0], 6)
            self.assertEqual(nonempty[0, 0], 52)
            self.assertLess(nonempty[0, 0], 70)  # A max before mixing would be looser.
            with self.assertRaises(ValueError):
                holder_branch_caps_arb(vectors, {**weights, 3: arb(1)}, arb(1))

    def test_both_positive_floors_match_independent_word_enumeration(self):
        empty, nonempty = [F(1), F(2), F(3)], [F(4), F(1), F(2)]
        exact = F(0)
        terms = [(x, y, coefficient) for y, poly in enumerate((empty, nonempty))
                 for x, coefficient in enumerate(poly)]
        for word in product(terms, repeat=4):
            if sum(term[0] for term in word) >= 5 and sum(term[1] for term in word) >= 2:
                exact += prod(term[2] for term in word)
        with ctx.workprec(128):
            actual, _ = bivariate_tail_arb(list(map(_rational, empty)), list(map(_rational, nonempty)),
                                          columns=4, x_floor=5, y_floor=2)
            self.assertTrue(actual.contains(_rational(exact)))

    def test_partition_receipt_witness_and_scope_fail_closed(self):
        valid = fake_receipt()
        self.assertEqual(len(checked_rows(valid, 256)), 12)
        self.assertEqual(sum(row['fixed_row_sets'] for row in valid['entries']), 219)
        variants = []
        for key, value in (('source_matching_counts', [0, 8]), ('sparse_pair_assignments', 28),
                           ('dense_nonempty_column_floor', 16), ('extra_dense_beta_factor', True),
                           ('directed', 1), ('whole_code_certificate', True)):
            changed = deepcopy(valid)
            changed[key] = value
            variants.append(changed)
        for key, value in (('mark', '99/100'), ('active_row_masks', [255]),
                           ('root_upper_checks_passed', False), ('exact_moments_sha256', 'wrong')):
            changed = deepcopy(valid)
            changed['entries'][0][key] = value
            variants.append(changed)
        changed = deepcopy(valid)
        changed['entries'][0] = changed['entries'][1]
        variants.append(changed)
        for changed in variants:
            with self.assertRaises(ValueError):
                checked_rows(changed, 256)
        with patch('hd_xor_dense_holder_certificate.arb_pair_regions') as regional:
            with self.assertRaises(ValueError):
                receipt(scope=dict(EXPECTED_SCOPE, cohort_counts=[0, 15, 2]))
            regional.assert_not_called()

    def test_five_banks_twelve_tails_and_once_only_prefix(self):
        namespace = 'hd_xor_dense_holder_certificate.'
        with ctx.workprec(256), patch(namespace+'arb_pair_regions', return_value=None) as regional, \
                patch(namespace+'ordered_tables', return_value=None), \
                patch(namespace+'mixed_vectors_arb', return_value={}) as vectors, \
                patch(namespace+'weights_arb', return_value=({0: arb(1)}, 'mock')), \
                patch(namespace+'holder_branch_caps_arb', return_value=({}, {})), \
                patch(namespace+'column_coefficients_arb', return_value=[arb(1)]), \
                patch(namespace+'bivariate_tail_arb', return_value=(arb(1), [])) as tail:
            result = receipt()
            self.assertEqual(regional.call_count, 5)
            self.assertEqual(vectors.call_count, 5)
            self.assertEqual(tail.call_count, 12)
            for call in tail.call_args_list:
                self.assertEqual(call.kwargs, dict(columns=32, x_floor=544, y_floor=17))
            for row in result['entries']:
                z, logp = F(row['mark']), F(row['log2_potential'])
                potential = arb(2)**_rational(logp)
                exact_prefix = SOURCE_COUNT*17*28**16*_rational(z)**(-734003)/(255*potential*255**512)
                self.assertEqual(F(row['failure_log2_upper']),
                                 rational_ceiling(exact_prefix.log()/arb(2).log()))

    def test_independent_rounded_union_no_extra_type_factor(self):
        value = rounded_union_512([F(0), F(1), F(3)])
        with ctx.workprec(1024):
            self.assertEqual(value, rational_ceiling(arb(11).log()/arb(2).log()))
            self.assertGreater(rational_ceiling(arb(12*11).log()/arb(2).log())-value, 3)

    def test_observed_dual_receipts_and_fresh_rounded_union(self):
        result = verify()
        self.assertEqual(result['raw_aggregate_log2_upper'], '-5617715193493/100000000000')
        self.assertEqual(result['rounded_cell_union_log2_upper'], '-56177151934929/1000000000000')
        self.assertTrue(result['all_per_type_ceilings_agree'])
        self.assertEqual([row['sha256'] for row in result['inputs']], [
            '33ad1a1b1550b640497a97b1d8e68101fe9bd618c8868bb2579e9bbb76f82ced',
            'b264538ea3e0ba2df0712b7d21e65dc4ed7af86bb960c0cb3648ce746ec4813c'])
        self.assertFalse(result['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
