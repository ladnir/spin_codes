"""Tiny algebra and inventory tests; no real regional provider run."""
from fractions import Fraction as F
from itertools import combinations
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_lambda_renyi_bridge as target


def inventory(moment=F(1)):
    return dict(target.MOMENT_SCOPE, complete=True,
        reviewed_source_sha256=target.MOMENT_PINS,
        source_sha256=target.MOMENT_PINS['hd_lambda_renyi_support.py'], source_pins_unchanged=True,
        rows=[dict(complement=list(c), moment=str(moment))
              for c in combinations(range(16), 4)])


class RenyiBridgeTests(unittest.TestCase):
    def test_complete_inventory_and_scope(self):
        data = inventory(F(5, 4))
        self.assertEqual(target.checked_moment_inventory(data), F(5, 4))
        data['rows'][12]['moment'] = '3/2'
        self.assertEqual(target.checked_moment_inventory(data), F(3, 2))
        for key, value in [('complete', False), ('field_order', 4), ('arity', 11),
                           ('moment_order', 2), ('exact_rational_arithmetic', False)]:
            wrong = dict(data, **{key: value})
            with self.assertRaises(ValueError):
                target.checked_moment_inventory(wrong)

    def test_missing_duplicate_and_bad_moment_rejected(self):
        data = inventory()
        for rows in (data['rows'][:-1], data['rows']+[data['rows'][0]]):
            with self.assertRaises(ValueError):
                target.checked_moment_inventory(dict(data, rows=rows))
        for value in ('99/100', 1):
            data = inventory()
            data['rows'][0]['moment'] = value
            with self.assertRaises(ValueError):
                target.checked_moment_inventory(data)

    def test_prefix_projective_and_density_norm(self):
        a = target.A
        with ctx.workprec(256):
            m = F(5, 4)
            expected = a.uniform.rational(F(3, 10))*(
                18*(1536-96*a.log2(arb(255))+a.log2(arb(1820))
                    + a.uniform.rational(F(16, 3))*a.log2(a.uniform.rational(m)))
                - a.log2(arb(255)))
            expected -= a.uniform.rational(F(1, 5))*734003*a.log2(a.uniform.rational(F(691, 700)))
            self.assertTrue(target.prefix(m).overlaps(expected))
            penalty = target.prefix(m)-target.prefix(1)
            self.assertTrue(penalty.overlaps(a.uniform.rational(F(144, 5))*a.log2(a.uniform.rational(m))))
        with self.assertRaises(ValueError):
            target.prefix(F(1, 2))

    def test_one_fresh_bank_and_source_count_once(self):
        a = target.A
        with ctx.workprec(256):
            terms = [(e, o, count, arb(count)*arb(2)**-20000)
                     for (e, o), count in a.profile_counts(18).items()]
            expected = target.prefix(1)+a.log2(arb(sum(t[2] for t in terms)))-20000
            expected = str(a.uniform.rational_ceiling(expected))
        with patch.object(a.uniform, 'regional_directed', return_value={}) as bank, \
             patch.object(a, 'source_profile_terms', return_value=terms):
            result = target.receipt(inventory(), 256)
        bank.assert_called_once_with(F(691, 700), F(1, 5), 18, 60)
        self.assertEqual(result['failure_log2_upper'], expected)
        self.assertEqual(result['profile_count'], 99)
        self.assertFalse(result['mixed_arities_covered'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertTrue(result['rowset_union_inside_theta'])

    def test_changed_dependencies_and_invalid_precision_rejected(self):
        with self.assertRaises(ValueError):
            target.receipt(inventory(), 128)
        a = target.A
        terms = [(e, o, count, arb(count)) for (e, o), count in a.profile_counts(18).items()]
        with patch.object(target, 'dependencies', side_effect=[{'pin': 'old'}, {'pin': 'new'}]), \
             patch.object(a.uniform, 'regional_directed', return_value={}), \
             patch.object(a, 'source_profile_terms', return_value=terms):
            with self.assertRaises(ValueError):
                target.receipt(inventory(), 256)


if __name__ == '__main__':
    unittest.main()
