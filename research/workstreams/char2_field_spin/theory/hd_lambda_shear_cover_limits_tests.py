from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb,prod
from unittest.mock import patch
import unittest

from flint import arb,ctx

import hd_lambda_shear_cover_limits as gate
import hd_lambda_shear_moment_neighbors as frozen


class CoverLimitsTests(unittest.TestCase):
    def test_max_level_sum_against_cartesian_gf4(self):
        for sizes in ((0,0),(1,2),(2,3),(0,1,3),(2,2,2),(0,0,1,2),(1,2,2,3)):
            p = len(sizes)
            expected = sum((F(prod(comb(s,k) for s,k in zip(sizes,ks)))
                *F(3)**(sum(ks)-p*max(ks))
                for ks in product(*(range(s+1) for s in sizes))),F(0))
            self.assertEqual(gate.size_sum(sizes,4),expected)

    def test_exact_support_sum_and_frozen_formula(self):
        for sizes in ((2,1),(2,1,1),(2,1,1,0)):
            p,offset,colors = len(sizes),0,[]
            for size in sizes:
                colors.append(set(range(offset,offset+size)))
                offset += size
            direct = sum((F(3)**(mask.bit_count()-p*max(
                sum((mask>>i)&1 for i in color) for color in colors))
                for mask in range(1<<5)),F(0))*F(4)**(2*p-5)
            self.assertEqual(gate.abstract_cover_expression(5,2,sizes,q=4),direct)
        for sizes in ((8,8),(8,8,8),(8,8,8,8),(10,9,7,6),(0,0,0,0)):
            self.assertEqual(gate.abstract_cover_expression(32,10,sizes),
                             frozen.moment_upper_from_sizes(32,10,sizes)[2])

    def test_profiles_are_complete_sorted_and_capacity_saturated(self):
        for width,checks,p in ((5,2,2),(5,2,3),(6,3,4),(4,0,3),(3,3,2)):
            expected = sorted({tuple(sorted(s)) for s in product(range(checks+1),repeat=p)
                               if sum(s)==min(width,p*checks)})
            actual = list(gate.size_profiles(width,checks,p))
            self.assertEqual(actual,expected)

    def test_adding_abstract_covered_coordinate_is_monotone(self):
        for sizes in ((0,0),(1,1),(0,1,1),(1,1,1,0)):
            old = gate.abstract_cover_expression(6,3,sizes,q=4)
            for i,s in enumerate(sizes):
                if s < 3 and sum(sizes)<6:
                    updated = list(sizes)
                    updated[i] += 1
                    self.assertLessEqual(gate.abstract_cover_expression(6,3,updated,q=4),old)

    def test_optimistic_minimum_is_no_larger_than_every_toy_profile(self):
        for p in (2,3,4):
            row = gate.optimistic_cover(6,3,p,q=4)
            actual = F(row['minimum_abstract_expression'])
            values = [gate.abstract_cover_expression(6,3,s,q=4)
                for s in product(range(4),repeat=p) if sum(s)<=6]
            self.assertEqual(actual,min(values))
            self.assertEqual(F(row['optimistic_nonzero_atom']),actual-F(row['zero_pair_moment']))
            self.assertFalse(row['actual_witness_claimed'])

    def test_raw_normalization_once_without_real_size_inventory(self):
        toy = dict(zero_pair_moment='1',optimistic_nonzero_atom='2')
        with ctx.workprec(256),patch.object(gate,'optimistic_cover',return_value=toy):
            raw,row = gate.raw_group_comparator(4,2)
            beta = F(256)**64/F(255)**32
            tail = gate.reference.positive_tail(F(1),F(2))
            expected = comb(16,4)*gate.reference.A.uniform.rational(beta)*gate.reference.A.uniform.rational(tail).sqrt()
            self.assertTrue((raw-expected).contains(0))
            self.assertEqual(row['labelled_rowsets'],comb(16,4))
            self.assertEqual(F(row['beta']),beta)

    def test_mock_controls_keep_full_order_choices(self):
        def raw(a,p,**kwargs):
            return arb(2)**(a+p),dict(arity=a,moment_order=p)
        with ctx.workprec(256),patch.object(gate,'raw_group_comparator',side_effect=raw):
            data = gate.controls_from_reference(arb(-7000))
        self.assertEqual(len(data['controls']),21)
        self.assertEqual(len(data['selected_homogeneous_arities']),7)
        for selected in data['selected_homogeneous_arities']:
            rows = [row for row in data['controls'] if row['arity']==selected['arity']]
            self.assertEqual(F(selected['method_expression_log2_lower']),
                             min(F(row['method_expression_log2_lower']) for row in rows))

    def test_mock_receipt_scope_and_no_provider(self):
        cache = b'tiny cached enclosure'
        bank = dict(file='unused.json',sha256=sha256(cache).hexdigest())
        with patch.object(gate,'dependencies',return_value={'test':'pin'}), \
                patch.object(gate.reference,'checked_reference',return_value=(arb(-7000),bank)), \
                patch.object(gate,'controls_from_reference',return_value={'controls':[],'selected_homogeneous_arities':[]}), \
                patch.object(gate.Path,'read_bytes',return_value=cache), \
                patch.object(gate.reference.A.uniform,'regional_directed',side_effect=AssertionError('No provider')):
            row = gate.receipt()
        self.assertTrue(row['complete'])
        self.assertTrue(row['comparator_is_lower_on_prescribed_method'])
        self.assertFalse(row['failure_probability_upper_bound'])
        self.assertFalse(row['comparator_is_failure_probability_lower_bound'])
        self.assertFalse(row['construction_rejected'])
        self.assertFalse(row['actual_pair_moment_bound'])

    def test_guards_and_timeout(self):
        for sizes in ((1,),(-1,2),(33,0)):
            with self.assertRaises(ValueError): gate.size_sum(sizes,4)
        with self.assertRaises(ValueError): gate.abstract_cover_expression(4,2,(2,2,1))
        with self.assertRaises(ValueError): gate.raw_group_comparator(11,2)
        with self.assertRaises(ValueError): list(gate.size_profiles(4,5,2))
        with self.assertRaises(TimeoutError):
            gate.optimistic_cover(6,3,2,check_time=lambda:(_ for _ in ()).throw(TimeoutError()))
        for budget in (0,float('inf'),float('nan')):
            with self.assertRaises(ValueError): gate.receipt(max_seconds=budget)


if __name__ == '__main__':
    unittest.main()
