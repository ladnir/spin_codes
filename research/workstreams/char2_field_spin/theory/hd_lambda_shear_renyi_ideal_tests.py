from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb
from unittest.mock import patch
import unittest

from flint import arb, ctx

import hd_lambda_shear_renyi_ideal as ideal

MUL = ((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2))


class ShearRenyiIdealTests(unittest.TestCase):
    def test_full_space_and_zero_atom(self):
        for q,p in product((4,256),(2,3,4)):
            z,n = ideal.ideal_atoms(32,32,p,q=q)
            self.assertEqual((z,n),(F(q)**-32,1-F(q)**-32))
            z,_ = ideal.ideal_atoms(18,32,p,q=q)
            self.assertEqual(z,F(q)**(14*p-32))

    def test_jensen_atom_against_exact_toy_code(self):
        # Four output forms on two GF4 variables, with independent D2 after.
        forms = ((1,2),(3,1),(1,3),(1,0))
        counts = Counter()
        for x,y in product(range(4),repeat=2):
            mask = sum((MUL[a][x]^MUL[b][y] != 0) << i for i,(a,b) in enumerate(forms))
            counts[mask] += 1
        self.assertEqual(counts[0],1)
        self.assertEqual(sum(counts.values()),16)
        for p in (2,3,4):
            z,n = ideal.ideal_atoms(2,4,p,q=4)
            actual = F(4)**(p*(4-2)-4)*sum(
                (F(count**p,3**((p-1)*mask.bit_count())) for mask,count in counts.items() if mask),F(0))
            self.assertLessEqual(n,actual)
            self.assertEqual(z,F(4)**(2*p-4))
            self.assertLessEqual(ideal.positive_tail(z,n,blocks=4,minimum=2),
                                 ideal.positive_tail(z,actual,blocks=4,minimum=2))

    def test_positive_tail_exact_enumeration(self):
        z,n = F(3,7),F(5,11)
        expected = sum((z**(4-sum(bits))*n**sum(bits)
                        for bits in product((0,1),repeat=4) if sum(bits)>=2),F(0))
        self.assertEqual(ideal.positive_tail(z,n,blocks=4,minimum=2),expected)
        self.assertEqual(ideal.positive_tail(z,n,blocks=4,minimum=4),n**4)

    def test_pair_floor_is_five_not_nine(self):
        # Exhaust the tiny binary activity geometry, not field amplitudes.
        minimum = 8
        for mask in range(1<<16):
            if mask.bit_count() >= 9:
                active_pairs = sum(bool(mask & (3<<(2*j))) for j in range(8))
                minimum = min(minimum,active_pairs)
        self.assertEqual(minimum,5)

    def test_raw_beta_is_not_doubled(self):
        with ctx.workprec(256):
            for a,p in ((9,2),(12,3)):
                raw,data = ideal.raw_group_comparator(a,p)
                z,n = ideal.ideal_atoms(2*a,32,p)
                tail = ideal.positive_tail(z,n)
                beta = F(256)**(16*a)/F(255)**(8*a)
                exact_power = F(comb(16,a))**p*beta**p*tail
                self.assertEqual(F(data['beta']),beta)
                self.assertTrue((raw**p).contains(ideal.A.uniform.rational(exact_power)))

    def test_expression_prefix_and_projective_once(self):
        with ctx.workprec(256):
            threshold = -ideal.A.uniform.rational(ideal.ALPHA)*734003*ideal.A.log2(
                ideal.A.uniform.rational(ideal.MARK))
            for p in (2,3,4):
                theta = ideal.A.uniform.rational(ideal.ALPHA*F(p,p-1))
                actual = ideal.method_expression(arb(2)**7,p,arb(-23))
                expected = theta*(18*7-ideal.A.log2(arb(255)))+threshold-23
                self.assertTrue(actual.overlaps(expected))

    def test_mock_controls_select_lower_method_not_failure_bound(self):
        def raw(a,p):
            return arb(1),dict(arity=a,moment_order=p)
        def expression(weight,p,reference):
            return arb({2:20,3:-50,4:10}[p])
        with ctx.workprec(256), \
             patch.object(ideal,'raw_group_comparator',side_effect=raw), \
             patch.object(ideal,'method_expression',side_effect=expression):
            data = ideal.controls_from_reference(arb(0))
        self.assertEqual(len(data['controls']),12)
        self.assertEqual(len(data['selected_homogeneous_arities']),4)
        self.assertTrue(all(row['selected_lower_order']==3 for row in data['selected_homogeneous_arities']))
        self.assertTrue(all(not row['this_fixed_method_cannot_reach_40_bits']
                            for row in data['selected_homogeneous_arities']))
        self.assertTrue(all('failure' not in key for row in data['controls'] for key in row))

    def test_mock_receipt_loads_cache_without_provider(self):
        fake = dict(controls=[],selected_homogeneous_arities=[])
        with patch.object(ideal,'controls_from_reference',return_value=fake), \
             patch.object(ideal.A.uniform,'regional_directed',side_effect=AssertionError('No provider call')):
            data = ideal.receipt(256)
        self.assertEqual(data['source_profiles'],99)
        self.assertEqual(data['minimum_active_pairs'],5)
        self.assertEqual(data['pair_count'],8)
        self.assertEqual(data['fresh_regional_banks'],0)
        self.assertFalse(data['ideal_is_failure_bound'])
        self.assertFalse(data['actual_pair_moment_upper_bound'])
        self.assertFalse(data['whole_code_certificate'])
        self.assertFalse(data['mixed_arities_covered'])

    def test_guards_and_pins(self):
        self.assertEqual(ideal.dependencies()['hd_lambda_renyi_limits.py'],ideal.PINS['hd_lambda_renyi_limits.py'])
        for args in ((0,32,3),(33,32,3),(2,33,3),(2,4,5)):
            with self.assertRaises(ValueError):
                ideal.ideal_atoms(*args)
        with self.assertRaises(ValueError):
            ideal.positive_tail(F(1),F(1),blocks=8,minimum=9)
        with self.assertRaises(ValueError):
            ideal.raw_group_comparator(8,3)
        with self.assertRaises(ValueError):
            ideal.receipt(128)


if __name__ == '__main__':
    unittest.main()
