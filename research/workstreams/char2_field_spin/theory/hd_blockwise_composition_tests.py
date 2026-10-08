"""Capped-support placement, route polynomial, and directed receipt tests."""

from fractions import Fraction as F
from itertools import combinations,permutations,product
from math import comb
import unittest

from flint import arb,ctx

import hd_blockwise_composition as proof
from hd_blockwise_maps import local_transfer


class BlockwiseCompositionTests(unittest.TestCase):
    def test_mds_exact_support_not_shell_count(self):
        # GF4 GRS[4,2]: one zero word, four supports of weight3 with3
        # words each, and one weight4 support with3 words.
        counts=[proof.prescribed_count(4,2,w,4) for w in range(5)]
        self.assertEqual(counts,[1,0,0,3,3])
        self.assertEqual(sum(comb(4,w)*n for w,n in enumerate(counts)),16)
        # Fractional power never absorbs the four distinct support events.
        self.assertEqual(comb(4,3)**2*counts[3],48)
        self.assertNotEqual(comb(4,3)**2*counts[3],comb(4,3)*counts[3])

    def test_marked_power_noncommuting(self):
        a=[[F(1,2),F(1,4)],[F(0),F(3,4)]]
        b=[[F(1,3),F(0)],[F(1,5),F(2,3)]]
        self.assertNotEqual(proof.matmul(a,b),proof.matmul(b,a))
        for n in (1,2,3,7):
            powers=[proof.identity_matrix(2)]
            for _ in range(n):
                powers.append(proof.matmul(powers[-1],a))
            insertion=proof.zero_matrix(2)
            for j in range(n):
                insertion=proof.add(insertion,proof.matmul(proof.matmul(powers[j],b),powers[n-1-j]))
            self.assertEqual(proof.marked_power(a,b,n),(powers[n],insertion))

    def test_entrywise_fractional_path_inequality(self):
        # These squared rational entries have exact rational square roots.
        root=[[F(1,2),F(1,3)],[F(1,4),F(1,5)]]
        matrix=[[v*v for v in row] for row in root]
        for n in range(1,6):
            normal=proof.identity_matrix(2)
            fractional=proof.identity_matrix(2)
            for _ in range(n):
                normal=proof.matmul(normal,matrix)
                fractional=proof.matmul(fractional,root)
            self.assertLessEqual(sum(normal[0]),sum(fractional[0])**2)
        for x in (F(1,16),F(1,4),F(1),F(4),F(16)):
            self.assertLessEqual(min(F(1),x)**2,x)

    def test_region_polynomial_equals_explicit_support_and_route_sum(self):
        old=ctx.prec
        try:
            ctx.prec=160
            q,t,s,b,packets=4,2,1,4,2
            theta,rho=F(2,3),F(1,2)
            quiet=proof.entry_power(local_transfer(0,theta,t,s,q),rho)
            marked=proof.entry_power(local_transfer(1,theta,t,s,q),rho)
            absent,insertion=proof.marked_power(quiet,marked,packets)
            present=[[v/packets for v in row] for row in insertion]
            actual=proof.region_coefficients(absent,present,b)
            for weight in range(b+1):
                total=arb(0)
                for support in combinations(range(b),weight):
                    for positions in product(range(packets),repeat=weight):
                        locations=dict(zip(support,positions))
                        row=[arb(1),arb(0),arb(0)]
                        for region in range(b):
                            for packet in range(packets):
                                row=proof.apply_row(row,marked if locations.get(region)==packet else quiet)
                        total+=sum(row,arb(0))/packets**weight
                self.assertTrue(actual[weight].overlaps(total))
        finally:
            ctx.prec=old

    def test_route_type_probability(self):
        values=(0,1,1,2)
        distinct=set(permutations(values))
        self.assertEqual(proof.type_probability(values),F(1,len(distinct)))
        # The type-only bound cannot replace actual MDS amplitude conditions.
        constant_columns=[(1,1)]*4
        self.assertEqual(proof.information_set_type_bound(constant_columns,2),1)
        distinct_columns=[(0,1,2,3)]*4
        self.assertEqual(proof.information_set_type_bound(distinct_columns,2),F(1,24**2))

    def test_directed_singleton_and_type_gap_receipts(self):
        result=proof.singleton_receipt()
        self.assertEqual(result['margin_bits_floor'],135)
        self.assertEqual(result['source_class'],'exactly one nonzero outer row')
        self.assertFalse(result['whole_code'])
        self.assertTrue(result['all_coordinate_supports'])
        self.assertTrue(result['support_pattern_multiplicities_unpowered'])
        self.assertEqual(result['prescribed_support_count_power'],'3/4')
        control=proof.singleton_receipt(state=8)
        self.assertEqual(control['margin_bits_floor'],135)
        larger=proof.singleton_receipt(state=8,length=256,dimension=128,rows=8192,theta=F(9987,10000))
        self.assertEqual(larger['margin_bits_floor'],222)
        narrower=proof.singleton_receipt(state=8,width=32)
        self.assertEqual(narrower['margin_bits_floor'],143)
        narrowest=proof.singleton_receipt(state=4,width=16,length=256,dimension=128,rows=8192,
                                          theta=F(9987,10000))
        self.assertEqual(narrowest['margin_bits_floor'],223)
        gap=proof.type_gap_receipt()
        self.assertTrue(arb(gap['unavoidable_log2_gap'])>2000000)
        self.assertFalse(gap['actual_failure_lower_bound'])

    def test_common_fractional_moment_still_needs_capping(self):
        result=proof.fractional_obstruction()
        self.assertTrue(result['rules_out_common_rho_bound'])
        self.assertTrue(arb(result['fractional_sum_log2_lower'])>13880)
        self.assertFalse(result['actual_failure_lower_bound'])


if __name__=='__main__':
    unittest.main()
