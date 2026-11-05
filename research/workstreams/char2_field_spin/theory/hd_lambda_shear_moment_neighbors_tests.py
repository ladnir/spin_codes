from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
from math import comb
from unittest.mock import patch
import unittest

from flint import arb,ctx

import hd_lambda_shear_moment_neighbors as gate
from hd_lambda_column_shear_ranks import paired_rows,rank_gf4

MUL=((0,0,0,0),(0,1,2,3),(0,2,3,1),(0,3,1,2))


def toy_counts(columns,checks):
    counts=Counter()
    for values in product(range(4),repeat=len(columns)):
        syndrome=[0]*checks
        for column,value in zip(columns,values):
            for row in range(checks): syndrome[row]^=MUL[(column>>(2*row))&3][value]
        if not any(syndrome):
            counts[sum((x!=0)<<i for i,x in enumerate(values))]+=1
    return counts


def grouped_tails():
    return {a:{p:Counter({F(2)**p:comb(16,a)}) for p in gate.ORDERS} for a in gate.ARITIES}


class ShearMomentNeighborsTests(unittest.TestCase):
    def test_actual_parity_small_constructor_checks(self):
        for a,column in ((11,2),(13,8),(16,0)):
            source=tuple(range(a))
            parity=gate.parity_columns(source,column)
            generator=sum(paired_rows(source,column,1),())
            b=32-2*a
            self.assertEqual(rank_gf4(parity,b),b)
            for input_index in range(2*a):
                for check in range(b):
                    value=0
                    for h,g in zip(parity,generator):
                        value^=MUL[(h>>(2*check))&3][(g>>(2*input_index))&3]
                    self.assertEqual(value,0)
        self.assertEqual(gate.parity_columns(tuple(range(12)),6),
                         gate.old.cover.parity_columns(tuple(range(12)),6))

    def test_incremental_basis_against_independent_rank(self):
        for sequence in ((1,2,4,8,16,32),(1,5,6,9,10,15),(0,0,1,1,4,5)):
            basis=(0,)*3
            accepted=[]
            for value in sequence:
                updated=gate._extend(basis,value,3)
                independent=rank_gf4(tuple(accepted+[value]),3)>len(accepted)
                self.assertEqual(updated is not None,independent)
                if updated is not None:
                    basis=updated
                    accepted.append(value)

    def test_p2_p3_p4_witnesses_and_toy_exact_moments(self):
        for columns,b in (((1,1,4,4),2),((1,4,5,6,9),2),((0,1,1,1),1)):
            counts=toy_counts(columns,b)
            n=len(columns)
            for p in gate.ORDERS:
                witness=gate.independent_partition(columns,b,p,max_attempts=2)
                checked=gate.verify_partition(columns,b,witness['colors'],p)
                self.assertEqual(checked['sizes'],witness['sizes'])
                zero,nonzero,upper=gate.moment_upper_from_sizes(n,b,witness['sizes'],q=4)
                actual=F(4)**(p*b-n)*sum((F(v**p,3**((p-1)*s.bit_count()))
                    for s,v in counts.items()),F(0))
                self.assertEqual(zero,F(4)**(p*b-n))
                self.assertEqual(nonzero,upper-zero)
                self.assertGreaterEqual(upper,actual)
                self.assertFalse(witness['maximum_coverage_claimed'])

    def test_general_support_sum_and_a12_overlap(self):
        colors=((0,1),(2,),(3,),())
        for p in gate.ORDERS:
            used=colors[:p]
            sizes=tuple(map(len,used))
            expected=F(0)
            for mask in range(1<<5):
                rank=max(sum((mask>>i)&1 for i in color) for color in used)
                expected+=F(3)**(mask.bit_count()-p*rank)
            expected*=F(4)**(2*p-5)
            self.assertEqual(gate.moment_upper_from_sizes(5,2,sizes,q=4)[2],expected)
        for sizes in ((8,8,8,8),(8,8,7,7),(5,4,3,2)):
            self.assertEqual(gate.moment_upper_from_sizes(32,8,sizes)[2],
                             gate.old.cover.moment_upper_from_sizes(32,8,sizes))

    def test_uniform_edge_and_zero_is_not_always_one(self):
        for p in gate.ORDERS:
            witness=gate.independent_partition((0,)*32,0,p,max_attempts=1)
            self.assertEqual(witness['sizes'],[0]*p)
            zero,nonzero,upper=gate.moment_upper_from_sizes(32,0,[0]*p)
            self.assertEqual((zero,nonzero,upper),(F(256)**-32,1-F(256)**-32,F(1)))
        self.assertEqual(gate.moment_upper_from_sizes(32,10,(8,)*4)[0],F(256)**8)
        self.assertEqual(gate.moment_upper_from_sizes(32,8,(8,)*2)[0],F(256)**-16)

    def test_positive_activity_tail_nonunit_zero(self):
        zero=F(3,7)
        atoms=tuple(F(i+1,13) for i in range(8))
        expected=F(0)
        for bits in product((0,1),repeat=8):
            if sum(bits)>=5:
                term=F(1)
                for bit,atom in zip(bits,atoms): term*=atom if bit else zero
                expected+=term
        self.assertEqual(gate.activity_tail(zero,atoms),expected)

    def test_uniform_record_validation_and_tamper(self):
        # Zero parity performs no rank search and tests all exact a16 atoms.
        row=gate.source_record(tuple(range(16)),max_attempts=1)
        tails=gate.validate_record(row,tuple(range(16)))
        self.assertEqual(set(tails),{2,3,4})
        for key in ('zero','nonzero','tail','source'):
            bad=deepcopy(row)
            if key=='zero': bad['pairs'][0]['orders']['2']['zero_moment']='1'
            if key=='nonzero': bad['pairs'][0]['orders']['4']['nonzero_moment_upper']='1'
            if key=='tail': bad['restricted_word_moment_uppers']['3']='1'
            if key=='source': bad['source'][0],bad['source'][1]=1,0
            with self.assertRaises(ValueError): gate.validate_record(bad,tuple(range(16)))

    def test_raw_weights_have_no_duplicate_binomial(self):
        with ctx.workprec(256):
            for p in gate.ORDERS:
                weights,diagnostics=gate.raw_weights(grouped_tails(),p)
                for a,weight,row in zip(gate.ARITIES,weights,diagnostics):
                    beta=F(256)**(16*a)/F(255)**(8*a)
                    self.assertTrue(weight.contains(gate.old.ideal.A.uniform.rational(2*comb(16,a)*beta)))
                    self.assertEqual(F(row['beta']),beta)
        bad=grouped_tails()
        bad[11][2][F(4)]-=1
        with self.assertRaises(ValueError): gate.raw_weights(bad,2)

    def test_mixed_coefficients_keep_p_fixed_and_range_index(self):
        # Positive mock weights1..6 make the total198 edge equal1 and
        # total288 equal6^18. The threshold/source formula is real Arb.
        def weights(grouped,p): return [arb(i) for i in range(1,7)],[]
        with ctx.workprec(256),patch.object(gate,'raw_weights',side_effect=weights):
            data=gate.controls_from_reference({},arb(-7000),total_min=287,total_max=288)
            threshold=-gate.old.ideal.A.uniform.rational(F(1,5))*734003*gate.old.ideal.A.log2(
                gate.old.ideal.A.uniform.rational(F(691,700)))
            for control in data['controls']:
                p=control['moment_order']
                expected=gate.old.ideal.limits.total_log_bound(arb(6)**18,p,arb(-7000),threshold)
                upper=F(control['totals'][1]['failure_log2_upper'])
                self.assertEqual(upper,gate.old.ideal.A.uniform.rational_ceiling(expected))
                self.assertEqual([r['total_arity'] for r in control['totals']],[287,288])
            for index,selected in enumerate(data['selected_totals']):
                self.assertEqual(F(selected['failure_log2_upper']),
                    min(F(c['totals'][index]['failure_log2_upper']) for c in data['controls']))

    def test_mock_inventory_scope_timeout_and_guards(self):
        with patch.object(gate,'source_record',side_effect=lambda source,**k:dict(source=list(source))):
            data=gate.inventory()
        self.assertTrue(data['complete'])
        self.assertEqual(len(data['rows']),6885)
        self.assertEqual(Counter(map(len,gate.SOURCES)),{a:comb(16,a) for a in gate.ARITIES})
        with patch.object(gate,'source_record',side_effect=TimeoutError('tiny mock timeout')):
            partial=gate.inventory()
        self.assertFalse(partial['complete'])
        with self.assertRaises(ValueError): gate.validate_inventory(partial)
        for args in ((197,288),(198,289),(250,249)):
            with self.assertRaises(ValueError): gate._total_range(*args)
        with self.assertRaises(ValueError): gate.moment_upper_from_sizes(32,10,(10,)*5)
        with self.assertRaises(ValueError): gate.parity_columns(tuple(range(10)),0)
        with self.assertRaises(ValueError): gate.verify_partition((1,1),1,((0,),(0,)),2)
        self.assertEqual(gate.dependencies()['hd_lambda_shear_moment_all_a12.py'],
                         gate.PINS['hd_lambda_shear_moment_all_a12.py'])


if __name__=='__main__': unittest.main()
