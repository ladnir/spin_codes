"""Independent aggregate tests; fixtures are not distance certificates."""
from copy import deepcopy
from fractions import Fraction as F
import unittest

from gf8_rate_climb_certificate import combine
from gf8_rate_climb_profile import metadata as fixed_metadata
from gf8_periodic_climb_profile import metadata as periodic_metadata


def fixtures(dimension=80,pattern=(6,)):
    g=(fixed_metadata(dimension,pattern[0]) if len(pattern)==1
       else periodic_metadata(dimension,pattern))
    n=g['input_n']
    first=(3*n+99)//100
    tail=(9*n+9)//10
    last=(first-1)//(129-dimension)
    base=dict(**g,complete=True,rigorous=True,whole_code=False,margin_bits_lower=80)
    small=['0x1',hex(1<<81)]
    sparse=dict(**base,active_range=[1,last],failure_upper_dyadic_hex=small)
    middle=dict(**base,weights=[first,tail-1],failure_upper=str(F(1,1<<81)))
    lost=n-g['output_n']
    dense=dict(**base,weights=[tail,n],failure_upper_dyadic_hex=small,
               omitted_symbols=lost,completion_factor_log2=8*lost,
               dense_inner_hypothesis='invertible full causal inner; independent nonzero input diagonals')
    return sparse,middle,dense


class RateClimbAggregateAuditTests(unittest.TestCase):
    def test_fixed_and_periodic_exact_union(self):
        for dimension,pattern,rate in ((80,(6,),'5/6'),(67,(5,),'67/80'),
                                      (84,(6,6,6,7),'21/25'),(74,(5,6),'37/44')):
            result=combine(*fixtures(dimension,pattern))
            self.assertTrue(result['whole_code'])
            self.assertEqual(result['nominal_rate'],rate)
            self.assertEqual(result['margin_bits_lower'],79)
            self.assertLess(10*(result['minimum_distance']-1),result['output_n'])
            self.assertGreaterEqual(10*result['minimum_distance'],result['output_n'])

    def test_reject_sparse_and_dense_gaps(self):
        for index,key in ((0,'active_range'),(1,'weights'),(2,'weights')):
            receipts=list(deepcopy(fixtures()))
            endpoints=receipts[index][key]
            if index==0:
                endpoints[1]-=1
            elif index==1:
                endpoints[0]=(129-receipts[0]['outer'][1])*(receipts[0]['active_range'][1]+1)+1
            else:
                endpoints[0]+=1
            with self.assertRaises(ValueError):
                combine(*receipts)

    def test_reject_missing_completion_factor(self):
        for dimension,pattern in ((80,(6,)),(74,(5,6))):
            for key,value in (('completion_factor_log2',0),('omitted_symbols',0),
                              ('dense_inner_hypothesis','unconditioned heuristic')):
                receipts=list(deepcopy(fixtures(dimension,pattern)))
                receipts[2][key]=value
                with self.assertRaises(ValueError):
                    combine(*receipts)

    def test_reject_wrong_phase_metadata(self):
        for dimension,pattern,changed in ((84,(6,6,6,7),[6,7,6,6]),
                                         (74,(5,6),[6,5])):
            receipts=list(deepcopy(fixtures(dimension,pattern)))
            receipts[1]['output_pattern']=changed
            with self.assertRaises(ValueError):
                combine(*receipts)

    def test_period_two74_canonical_geometry_and_joins(self):
        sparse,middle,dense=fixtures(74,(5,6))
        expected=periodic_metadata(74,(5,6))
        for receipt in (sparse,middle,dense):
            for key,value in expected.items():
                self.assertEqual(receipt[key],value)
        self.assertEqual(sparse['active_range'],[1,991])
        self.assertEqual(middle['weights'],[54559,1636761])
        self.assertEqual(dense['weights'],[1636762,1818624])
        self.assertEqual(dense['omitted_symbols'],568320)
        self.assertEqual(dense['completion_factor_log2'],4546560)
        result=combine(sparse,middle,dense)
        self.assertEqual(result['padded_k'],1051392)
        self.assertEqual(result['requested_k'],1048576)
        self.assertEqual(result['rows'],14208)
        self.assertEqual(result['output_n'],1250304)
        self.assertEqual(result['minimum_distance'],125031)
        self.assertEqual(result['nominal_rate'],'37/44')
        self.assertEqual(result['requested_rate'],'1024/1221')
        self.assertEqual(result['first_remaining_outer_weight'],54560)

    def test_period_two74_reject_actual_boundary_gaps(self):
        for index,key,value in ((0,'active_range',[1,990]),
                                (1,'weights',[54561,1636761]),
                                (2,'weights',[1636763,1818624]),
                                (2,'completion_factor_log2',4546552)):
            receipts=list(deepcopy(fixtures(74,(5,6))))
            receipts[index][key]=value
            with self.assertRaises(ValueError):
                combine(*receipts)

    def test_reject_unsupported_margin_and_weak_union(self):
        receipts=list(deepcopy(fixtures()))
        receipts[1]['failure_upper']='1'
        with self.assertRaises(ValueError):
            combine(*receipts)
        receipts=list(deepcopy(fixtures()))
        for index,receipt in enumerate(receipts):
            receipt['margin_bits_lower']=40
            if index==1:
                receipt['failure_upper']=str(F(1,1<<40))
            else:
                receipt['failure_upper_dyadic_hex']=['0x1',hex(1<<40)]
        with self.assertRaises(ArithmeticError):
            combine(*receipts)


if __name__=='__main__':
    unittest.main()
