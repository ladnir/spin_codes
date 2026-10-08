"""Exact aggregation checks; the full numerical replay is a separate CLI."""

from copy import deepcopy
from fractions import Fraction as F
import unittest

from field256_whole_certificate import combine, margin_bits


def receipts():
    # Synthetic upper bounds isolate aggregation from the proof engines.
    common = dict(outer=[256,176],q=256,width=16,state=8,k_log=20,rows=6016)
    geometry = dict(k=1058816,n=1540096,output_threshold=154009,rigorous=True)
    sparse = dict(**common,**geometry,complete_sparse_range=True,active_range=[1,600],
                  interval_count=80,failure_upper_dyadic=['1',str(2**120)])
    middle = dict(**common,**geometry,minimum_input_weight=46203,maximum_input_weight=369623,
                  covered_weights=323421,failure_upper='1e-104',checked_cells=216,
                  partition_sha256='synthetic')
    dense = dict(**common,dense_rigorous=True,cutoff_cover_complete=True,
                 padded_source=1058816,output_length=1540096,forbidden_output_threshold=154009,
                 dense_weights=[369624,1540096],dense_margin_bits_lower=12072)
    return sparse,middle,dense


class WholeCertificateTests(unittest.TestCase):
    def test_exact_integer_margin(self):
        for bits in (0,40,119,344,12072):
            if bits:
                self.assertEqual(margin_bits(F(1,2**bits)),bits)
            self.assertEqual(margin_bits(F(3,2**(bits+2))),bits)

    def test_union_is_summed_not_minimum_of_bits(self):
        result = combine(*receipts())
        self.assertEqual(result['component_margin_bits_lower'],[120,345,12072])
        self.assertEqual(result['whole_code_margin_bits_lower'],119)
        exact = F(1,2**120)+F('1e-104')+F(1,2**12072)
        self.assertGreaterEqual(F(result['failure_upper']),exact)
        self.assertTrue(result['whole_code'])
        self.assertFalse(result['native_performance_measured'])

    def test_geometry_and_missing_weights_rejected(self):
        for index,key,value in ((0,'state',4),(1,'width',32),(2,'output_length',1540095),
                                (0,'active_range',[1,599]),(1,'minimum_input_weight',46204),
                                (1,'maximum_input_weight',369622),(2,'dense_weights',[369625,1540096])):
            values = deepcopy(receipts())
            values[index][key] = value
            with self.subTest(index=index,key=key):
                with self.assertRaises(ValueError):
                    combine(*values)

    def test_incomplete_or_nondirected_replay_rejected(self):
        for index,key in ((0,'complete_sparse_range'),(0,'rigorous'),(1,'rigorous'),
                          (2,'dense_rigorous'),(2,'cutoff_cover_complete')):
            values = deepcopy(receipts())
            values[index][key] = False
            with self.assertRaises(ValueError):
                combine(*values)

    def test_invalid_or_insufficient_bounds_rejected(self):
        for numerator,denominator in (('0','1'),('-1','2'),('1','2'),('1',str(2**40))):
            values = receipts()
            values[0]['failure_upper_dyadic'] = [numerator,denominator]
            with self.assertRaises(ArithmeticError):
                combine(*values)


if __name__ == '__main__':
    unittest.main()
