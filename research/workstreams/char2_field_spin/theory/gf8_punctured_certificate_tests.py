from copy import deepcopy
from fractions import Fraction as F
import unittest

from gf8_punctured_profile import metadata
from gf8_punctured_certificate import combine,validate_maps,MATRIX_SHA256


def fixtures():
    base=dict(**metadata(),complete=True,rigorous=True,whole_code=False,margin_bits_lower=80)
    small=['0x1',hex(1<<81)]
    sparse=dict(**base,active_range=[1,713],failure_upper_dyadic_hex=small)
    middle=dict(**base,weights=[32113,722534],failure_upper=str(F(1,1<<81)))
    dense=dict(**base,weights=[722535,1605632],failure_upper_dyadic_hex=small,
               omitted_symbols=200704,completion_factor_log2=1605632,
               dense_inner_hypothesis='invertible full causal inner; independent nonzero input diagonals')
    return sparse,middle,dense


class PuncturedWholeTests(unittest.TestCase):
    def test_exact_packet_binding(self):
        result=validate_maps()
        self.assertEqual(result['matrix_sha256'],MATRIX_SHA256)
        self.assertTrue(result['full_feedthrough_invertible'])
        self.assertTrue(result['retained_output_pivot_invertible'])

    def test_join_and_probability_union(self):
        result=combine(*fixtures())
        self.assertTrue(result['whole_code'])
        self.assertEqual(result['minimum_distance'],140493)
        self.assertEqual(result['first_remaining_outer_weight'],32130)
        self.assertEqual(result['margin_bits_lower'],79)

    def test_reject_gaps_and_missing_preimage_factor(self):
        for index,key,value in ((0,'active_range',[1,712]),
                                (1,'weights',[32131,722534]),
                                (2,'weights',[722536,1605632]),
                                (2,'completion_factor_log2',0)):
            receipts=list(deepcopy(fixtures()))
            receipts[index][key]=value
            with self.assertRaises(ValueError):
                combine(*receipts)

    def test_reject_metadata_and_probability_mismatch(self):
        for index,key,value in ((0,'outputs',8),(1,'output_n',1605632),
                                (2,'matrix_sha256','changed'),(0,'complete',False),
                                (1,'failure_upper','1'),(2,'rigorous',False)):
            receipts=list(deepcopy(fixtures()))
            receipts[index][key]=value
            with self.assertRaises(ValueError):
                combine(*receipts)


if __name__=='__main__':
    unittest.main()
