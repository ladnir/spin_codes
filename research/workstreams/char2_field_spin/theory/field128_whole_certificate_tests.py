import unittest

from field128_whole_certificate import combine


def fixtures():
    common=dict(outer=[128,87],q=256,width=16,state=8,k_log=20,rows=12096,
                rigorous=True,k=1052352,n=1548288,output_threshold=154828)
    sparse=dict(common,complete_sparse_range=True,active_range=[1,1105],interval_count=129,
                failure_upper_dyadic=['1',str(2**68)],partition_sha256='fixture')
    middle=dict(common,minimum_input_weight=46449,maximum_input_weight=371589,
                covered_weights=325141,failure_upper='1e-123',checked_cells=210,
                partition_sha256='fixture')
    dense=dict(common,dense_rigorous=True,cutoff_cover_complete=True,padded_source=1052352,
               output_length=1548288,forbidden_output_threshold=154828,
               dense_weights=[371590,1548288],dense_margin_bits_lower=13690,
               failure_upper_power_of_two=-13690)
    return sparse,middle,dense


class JoinTests(unittest.TestCase):
    def test_sum_not_minimum(self):
        r=combine(*fixtures())
        self.assertEqual(r['whole_code_margin_bits_lower'],67)
        self.assertEqual(r['minimum_distance'],154829)
        self.assertEqual(r['effective_requested_rate'],'128/189')
        self.assertEqual(r['binary_field_degrees'],[8,32,64,128])
        self.assertEqual(r['extension_degrees_over_gf256'],[1,4,8,16])

    def test_geometry_mismatch(self):
        for i in range(3):
            r=list(fixtures());r[i]['rows']=6016
            with self.assertRaises(ValueError):combine(*r)

    def test_missing_classes(self):
        for i,key,value in ((0,'active_range',[1,1104]),(1,'minimum_input_weight',46453),
                            (1,'maximum_input_weight',371588),(2,'dense_weights',[371591,1548288])):
            r=list(fixtures());r[i][key]=value
            with self.assertRaises(ValueError):combine(*r)

    def test_nonrigorous(self):
        for i,key in ((0,'rigorous'),(1,'rigorous'),(2,'dense_rigorous'),(2,'cutoff_cover_complete')):
            r=list(fixtures());r[i][key]=False
            with self.assertRaises(ValueError):combine(*r)

    def test_bad_probability(self):
        r=list(fixtures());r[0]['failure_upper_dyadic']=['1','1']
        with self.assertRaises(ArithmeticError):combine(*r)
        r=list(fixtures());r[1]['failure_upper']='0.01'
        with self.assertRaises(ArithmeticError):combine(*r)


if __name__=='__main__':unittest.main()
