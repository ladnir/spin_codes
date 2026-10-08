import unittest
from fractions import Fraction as F

from hd_blockwise_shared_transfer import (multiply,power,rank,geometry,transition_bounds,
    successors_exact,all_subspaces,beta64_line_receipt)


class SharedTransferTests(unittest.TestCase):
    def test_exhaustive_joint_subspaces(self):
        checked=0
        for q,expected in ((2,67),(4,529)):
            count=0
            for joint in all_subspaces(4,q):
                profile=geometry(joint,q)
                exact=successors_exact(joint,q)
                bound=transition_bounds(profile,tau=F(3,2))
                r=profile['difference_rank']
                d=profile['joint_dimension']
                self.assertEqual(exact.get(r+1,F(0)),bound['up'])
                self.assertLessEqual(bound['down_lower'],exact.get(r-1,F(0)))
                self.assertLessEqual(exact.get(r-1,F(0)),bound['down_upper'])
                self.assertTrue(set(exact)<={r-1,r,r+1})
                for key,value in exact.items():
                    self.assertLessEqual(value,bound['positive_rank_weights'].get(key,F(0)))
                actual_return=sum((v for key,v in exact.items() if key<d),F(0))
                self.assertLessEqual(actual_return,bound['capped_live_return_upper'])
                if profile['diagonal_dimension']:
                    self.assertEqual(actual_return,bound['capped_live_return_upper'])
                moment=sum((v*F(3,2)**(-key) for key,v in exact.items()),F(0))
                self.assertLessEqual(moment,bound['rank_moment_upper'])
                count+=1
            self.assertEqual(count,expected)
            checked+=count
        self.assertEqual(checked,596)

    def test_two_dimensional_diagonal_always_returns(self):
        joint=((1,0),(0,1),(0,0),(0,0))*2
        result=transition_bounds(geometry(joint),tau=65536)
        self.assertEqual(result['profile']['diagonal_dimension'],2)
        self.assertEqual(result['capped_live_return_upper'],1)
        self.assertEqual(set(result['positive_rank_weights']),{0,1})

    def test_saved_beta64_line(self):
        r=beta64_line_receipt()
        self.assertEqual(r['profile']['diagonal_dimension'],0)
        self.assertEqual(r['profile']['difference_rank'],1)
        self.assertEqual(r['down_lower'],r['exact_return_probability'])
        self.assertEqual(r['down_upper'],F(1,65536**4-1))
        self.assertFalse(r['source_word_union'])

    def test_true_quadratic_field_arithmetic_and_geometry(self):
        self.assertEqual(multiply(256,256,65536),256^0x20)
        for x in (1,2,255,256,257,0x1234,65535):
            self.assertEqual(multiply(x,power(x,65534,65536),65536),1)
        # Genuine extension coefficients, not only embedded GF256 values.
        q=(256,1,0,0)
        c=(1,256,1,0)
        result=transition_bounds(geometry(tuple((x,) for x in q+c)))
        self.assertEqual(result['profile']['joint_dimension'],1)
        self.assertEqual(rank((q,c),65536),2)
        self.assertEqual(result['capped_live_return_upper'],F(1,65536**4-1))

    def test_duplicate_source_columns_do_not_multiply_event(self):
        joint=((1,2),(0,0),(0,0),(0,0),(0,0),(1,2),(0,0),(0,0))
        single=tuple((row[0],) for row in joint)
        self.assertEqual(geometry(joint),geometry(single))
        self.assertEqual(transition_bounds(geometry(joint)),transition_bounds(geometry(single)))


if __name__=='__main__':
    unittest.main()
