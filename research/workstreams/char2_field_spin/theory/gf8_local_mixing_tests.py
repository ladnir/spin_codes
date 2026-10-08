"""Independent small-field enumeration for the local-mixer laws."""
from fractions import Fraction as F
from itertools import product
import unittest

from gf8_local_mixing_screen import (accounting, archived_flat_witness,
    fixed_mds_pair_moment, gl_pair_nonzero_moment, gl_pair_rank_moment,
    matching_moment, matching_zero_coefficients, matchings)


def mul(a,b):
    result=0
    while b:
        if b&1:result^=a
        a<<=1
        if a&4:a^=7
        b>>=1
    return result


def gl2():
    return tuple((a,b,c,d) for a,b,c,d in product(range(4),repeat=4)
                 if mul(a,d)^mul(b,c))


class LocalMixingTests(unittest.TestCase):
    def test_uniform_gl2_pair_law(self):
        matrices=gl2()
        self.assertEqual(len(matrices),(16-1)*(16-4))
        y=F(3,5)
        for u,v in product(range(4),repeat=2):
            actual=sum(y**(bool(mul(a,u)^mul(b,v))+bool(mul(c,u)^mul(d,v)))
                       for a,b,c,d in matrices)/len(matrices)
            self.assertEqual(actual,gl_pair_nonzero_moment(y,4) if u or v else 1)

    def test_fixed_mds_pair_law(self):
        y=F(3,5)
        for u,v in product(range(4),repeat=2):
            actual=F(0)
            for d0,d1 in product(range(1,4),repeat=2):
                a,b=mul(d0,u),mul(d1,v)
                actual+=y**(bool(a^b)+bool(a^mul(2,b)))
            actual/=9
            self.assertEqual(actual,fixed_mds_pair_moment(bool(u)+bool(v),y,4))
        self.assertGreater(fixed_mds_pair_moment(2,y),fixed_mds_pair_moment(1,y))

    def test_all_local_matchings(self):
        width=8
        pairings=tuple(matchings(range(width)))
        self.assertEqual(len(pairings),105)
        for q,y in ((4,F(3,5)),(256,F(9,10))):
            g=gl_pair_nonzero_moment(y,q)
            for mask in range(1<<width):
                moment=sum(g**sum(bool(mask&(1<<a) or mask&(1<<b)) for a,b in pairs)
                           for pairs in pairings)/len(pairings)
                self.assertEqual(moment,matching_moment(mask.bit_count(),y,q,width))

    def test_positive_zero_expansion(self):
        for width in (2,4,8,16):
            for y in (F(1,10),F(1,2),F(9,10),F(1)):
                coefficients=matching_zero_coefficients(y,256,width)
                self.assertTrue(all(c>=0 for c in coefficients))
                self.assertTrue(all(coefficients[z]==0 for z in range(1,width+1,2)))
                self.assertEqual(matching_moment(0,y,256,width),1)

    def test_rank_moment(self):
        # L=zero, span(e1), or identity. Enumerate the fresh selected row.
        matrices=gl2()
        tau=F(3)
        for form_rank in (0,1,2):
            for selected in (0,1,2):
                total=F(0)
                for a,b,c,d in matrices:
                    if not form_rank or not selected:r=0
                    elif selected==2:r=form_rank
                    elif form_rank==1:r=int(bool(a))
                    else:r=1
                    total+=tau**(-r)
                self.assertEqual(total/len(matrices),gl_pair_rank_moment(form_rank,selected,tau,4))

    def test_archived_zero_flat_survives(self):
        witness=archived_flat_witness()
        self.assertEqual(witness['restricted_rank'],1)
        self.assertEqual(witness['zero_outputs'],[0,1,2,3])
        for matrix in gl2():
            a,b,c,d=matrix
            self.assertEqual((mul(a,0)^mul(b,0),mul(c,0)^mul(d,0)),(0,0))

    def test_accounting_scope(self):
        cost=accounting()
        self.assertEqual(cost['postmix_gfni'],78)
        self.assertFalse(cost['native_measurement'])
        self.assertFalse(cost['full_certificate'])


if __name__=='__main__':
    unittest.main()
