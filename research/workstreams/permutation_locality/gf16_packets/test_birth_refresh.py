import unittest
from fractions import Fraction as Q
from flint import arb,ctx
import birth_refresh as birth


def upper_q(value):
    m,e=value.upper().man_exp()
    return Q(int(m))*Q(2)**int(e)


class BirthOutwardTests(unittest.TestCase):
    def test_all_cutoffs_dominate_exact_weighted_feedback(self):
        ctx.prec=192;images=[sum(((a>>i)&1)*v for i,v in enumerate((1,6,120))) for a in range(8)]
        for columns in ([1,2,4,3,5,7,6,1],[0]*8,[1]*8):
            data=birth.prepare(images,columns,3)
            for p in (Q(0),Q(1,5),Q(15,16),Q(1)):
                z=Q(3,4);actual=[Q(0)]*8
                for x in range(256):
                    mass=(p/15 if x&15 else 1-p)*(p/15 if x>>4 else 1-p)
                    feedback=0
                    for b,c in enumerate(columns):
                        if x>>b&1:feedback^=c
                    actual[feedback]+=mass*z**x.bit_count()
                probabilities=birth.sc.probabilities(p)
                original=birth.base.outward_at_z(data,probabilities,arb(3)/4)
                for cutoff in (1,2,3):
                    matrix=birth.outward_candidate(data,p,arb(3)/4,cutoff,original)
                    self.assertLessEqual(actual[0],upper_q(matrix[0,0]))
                    residual=sum((max(Q(0),value-upper_q(matrix[0,2])/7) for value in actual[1:]),Q(0))
                    self.assertLessEqual(residual,upper_q(matrix[0,1]))
                    for i in (1,2):
                        for j in range(3):self.assertEqual(matrix[i,j],original[i,j])
                selected=birth.outward_at_z(data,probabilities,arb(3)/4)
                self.assertTrue(all(selected[i,j]>=0 for i in range(3) for j in range(3)))

    def test_invalid_cutoff(self):
        data=birth.prepare(list(range(8)),[1,2,4,3,5,7,6,1],3)
        original=birth.base.outward_at_z(data,birth.sc.probabilities(Q(1,5)),arb(3)/4)
        for cutoff in (0,4,1.5):
            with self.assertRaises(ValueError):birth.outward_candidate(data,Q(1,5),arb(3)/4,cutoff,original)


if __name__=='__main__':unittest.main()
