from fractions import Fraction as Q
import unittest

from flint import arb,ctx

from certificate_slack import output_tilt,upper_at,frontier


class CertificateSlackTests(unittest.TestCase):
    def test_all_witness_encodings_have_the_same_tilt(self):
        for witness in (['1/10',0,0,0],{'parameters':['1/10',0,0,0]},
                        {'old':{'base':{'parameters':['1/10',0,0,0]}}},
                        {'plane':{'tilt':'1/10'}}):
            self.assertEqual(output_tilt(witness),Q(1,10))
        with self.assertRaises(ValueError):output_tilt(['0'])

    def test_cutoff_search_matches_every_small_integer(self):
        ctx.prec=192
        terms=[(Q(1,3),arb(2)**-90),(Q(1,7),arb(2)**-70)]
        sparse=arb(2)**-60;maximum=200;target=arb(2)**-40
        expected=max(d for d in range(maximum+1) if sparse+upper_at(terms,d)<target)
        cutoff,dense,total=frontier(terms,sparse,maximum,40)
        self.assertEqual(cutoff,expected)
        self.assertTrue(total<target)
        self.assertTrue(sparse+upper_at(terms,cutoff+1)>target)
        self.assertTrue(total>=dense)
        self.assertEqual(frontier(terms,sparse,5,40)[0],5)
        with self.assertRaises(ArithmeticError):frontier(terms,arb(1),maximum,40)
        with self.assertRaises(ValueError):upper_at([(Q(0),arb(1))],10)
        with self.assertRaises(ValueError):upper_at(terms,-1)


if __name__=='__main__':unittest.main()
