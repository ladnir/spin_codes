"""Exhaustive small analogue for one omitted mixer and two-mixer obstruction.

Two parallel GF4 RS[4,2] rows have16-bit outputs and8-bit messages. These
checks enumerate all256messages; no random sampling supports the lemma.
"""
from collections import Counter
from fractions import Fraction as F
import unittest


def mul(a,b):
    value=0
    for _ in range(2):
        if b&1:value^=a
        b>>=1;a<<=1
        if a&4:a^=7
    return value


def words():
    for message in range(1,256):
        rows=[[(message>>(4*r+2*i))&3 for i in range(2)] for r in range(2)]
        code=[[mul(x^1,row[0])^mul(x,row[1]) for x in range(4)] for row in rows]
        yield tuple(code[0][i]|(code[1][i]<<2) for i in range(4))


class OmissionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.code=list(words())
        cls.support=Counter(tuple(i for i,x in enumerate(w) if x) for w in cls.code)

    def test_parallel_code_mds_spectrum(self):
        weights=Counter(sum(bool(x) for x in w) for w in self.code)
        self.assertEqual(weights,{3:60,4:195})
        self.assertEqual(len(set(self.code)),255)

    def test_all_single_omissions_preserve_expected_pointwise_measure(self):
        for omitted in range(4):
            joint=Counter((tuple(i for i,x in enumerate(w) if x),w[omitted]) for w in self.code)
            for support,count in self.support.items():
                h=len(support)
                for a in range(1,16) if omitted in support else (0,):
                    random_positions=h-int(omitted in support)
                    omitted_density=F(joint[support,a],15**random_positions)
                    full_density=F(count,15**h)
                    self.assertEqual(omitted_density,full_density)

    def test_two_omissions_do_not_preserve_joint_measure(self):
        support=(0,1,2)
        joint=Counter((w[0],w[1]) for w in self.code
            if tuple(i for i,x in enumerate(w) if x)==support)
        self.assertEqual(len(joint),15)
        self.assertEqual(set(joint.values()),{1})
        self.assertEqual(F(max(joint.values()),15),15*F(self.support[support],15**3))


if __name__=='__main__':unittest.main()
