import unittest
import math
from fractions import Fraction as F
from collections import Counter
import test_general_occupancy as general_tests
import syndrome_activation as model
from test_activation_bridge import word


class SyndromeTests(unittest.TestCase):
    def test_exact_syndrome_probabilities(self):
        rows=[0xff,0xaa,0xcc,0xf0];t=8;s=4
        columns=[sum(((r>>j)&1)<<i for i,r in enumerate(rows)) for j in range(t)]
        spectrum=Counter(word(q,rows).bit_count() for q in range(1,16))
        counts=Counter()
        for x in range(1<<t):
            syndrome=0
            for j in range(t):
                if x>>j&1:syndrome^=columns[j]
            counts[x.bit_count(),syndrome]+=1
        kernel={j:counts[j,0] for j in range(t+1)}
        for j in range(t+1):
            rho=model.syndrome_density(t,s,spectrum,kernel,j)
            for syndrome in range(1,1<<s):
                self.assertLessEqual(F(counts[j,syndrome],math.comb(t,j)),rho/((1<<s)-1))


if __name__=='__main__':
    general_tests.model=model
    suite=unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromModule(general_tests),
                             unittest.defaultTestLoader.loadTestsFromTestCase(SyndromeTests)])
    result=unittest.TextTestRunner().run(suite)
    raise SystemExit(not result.wasSuccessful())
