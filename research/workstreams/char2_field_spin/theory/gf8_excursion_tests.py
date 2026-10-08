from fractions import Fraction as F
from itertools import product
import math
import unittest
import numpy as np
from gf8_excursion_screen import marked_local, path_mass


class ExcursionTests(unittest.TestCase):
    def test_every_formal_path_already_obeys_excursion_constraints(self):
        cases=0
        for weights in product(range(3),repeat=4):
            for following in product(range(2),repeat=4):
                previous=0
                births=returns=0
                legal=True
                for w,nxt in zip(weights,following):
                    if w==0 and previous!=nxt:
                        legal=False
                        break
                    if w and previous==0 and nxt==0:
                        legal=False
                        break
                    births += previous==0 and nxt==1
                    returns += previous==1 and nxt==0
                    previous=nxt
                if legal:
                    self.assertEqual(births-returns,following[-1])
                    self.assertLessEqual(births+returns,sum(weights))
                    cases+=1
        self.assertGreater(cases,100)

    def test_balance_mark_is_similarity_transform(self):
        y,tau=.97,6.
        original=marked_local(y,tau)
        amount=.3
        diagonal=np.diag([1.,math.exp(amount)])
        inverse=np.diag([1.,math.exp(-amount)])
        marked=marked_local(y,tau,'balance',amount)
        for w in range(9):
            np.testing.assert_allclose(marked[w],inverse@original[w]@diagonal,rtol=1e-14)

    def test_regional_balance_has_only_terminal_effect(self):
        point=(6.2,math.log(6))
        plain=path_mass(3,point,k_log=12,terminal=True)
        marked=path_mass(3,point,'balance',.3,k_log=12,terminal=True)
        np.testing.assert_allclose(marked,plain+[0.,.3],rtol=1e-12,atol=1e-11)


if __name__=='__main__':
    unittest.main()
