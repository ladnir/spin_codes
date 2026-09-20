import itertools
from fractions import Fraction as F
import unittest

import certify_refresh_q1 as q
from dual_track_q1 import even_binomial
from test_activation_q1 import actual_epochs, exact_coefficients, positive_regions


class RefreshTest(unittest.TestCase):
    def test_actual_gf4_encoder_is_dominated(self):
        for z in (F(1,3), F(1,2), F(3,4), F(1)):
            zero, one = q.epoch(4,2,{2:2,4:1},z,F)
            for epochs in (1,2,3):
                actual = exact_coefficients(*positive_regions(*actual_epochs(z),epochs),4)
                bound = q.coefficients(*q.region(zero,one,epochs,F),4,F)
                self.assertTrue(all(a <= b for a,b in zip(actual,bound)))

    def test_positive_epoch_power_and_support_average(self):
        zero,one=q.epoch(4,2,{2:2,4:1},F(2,3),F)
        for epochs in (1,2,3,7):
            rz,ra=q.region(zero,one,epochs,F)
            self.assertEqual((rz,ra),q.region(zero,one,epochs,F,linear=True))
            values=q.coefficients(rz,ra,4,F)
            for weight in range(5):
                total=F(0)
                supports=list(itertools.combinations(range(4),weight))
                for support in supports:
                    row=(F(1),F(0),F(0),F(0))
                    for region in range(4):
                        row=q.row_product(row,ra if region in support else rz)
                    total+=sum(row)
                self.assertEqual(values[weight],total/len(supports))

    def test_model_mass_and_symmetry(self):
        counts=even_binomial()
        self.assertEqual(sum(counts.values()),1<<128)
        self.assertEqual(counts[0],1)
        self.assertEqual(counts[256],1)
        self.assertTrue(all(counts[w]==counts[256-w] for w in counts))
        self.assertEqual(min(w for w in counts if w),38)


if __name__=='__main__':unittest.main()
