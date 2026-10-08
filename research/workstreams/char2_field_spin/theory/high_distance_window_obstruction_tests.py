from fractions import Fraction as F
import unittest
from high_distance_window_obstruction import lower_bound, enumerate_small, screen


class WindowTests(unittest.TestCase):
    def test_exact_small_joint_enumeration(self):
        for b,k,r,w,s in ((2,1,2,1,0),(4,2,4,3,1),(4,2,4,2,0),(3,1,4,3,0)):
            actual_mean,actual_probability=enumerate_small(b,k,r,w,s)
            mean,lower=lower_bound(b,k,r,w,s)
            self.assertEqual(mean,actual_mean)
            self.assertLessEqual(lower,actual_probability)

    def test_k20_half_rate_obstruction(self):
        result=screen()
        self.assertEqual(result['rate'],'1/2')
        self.assertLess(F(result['window_fraction']),F(result['target_distance']))
        self.assertEqual(result['window']%4,0)
        self.assertGreater(F(result['failure_probability_lower_bound']),F(1,1024))
        self.assertLess(result['margin_upper_bits'],10)

    def test_bigger_outer_not_implicitly_certified(self):
        self.assertEqual(F(screen(256,128)['failure_probability_lower_bound']),0)
        # Zero means this event did not prove failure, not a distance proof.
        self.assertIsNone(screen(256,128)['margin_upper_bits'])


if __name__=='__main__':unittest.main()
