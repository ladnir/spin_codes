import unittest
from gf8_random_reference import finite_margin,reference_rate,smallest_sufficient_length


class RandomReferenceTests(unittest.TestCase):
    def test_reference(self):
        self.assertAlmostEqual(reference_rate(),.8414461328406041,places=14)

    def test_exact_small_length_search(self):
        for k in (1,17,100):
            result=smallest_sufficient_length(k,40)
            n=result['n']
            self.assertGreaterEqual(finite_margin(k,n).lower(),40)
            for previous in range(k,n):
                self.assertLess(finite_margin(k,previous).upper(),40)


if __name__=='__main__':unittest.main()
