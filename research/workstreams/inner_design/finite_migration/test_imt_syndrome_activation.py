"""IMT syndrome tests; unique module name avoids historical test collisions."""
import itertools
import unittest
from fractions import Fraction
import numpy as np
import parameter_syndrome_activation as p


class SyndromeTests(unittest.TestCase):
    def test_transform_against_exact_supports(self):
        columns = [1,3,5,2,7,6]
        weights = p.image_weights(columns,3)
        for probability in (Fraction(0),Fraction(1,8),Fraction(1,2),Fraction(7,8),Fraction(1)):
            exact = [Fraction() for _ in range(8)]
            for bits in itertools.product((0,1),repeat=len(columns)):
                syndrome = 0
                for bit,column in zip(bits,columns):
                    if bit: syndrome ^= column
                j = sum(bits)
                exact[syndrome] += probability**j*(1-probability)**(len(columns)-j)
            got = p.syndrome_probabilities(weights,float(probability))
            np.testing.assert_allclose(got,list(map(float,exact)),atol=1e-15,rtol=0)

    def test_image_order(self):
        columns = [1,3,5,2,7,6]
        expected = [sum((q&c).bit_count()%2 for c in columns) for q in range(8)]
        self.assertEqual(p.image_weights(columns,3).tolist(),expected)


if __name__ == '__main__':
    unittest.main()
