"""Independent exhaustive toy check of the three-state GL envelope."""

from fractions import Fraction as F
import unittest


class SharedInnerTransferTest(unittest.TestCase):
    def test_exact_gl2_envelope(self):
        # A(q0,q1)=(q0,q1,q0+q1), with C=A^T.
        states = (1, 2, 3)
        maps = [(a, b) for a in states for b in states if a != b]
        self.assertEqual(len(maps), 6)
        beta, z = F(23, 100), F(61, 100)

        def expand(q):
            a, b = q & 1, (q >> 1) & 1
            return a | (b << 1) | ((a ^ b) << 2)

        def feedback(x):
            return ((x & 1) ^ ((x >> 2) & 1)) | (
                (((x >> 1) & 1) ^ ((x >> 2) & 1)) << 1)

        def apply(matrix, q):
            return (matrix[0] if q & 1 else 0) ^ (
                matrix[1] if q & 2 else 0)

        input_probability = {
            x: beta**x.bit_count() * (1 - beta)**(3 - x.bit_count())
            for x in range(8)
        }
        exact = [[F(0) for _ in range(4)] for _ in range(4)]
        for q in range(4):
            for x, probability in input_probability.items():
                moment = probability * z**(x ^ expand(q)).bit_count()
                for matrix in maps:
                    exact[q][apply(matrix, q) ^ feedback(x)] += moment / len(maps)

        g0 = 1 - beta + beta * z
        g1 = beta + (1 - beta) * z
        live_moments = [g0**(3 - expand(q).bit_count())
                        * g1**expand(q).bit_count() for q in states]
        fresh = max(live_moments)
        uniform = sum(live_moments, F(0)) / len(states)
        zero_zero = sum((probability * z**x.bit_count()
                         for x, probability in input_probability.items()
                         if feedback(x) == 0), F(0))
        active = sum((probability for x, probability in input_probability.items()
                      if feedback(x) != 0), F(0))

        self.assertLessEqual(exact[0][0], zero_zero)
        self.assertLessEqual(sum(exact[0][1:]), g0**3 - zero_zero)
        for q in states:
            self.assertLessEqual(exact[q][0], min(fresh, active) / 3)
            for target in states:
                self.assertLessEqual(exact[q][target], fresh / 3)
        self.assertLessEqual(sum(exact[q][0] for q in states) / 3,
                             min(uniform, active) / 3)
        for target in states:
            self.assertLessEqual(sum(exact[q][target] for q in states) / 3,
                                 uniform / 3)


if __name__ == "__main__":
    unittest.main()
