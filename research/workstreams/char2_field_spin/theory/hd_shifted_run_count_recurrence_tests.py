from fractions import Fraction as F
import unittest

from flint import arb, arb_poly, ctx

from hd_shifted_alias_certificate import rational
from hd_shifted_pairs_count import circular_run_count
from hd_shifted_run_potential import circular_run_counts, run_sum_arb
from hd_shifted_swap_potential import contraction_arb, full_counts


class RunCountRecurrenceTests(unittest.TestCase):
    def test_every_tiny_weight_and_run_count(self):
        for length in range(1, 41):
            for weight in range(length+1):
                actual = circular_run_counts(length, weight)
                expected = [circular_run_count(length, weight, runs)
                            for runs in range(weight+1)]
                self.assertEqual(actual, expected, (length, weight))

    def test_production_weights_exactly(self):
        for weight in (1280, 1536, 2048, 3072):
            actual = circular_run_counts(8192, weight)
            self.assertEqual(actual, [circular_run_count(8192, weight, runs)
                                      for runs in range(weight+1)])

    def test_constants_and_impossible_runs(self):
        self.assertEqual(circular_run_counts(8, 0), [1])
        self.assertEqual(circular_run_counts(8, 8), [1]+[0]*8)
        self.assertEqual(circular_run_counts(8, 7), [0, 8]+[0]*6)

    def test_directed_operation_order_is_unchanged(self):
        with ctx.workprec(128):
            for weight in range(9):
                counts = full_counts(4, weight)
                regions = {f: tuple(map(arb, (2+f, 1+2*f, 3+weight-f, 4+f)))
                           for f in counts}
                vector = (F(1), F(3, 7))
                coefficients = [arb(0) for _ in range(weight//2+1)]
                for full in counts:
                    coefficients[full] = (contraction_arb(regions[full], vector)**4).upper()
                convolution = arb_poly(coefficients)**2
                expected = arb(0)
                for degree in range(2*counts.start, 2*(counts.stop-1)+1):
                    multiplicity = circular_run_count(8, weight, weight-degree)
                    if multiplicity:
                        expected += multiplicity*convolution[degree].upper()
                expected = (expected*rational(vector[0]/min(vector))).upper()
                actual = run_sum_arb(4, weight, regions, vector, repetitions=8)
                self.assertEqual(actual, expected)

    def test_guards(self):
        for length, weight in ((0, 0), (8, -1), (8, 9), (8., 1), (True, 1)):
            with self.assertRaises(ValueError):
                circular_run_counts(length, weight)


if __name__ == '__main__':
    unittest.main()
