"""Independent exact small-support checks; no numerical tolerances."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product, accumulate, permutations
from math import comb, prod
import unittest

from support import (union_shells, lowest_weight_spectrum, union_cdf_upper,
                     fixed_weight_union, expected_union, weighted_union_shells,
                     weighted_cdf_upper, integer_cdf, tilted_union_shells,
                     tilted_support_caps)


def enumerate_supports(spectra, full_weight=Q(1)):
    n = len(spectra[0])-1
    choices = [[(mask, Q(row[mask.bit_count()], comb(n, mask.bit_count())))
                for mask in range(1 << n) if row[mask.bit_count()]]
               for row in spectra]
    result = [Q(0)]*(n+1)
    for selected in product(*choices):
        union = 0
        intersection = (1 << n)-1
        for mask, _ in selected:
            union |= mask
            intersection &= mask
        result[union.bit_count()] += prod(value for _, value in selected)*full_weight**intersection.bit_count()
    return result


class SupportTests(unittest.TestCase):
    def test_upward_integer_cdf(self):
        for values in ([Q(0),Q(1,3),Q(2,3),Q(4,3)], [Q(1),Q(7,2),Q(9)]):
            rounded=integer_cdf(values)
            self.assertTrue(all(a<=b<a+1 for a,b in zip(values,rounded)))
        with self.assertRaises(ValueError):
            integer_cdf([Q(2),Q(1)])

    def test_weighted_support_enumeration(self):
        for spectra in ([[1,0,3,0,1]]*4,
                        [[1,2,1,0,0], [1,0,0,2,1]],
                        [[0,3,0,0], [0,0,5,0], [1,0,0,0]],
                        [[0,0], [1,1]], [[1]]*4):
            self.assertEqual(union_shells(spectra), enumerate_supports(spectra))

    def test_all_weight_triples(self):
        for n in range(1,5):
            for weights in product(range(n+1), repeat=3):
                spectra = [[int(v == w) for v in range(n+1)] for w in weights]
                actual = fixed_weight_union(n, weights)
                self.assertEqual(actual, enumerate_supports(spectra))
                self.assertEqual(sum(u*p for u,p in enumerate(actual)),
                                 expected_union(n, weights))

    def test_identical_words_are_not_shared_supports(self):
        for n in range(1,9):
            for w in range(n+1):
                distribution = fixed_weight_union(n, [w]*4)
                self.assertEqual(distribution[w], Q(1, comb(n,w)**3))
                if 0 < w < n:
                    self.assertLess(distribution[w], 1)

    def test_shell_cap_monotonicity(self):
        exact = [1,0,3,0,1]
        caps = [1,2,5,3,1]
        true_shells = union_shells([exact]*4)
        upper_shells = union_shells([caps]*4)
        self.assertTrue(all(a <= b for a,b in zip(true_shells, upper_shells)))

    def test_total_constrained_cdf(self):
        caps = [1,2,2,2,1]
        total = 5
        upper = union_cdf_upper(caps, total)
        cases = 0
        for tail in product(*(range(c+1) for c in caps[1:])):
            exact = [1,*tail]
            if sum(exact) != total:
                continue
            shells = union_shells([exact]*4)
            shells[0] -= 1
            self.assertTrue(all(a <= b for a,b in zip(accumulate(shells), upper)))
            cases += 1
        self.assertGreater(cases, 10)
        self.assertEqual(sum(lowest_weight_spectrum(caps, total)), total)

    def test_shared_and_independent_differ(self):
        # A tiny linear code. Sharing a coordinate permutation preserves
        # each tuple's original union; independent permutations need not.
        words = [0, 3, 12, 15]
        spectrum = [sum(x.bit_count()==w for x in words) for w in range(5)]
        shared = Counter((a|b|c|d).bit_count() for a,b,c,d in product(words,repeat=4))
        independent = union_shells([spectrum]*4)
        self.assertNotEqual(independent, [shared[u] for u in range(5)])
        self.assertEqual(sum(independent), len(words)**4)

    def test_permutation_invariant_code_average(self):
        # If the code contains every support in each of its shells, the
        # tuple-averaged counts coincide even though fixed tuples differ.
        words = [0, 3, 5, 6]
        spectrum = [sum(x.bit_count()==w for x in words) for w in range(4)]
        shared = Counter((a|b|c|d).bit_count() for a,b,c,d in product(words,repeat=4))
        self.assertEqual(union_shells([spectrum]*4), [shared[u] for u in range(4)])

    def test_weighted_enumeration(self):
        for spectrum in ([1,0,3,0,1], [1,2,0,1], [0,0,1]):
            for weight in (Q(1), Q(4,3), Q(2)):
                self.assertEqual(weighted_union_shells(spectrum,4,weight),
                                 enumerate_supports([spectrum]*4,weight))

    def test_weighted_cdf_cap(self):
        exact, caps = [1,0,3,0,1], [1,2,5,3,1]
        for weight in (Q(1), Q(4,3), Q(2)):
            shells = weighted_union_shells(exact,4,weight)
            shells[0] -= 1
            upper = weighted_cdf_upper(caps,sum(exact),4,weight)
            self.assertTrue(all(a <= b for a,b in zip(accumulate(shells),upper)))
            shell_caps = weighted_union_shells(caps,4,weight)
            shell_caps[0] -= 1
            self.assertTrue(all(a <= b for a,b in zip(shells,shell_caps)))

    def test_conditional_histogram_exchangeability(self):
        for n,weights in ((3,(1,1,2,2)), (4,(0,1,2,4)), (4,(1,1,2,2))):
            choices = [[x for x in range(1 << n) if x.bit_count()==w] for w in weights]
            sequences = Counter(tuple(sum((x >> c)&1 for x in rows) for c in range(n))
                                for rows in product(*choices))
            histograms = {tuple(sorted(sequence)) for sequence in sequences}
            for histogram in histograms:
                # Conditional on the histogram, all column-weight sequences
                # occur equally often, including every ordering.
                masses = {sequences[sequence] for sequence in set(permutations(histogram))}
                self.assertEqual(len(masses),1)
                self.assertGreater(next(iter(masses)),0)

    def test_input_tilt_exact_enumeration(self):
        for spectrum,rows in (([1,1,0,1],4), ([1,0,2,1],3),
                              ([Q(1,2),Q(2,3),Q(1,5)],2), ([1],4)):
            for input_weight in (Q(1,2),Q(49,50),Q(1),Q(5,4),Q(2)):
                tilted = [a/input_weight**w for w,a in enumerate(spectrum)]
                for full_weight in (Q(1),Q(4,3),Q(8)):
                    actual = enumerate_supports([tilted]*rows,full_weight)
                    self.assertEqual(tilted_union_shells(spectrum,rows,full_weight,input_weight),actual)

    def test_input_tilt_caps_exhaustive_small_spectra(self):
        caps,total = [1,2,2,1],4
        actual_spectra = [[1,*tail] for tail in product(*(range(c+1) for c in caps[1:]))
                          if 1+sum(tail)==total]
        self.assertGreater(len(actual_spectra),3)
        for rows in (2,3):
            for input_weight in (Q(1,2),Q(49,50),Q(1),Q(5,4),Q(2)):
                for full_weight in (Q(1),Q(4,3),Q(8)):
                    shells,cdf = tilted_support_caps(caps,total,rows,full_weight,input_weight)
                    for spectrum in actual_spectra:
                        tilted = [a/input_weight**w for w,a in enumerate(spectrum)]
                        actual = enumerate_supports([tilted]*rows,full_weight)
                        actual[0] -= 1
                        self.assertTrue(all(a<=b for a,b in zip(actual,shells)))
                        self.assertTrue(all(a<=b for a,b in zip(accumulate(actual),cdf)))

    def test_input_tilt_identity_preserves_old_bounds(self):
        caps,total = [1,2,5,3,1],5
        for rows in (1,2,4):
            for full_weight in (Q(1),Q(4,3),Q(2)):
                shells,cdf = tilted_support_caps(caps,total,rows,full_weight,Q(1))
                old_shells = weighted_union_shells(caps,rows,full_weight)
                old_shells[0] -= 1
                self.assertEqual(shells,old_shells)
                self.assertEqual(cdf,weighted_cdf_upper(caps,total,rows,full_weight))

    def test_joint_tilt_is_not_greedy_monotone(self):
        # Lowering row weights increases a^-W but can decrease gamma^J.
        # Directly using the greedy spectrum for the joint factor is unsafe.
        actual,caps,total = [1,0,1],[1,1,1],2
        rows,full_weight,input_weight = 2,Q(8),Q(2)
        greedy = lowest_weight_spectrum(caps,total)
        actual_shells = tilted_union_shells(actual,rows,full_weight,input_weight)
        greedy_shells = tilted_union_shells(greedy,rows,full_weight,input_weight)
        self.assertEqual(sum(actual_shells)-1,Q(9,2))
        self.assertEqual(sum(greedy_shells)-1,Q(17,8))
        self.assertGreater(sum(actual_shells),sum(greedy_shells))
        shells,cdf = tilted_support_caps(caps,total,rows,full_weight,input_weight)
        actual_shells[0] -= 1
        self.assertTrue(all(a<=b for a,b in zip(actual_shells,shells)))
        self.assertTrue(all(a<=b for a,b in zip(accumulate(actual_shells),cdf)))

    def test_input_tilt_validation(self):
        for kwargs in ({'input_weight':Q(0)}, {'input_weight':Q(-1)},
                       {'full_weight':Q(1,2)}, {'rows':0}):
            with self.assertRaises(ValueError):
                tilted_union_shells([1,1],**kwargs)
            with self.assertRaises(ValueError):
                tilted_support_caps([1,1],2,**kwargs)
        for caps,total in (([2,1],3), ([1,-1],1), ([1,1],3), ([1,Q(1,2)],1)):
            with self.assertRaises(ValueError):
                tilted_support_caps(caps,total)

    def test_input_tilt_below_one_breaks_greedy_coupling(self):
        actual,caps,total = [1,0,1],[1,1,1],2
        greedy = lowest_weight_spectrum(caps,total)
        actual_shells = tilted_union_shells(actual,2,Q(1),Q(1,2))
        greedy_shells = tilted_union_shells(greedy,2,Q(1),Q(1,2))
        self.assertEqual(sum(actual_shells)-1,24)
        self.assertEqual(sum(greedy_shells)-1,8)
        self.assertGreater(sum(actual_shells),sum(greedy_shells))
        shells,cdf = tilted_support_caps(caps,total,2,Q(1),Q(1,2))
        actual_shells[0] -= 1
        self.assertTrue(all(a<=b for a,b in zip(actual_shells,shells)))
        self.assertTrue(all(a<=b for a,b in zip(accumulate(actual_shells),cdf)))


if __name__ == '__main__':
    unittest.main()
