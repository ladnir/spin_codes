"""Small exhaustive checks and exact regression for the window-box gate."""

from fractions import Fraction
import itertools
import math
import unittest

import rate_design_window_boxes as boxes


class WindowBoxTests(unittest.TestCase):
    def test_alternating_tail_against_subset_enumeration(self):
        for rows, window, regions in ((4, 2, 1), (4, 2, 3), (5, 3, 2), (6, 3, 3)):
            subsets = [set(value) for value in itertools.combinations(range(rows), window)]
            counts = [0]*(window+1)
            for selected in itertools.product(subsets, repeat=regions):
                counts[len(set.intersection(*selected))] += 1
            for minimum in range(1, window+1):
                exact = Fraction(sum(counts[minimum:]), len(subsets)**regions)
                for extra in range(window-minimum+2):
                    partial = boxes.tail_partial(rows, window, regions, minimum, extra)
                    if extra % 2:
                        self.assertLessEqual(partial, exact)
                    else:
                        self.assertGreaterEqual(partial, exact)
                self.assertEqual(boxes.tail_partial(rows, window, regions, minimum, window), exact)

    def test_shared_region_bounds_against_injective_assignments(self):
        for rows, window, overlap in ((4, 2, 1), (7, 3, 1)):
            first = set(range(window))
            windows = (first, set(range(window-overlap, 2*window-overlap)),
                       set(range(rows-window, rows)))
            for degree in (1, 2):
                for common in range(degree+1):
                    size = 2*degree-common
                    counts = [0, 0, 0]
                    for positions in itertools.permutations(range(rows), size):
                        left = positions[:degree]
                        right = positions[:common]+positions[degree:]
                        if all(position in first for position in left):
                            for category, second in enumerate(windows):
                                counts[category] += all(position in second for position in right)
                    same, adjacent, disjoint = boxes.shared_region_numerators(rows, window, overlap, degree, common)
                    self.assertEqual(same, counts[0])
                    self.assertGreaterEqual(adjacent, counts[1])
                    self.assertEqual(disjoint, counts[2])

    def test_all_box_pairs_against_full_small_permutations(self):
        rows, window, length, regions, minimum = 4, 2, 3, 2, 2
        windows = [set(range(start, start+window)) for start in range(3)]
        family = [(selected, choices)
                  for selected in itertools.combinations(range(length), regions)
                  for choices in itertools.product(range(3), repeat=regions)]
        permutations = list(itertools.permutations(range(rows)))
        exact_second = 0
        moments = [0, 0]
        for route in itertools.product(permutations, repeat=length):
            occupancies = [sum(all(route[region][row] in windows[choice]
                                   for region, choice in zip(selected, choices))
                               for row in range(rows)) for selected, choices in family]
            exact_second += math.comb(sum(value >= minimum for value in occupancies), 2)
            for degree in (1, 2):
                values = [boxes.choose(value, degree) for value in occupancies]
                moments[degree-1] += (sum(values)**2-sum(value*value for value in values))//2
        denominator = len(permutations)**length
        for degree in (1, 2):
            moment = Fraction(moments[degree-1], denominator*math.comb(minimum, degree)**2)
            upper = boxes.pair_upper(rows, window, 1, length, regions, minimum, degree)
            self.assertLessEqual(Fraction(exact_second, denominator), moment)
            self.assertLessEqual(moment, upper['total'])
            self.assertEqual(sum(upper['by_shared_regions']), upper['total'])

    def test_canonical_exact_rejection(self):
        result = boxes.window_box_lower()
        self.assertEqual(result['kernel_dimension_lower'], 1)
        self.assertEqual(result['box_count'], 552084530400)
        self.assertEqual(result['total_window'], 139680)
        self.assertLess(result['total_window'], result['bad_weight_threshold'])
        self.assertLessEqual(result['tail_lower'], result['tail_upper'])
        self.assertLess(result['tail_upper']-result['tail_lower'], result['tail_lower']/10**14)
        self.assertGreater(result['first_bonferroni_lower'], Fraction(63879, 10**16))
        self.assertLess(result['second_bonferroni_upper'], Fraction(28870, 10**20))
        self.assertEqual(result['failure_lower'], result['first_bonferroni_lower']-result['second_bonferroni_upper'])
        self.assertGreater(result['failure_lower'], Fraction(1, 1 << 38))
        self.assertTrue(result['exceeds_2_minus_40'])


if __name__ == '__main__':
    unittest.main()
