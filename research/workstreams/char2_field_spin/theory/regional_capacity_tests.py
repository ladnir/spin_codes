"""Exact and exhaustive tests for the generalized structural capacity gates."""

from fractions import Fraction
import itertools
import math
import unittest

import rate_design_window_boxes as original
import regional_capacity_screen as capacity


class RegionalCapacityTests(unittest.TestCase):
    def test_kernel_geometry_and_suffix_saving(self):
        for b, k, state in ((33, 24, 4), (32, 23, 4), (32, 22, 4), (34, 24, 4), (32, 24, 5), (32, 24, 6)):
            for log_size in (16, 18, 20):
                for h in range(b-k+1, b+1):
                    for suffix in (False, True):
                        result = capacity.geometry(log_size, b, k, state, h, 3, suffix)
                        a, d, returns = result['minimum'], h-(b-k), h-int(suffix)
                        self.assertGreater(a*d-state*returns, 0)
                        self.assertLessEqual((a-1)*d-state*returns, 0)
                        self.assertEqual(result['rows'] % 64, 0)
                        self.assertLessEqual(result['total_window'], result['threshold'])
                        self.assertEqual(result['starts'][-1]+result['window'], result['rows'])
                        self.assertTrue(all(start % 16 == 0 for start in result['starts']))

    def test_window_assignments_exact(self):
        for size, width, overlap in ((4, 2, 1), (6, 3, 1), (7, 3, 0)):
            left, right = set(range(width)), set(range(width-overlap, 2*width-overlap))
            for degree in (1, 2):
                for common in range(degree+1):
                    count = 0
                    for positions in itertools.permutations(range(size), 2*degree-common):
                        a, b = positions[:degree], positions[:common]+positions[degree:]
                        count += all(x in left for x in a) and all(x in right for x in b)
                    self.assertEqual(capacity.joint_window_assignments(width, overlap, degree, common), count)

    def test_original_coarse_bound_reproduced(self):
        record = capacity.geometry(20, 32, 24, 4, 9, 3)
        record.update(window=15520, starts=(0, 14096, 28192))
        old = original.pair_upper(43712, 15520, 1424, 32, 9, 37, 29)
        new = capacity.mixed_pair_upper(record, 29, exact_overlap=False)
        self.assertEqual(new['pair_upper'], old['total'])
        self.assertLessEqual(capacity.mixed_pair_upper(record, 29)['pair_upper'], old['total'])

    def test_suffix_and_checksum_pairs_by_full_enumeration(self):
        rows, length, h, width, minimum = 4, 3, 2, 2, 2
        starts, count = (0, 1, 2), 3
        windows = [set(range(start, start+width)) for start in starts]
        permutations = list(itertools.permutations(range(rows)))
        for suffix, checksum in itertools.product((False, True), repeat=2):
            free = h-int(suffix)
            record = dict(rows=rows, length=length, regions=h, free_regions=free,
                          window=width, minimum=minimum, starts=starts, suffix=suffix, checksum=checksum)
            family = []
            for selected in itertools.combinations(range(length-int(suffix)), free):
                for choices in itertools.product(range(count), repeat=free):
                    if checksum and sum(choices) % count:
                        continue
                    family.append((selected+((length-1,) if suffix else ()), choices+((count-1,) if suffix else ())))
            self.assertEqual(len(family), math.comb(length-int(suffix), free)*count**(free-int(checksum)))
            moments = [0, 0]
            for route in itertools.product(permutations, repeat=length):
                occupancies = [sum(all(route[region][row] in windows[choice]
                                       for region, choice in zip(selected, choices))
                                   for row in range(rows)) for selected, choices in family]
                for degree in (1, 2):
                    values = [capacity.choose(value, degree) for value in occupancies]
                    moments[degree-1] += (sum(values)**2-sum(value*value for value in values))//2
            for degree in (1, 2):
                exact = Fraction(moments[degree-1], len(permutations)**length*math.comb(minimum, degree)**2)
                upper = capacity.mixed_pair_upper(record, degree)['pair_upper']
                self.assertLessEqual(exact, upper)
                if not checksum:
                    self.assertEqual(exact, upper)

    def test_sliding_interval_probability(self):
        for rows in range(2, 9):
            for width in range(1, rows+1):
                for minimum in range(1, width+1):
                    record = dict(rows=rows, window=width, minimum=minimum, regions=2, length=3, suffix=False)
                    exact = Fraction(sum(max(subset)-min(subset) < width
                                         for subset in itertools.combinations(range(rows), minimum)),
                                     math.comb(rows, minimum))
                    self.assertEqual(capacity.sliding_event_upper(record)['sliding_probability'], exact)

    def test_unequal_spans_against_exhaustive_subsets(self):
        for rows, minimum, regions in ((4, 2, 2), (5, 2, 3), (6, 3, 2)):
            subsets = list(itertools.combinations(range(rows), minimum))
            for suffix in (False, True):
                histogram = [0]*(rows*regions+1)
                for positions in itertools.product(subsets, repeat=regions):
                    total = sum(max(part)-min(part)+1 for part in positions[:regions-int(suffix)])
                    if suffix:
                        total += rows-min(positions[-1])
                    histogram[total] += 1
                running = 0
                for threshold, count in enumerate(histogram):
                    running += count
                    record = dict(rows=rows, minimum=minimum, regions=regions, length=regions,
                                  threshold=threshold, suffix=suffix)
                    bound = capacity.unequal_window_upper(record)['unclipped_first_moment']/len(subsets)
                    self.assertGreaterEqual(bound, Fraction(running, len(subsets)**regions))

    def test_span_generating_coefficient(self):
        for minimum, regions, suffix in itertools.product((2, 3), (2, 3), (False, True)):
            for threshold in range(regions*minimum, regions*minimum+5):
                total = 0
                for spans in itertools.product(range(minimum, threshold+1), repeat=regions):
                    if sum(spans) <= threshold:
                        total += math.prod(capacity.choose(span-2, minimum-2) for span in spans[:regions-int(suffix)])*\
                                 (capacity.choose(spans[-1]-1, minimum-1) if suffix else 1)
                self.assertEqual(total, capacity.choose(threshold-regions+int(suffix), regions*(minimum-1)+int(suffix)))

    def test_finite_comparison_and_scope(self):
        for b, k, ceiling in ((33, 24, 58), (32, 23, 65)):
            result = capacity.all_support_sliding_upper(20, b, k, 4)
            self.assertLess(result['event_union_upper'], Fraction(1, 1 << ceiling))
            self.assertIn('equal-length-window-kernel-events-only', result['scope'])
        for b, k, floor_bits in ((33, 24, 19), (32, 23, 27), (32, 22, 90), (34, 24, 74)):
            unequal = capacity.all_support_unequal_upper(20, b, k, 4)
            self.assertLess(unequal['event_union_upper'], Fraction(1, 1 << floor_bits))
            self.assertGreater(unequal['event_union_upper'], Fraction(1, 1 << (floor_bits+1)))
            self.assertIn('unequal-window-kernel-events-only', unequal['scope'])
        strong = capacity.exact_check(capacity.geometry(18, 28, 21, 4, 8, 5, True))
        self.assertTrue(strong['rejects_40_bits'])
        near = capacity.exact_check(capacity.geometry(20, 32, 24, 5, 9, 14, True, checksum=True), degrees=[30])
        self.assertGreater(near['failure_lower'], Fraction(1, 1 << 41))
        self.assertLess(near['failure_lower'], Fraction(1, 1 << 40))
        self.assertFalse(near['rejects_40_bits'])


if __name__ == '__main__':
    unittest.main()
