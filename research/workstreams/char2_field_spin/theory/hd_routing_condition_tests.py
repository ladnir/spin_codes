"""Exact finite tests of the component-containment conditioning gate."""
import unittest
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, factorial

from hd_routing_condition import component_intervals, contained_probability, mgf_upper, receipt, interval_receipt


class RoutingConditionTests(unittest.TestCase):
    def test_components_are_disjoint_nonempty(self):
        supports = component_intervals([[1, 3], [0, 2]], 4)
        self.assertEqual(supports, [(1, 2), (3,), (4, 5), (6, 7)])
        crossing = component_intervals([[0, 2], [1, 3]], 4)
        self.assertEqual(crossing[1], (2, 3, 4))
        self.assertEqual(len({x for support in crossing for x in support}), 8)

    def test_exact_containment_and_mgf(self):
        # Cross-region components and arbitrary row placements are included.
        fixtures = [(3, [[0, 2], [1, 2]], [1, 1]),
                    (4, [[0, 2], [1, 3]], [2, 2]),
                    (4, [[1, 2], [0, 3]], [2, 1])]
        for rows, cuts, windows in fixtures:
            supports = component_intervals(cuts, rows)
            histogram = [0]*(len(supports)+1)
            joint = [0]*len(supports)
            total = factorial(rows)**len(cuts)
            for maps in product(permutations(range(rows)), repeat=len(cuts)):
                contained = [all(maps[x//rows][x % rows] < windows[x//rows]
                                 for x in support) for support in supports]
                histogram[sum(contained)] += 1
                for j, flag in enumerate(contained):
                    joint[j] += flag
            for j in range(len(supports)):
                self.assertEqual(F(joint[j], total),
                                 contained_probability(supports, [j], rows, windows))
            p = F(max(windows), rows)
            for degree in range(len(supports)+1):
                exact = sum(F(comb(h, degree)*count, total) for h, count in enumerate(histogram))
                selected = sum(contained_probability(supports, subset, rows, windows)
                               for subset in combinations(range(len(supports)), degree))
                self.assertEqual(exact, selected)
                self.assertLessEqual(exact, comb(len(supports), degree)*p**degree)
            tilt = F(3, 2)
            moment = sum(F(count, total)*tilt**h for h, count in enumerate(histogram))
            self.assertLessEqual(moment, mgf_upper(len(supports), p, tilt))

    def test_geometric_tail(self):
        for rows in range(3, 20):
            ratio_base = F(24, 25)**32
            start = rows//2+1
            first = comb(rows, start)*ratio_base**start
            ratio = F(rows-start, start+1)*ratio_base
            exact = sum(comb(rows, a)*ratio_base**a for a in range(start, rows+1))
            self.assertLessEqual(exact, first/(1-ratio))

    def test_run_set_count(self):
        for n in range(1, 12):
            counts = [0]*(n+1)
            for word in product((0, 1), repeat=n):
                runs = sum(bit and (i == 0 or not word[i-1]) for i, bit in enumerate(word))
                counts[runs] += 1
            for runs in range(1, (n+1)//2+1):
                self.assertEqual(counts[runs], comb(n+1, 2*runs))

    def test_start_mgf_unequal_regions(self):
        for rows, active in product(range(2, 7), (1, 2)):
            if active > rows:
                continue
            tilt = F(3, 2)
            for left, right in product(range(rows+1), repeat=2):
                moments = []
                for size in (left, right):
                    exact = sum(F(comb(size, hit)*comb(rows-size, active-hit), comb(rows, active))
                                *tilt**hit for hit in range(max(0, active-(rows-size)), min(size, active)+1))
                    moments.append(exact)
                average = F(left+right, 2*rows)
                self.assertLessEqual(moments[0]*moments[1], mgf_upper(2*active, average, tilt))

    def test_k20_gate(self):
        for precision in (128, 192):
            result = receipt(precision=precision)
            self.assertTrue(result['below_2_to_minus_100'])
            self.assertFalse(result['whole_distance_certificate'])
            general = interval_receipt(precision=precision)
            self.assertTrue(general['below_2_to_minus_100'])


if __name__ == '__main__':
    unittest.main()
