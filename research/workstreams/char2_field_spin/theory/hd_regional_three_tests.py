"""Bounded exact regressions for the regional-scan distance obstructions."""

from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb
import random
import unittest

from hd_regional_three_obstruction import (
    block_degree_ceiling, block_placement_failure_lower, cyclic_ceiling,
    regional_failure_lower, regional_prefix_ceiling,
    source_group_supports, window_count,
)


def matmul(a, b, prime):
    return [[sum(x * y for x, y in zip(row, col)) % prime
             for col in zip(*b)] for row in a]


def inverse_prefix(weights, prime):
    n = len(weights)
    d = [[0] * n for _ in range(n)]
    for i, weight in enumerate(weights):
        inv = pow(weight, -1, prime)
        d[i][i] = inv
        if i:
            d[i][i - 1] = -inv % prime
    return d


def prefix(weights, prime):
    return [[weights[j] if j <= i else 0 for j in range(len(weights))]
            for i in range(len(weights))]


class RegionalThreeTests(unittest.TestCase):
    def test_every_binary_window_bound(self):
        for n in range(2, 10):
            for bits in product((0, 1), repeat=n):
                a = sum(bits)
                for q in range(1, min(n, 5)):
                    z = window_count(bits, q)
                    self.assertGreaterEqual(z, max(0, (q + 1) * a - q * n - q))
                    self.assertLessEqual(z, max(0, a - q))

    def test_exact_uniform_window_mean_and_probability_lower(self):
        for n in range(3, 9):
            for a in range(1, n + 1):
                for q in range(1, min(n, 4)):
                    counts = []
                    for selected in combinations(range(n), a):
                        safe = set(selected)
                        counts.append(window_count([i in safe for i in range(n)], q))
                    expected = Fraction(sum(counts), len(counts))
                    formula = (n - q) * Fraction(comb(a, q + 1), comb(n, q + 1))
                    self.assertEqual(expected, formula)
                    maximum = max(0, a - q)
                    self.assertEqual(max(counts), maximum)
                    for threshold in range(maximum):
                        lower = max(Fraction(0),
                                    (expected - threshold) / (maximum - threshold))
                        actual = Fraction(sum(z > threshold for z in counts), len(counts))
                        self.assertGreaterEqual(actual, lower)

    def test_cyclic_cluster_containment_including_wrap(self):
        for rows in range(3, 10):
            for p in range(rows - 1):
                for a in range(p + 1, rows + 1):
                    for start in range(rows):
                        cluster = {(start + j) % rows for j in range(p + 1)}
                        count = sum(cluster <= {(left + j) % rows for j in range(a)}
                                    for left in range(rows))
                        self.assertGreaterEqual(count, a - p)

    def test_all_small_global_routes_cyclic_average(self):
        b, rows, p, q = 2, 3, 1, 1
        for route in permutations(range(b * rows)):
            supports = source_group_supports(b, rows, p, q, route)
            for a in (2, 3):
                total = sum(sum(support <= {(left + j) % rows for j in range(a)}
                                for support in supports)
                            for left in range(rows))
                self.assertGreaterEqual(Fraction(total, rows),
                                        b * ((q + 1) * (a - p) - q * rows))

    def test_weighted_matrix_support_and_global_carry(self):
        rng = random.Random(12391)
        b, rows, p, q, prime = 4, 9, 2, 3, 17
        n = b * rows
        identity = [[int(i == j) for j in range(n)] for i in range(n)]
        for _ in range(8):
            weights = [[rng.randrange(1, prime) for _ in range(n)] for _ in range(p + q)]
            for w in weights:
                self.assertEqual(matmul(inverse_prefix(w, prime), prefix(w, prime), prime),
                                 identity)
            dpre = identity
            for w in weights[:p]:
                dpre = matmul(dpre, inverse_prefix(w, prime), prime)
            dpost = identity
            for w in weights[p:]:
                dpost = matmul(dpost, inverse_prefix(w, prime), prime)
            route = list(range(n))
            rng.shuffle(route)
            pinv = [[int(i == route[j]) for j in range(n)] for i in range(n)]
            exact = matmul(matmul(dpre, pinv, prime), dpost, prime)
            conservative = source_group_supports(b, rows, p, q, route)
            for col in range(n):
                actual = {row % rows for row in range(n) if exact[row][col]}
                self.assertLessEqual(actual, conservative[col])
            # The global inverse really crosses a coordinate-region boundary.
            self.assertNotEqual(inverse_prefix(weights[0], prime)[rows][rows - 1], 0)
            self.assertNotEqual(inverse_prefix(weights[-1], prime)[rows][rows - 1], 0)

    def test_fixed_prefix_every_sampled_regional_route(self):
        rng = random.Random(80331)
        b, rows, p, q = 4, 40, 2, 1
        bound = regional_prefix_ceiling(b, rows, p, q)
        safe_groups = set(range(bound["a"]))
        for _ in range(100):
            route = []
            for region in range(b):
                labels = list(range(region * rows, (region + 1) * rows))
                rng.shuffle(labels)
                route.extend(labels)
            supports = source_group_supports(b, rows, p, q, route)
            # Count only within-region post windows, just as the proof does.
            confined = sum(supports[j] <= safe_groups for j in range(b * rows)
                           if j % rows < rows - q)
            self.assertGreaterEqual(confined, bound["confined_lower"])
            self.assertGreater(confined, bound["capacity"])

    def test_k20_exact_finite_conclusions(self):
        b, rows, p = 64, 32768, 2
        prefix_bound = regional_prefix_ceiling(b, rows, p, 1)
        cyclic_bound = cyclic_ceiling(b, rows, p, 1)
        self.assertEqual(prefix_bound["a"], 21849)
        self.assertEqual(prefix_bound["confined_lower"], 699200)
        self.assertEqual(prefix_bound["distance_upper"], 699169)
        self.assertEqual(cyclic_bound["distance_upper"], 699169)
        self.assertEqual(cyclic_bound["confined_lower"], 699264)
        for q, lower_decimal in ((1, "0.59980"), (2, "0.27964"), (3, "0.02350")):
            result = regional_failure_lower(b, rows, p, q)
            self.assertEqual(result["a"], 26214)
            self.assertEqual(result["capacity"], 838848)
            self.assertEqual(result["distance_upper_on_event"], 838849)
            self.assertLessEqual(result["distance_upper_on_event"], result["bad_max"])
            self.assertGreater(result["failure_lower"], Fraction(lower_decimal))
            self.assertGreater(result["failure_lower"], Fraction(1, 1 << 40))
        self.assertEqual(regional_failure_lower(b, rows, p, 4)["failure_lower"], 0)

    def test_block_preserving_placement_gate(self):
        rng = random.Random(71892)
        rows, block, p, selected = 32, 8, 2, {0, 2}
        safe_groups = {j for k in selected for j in range(k * block, (k + 1) * block)}
        for _ in range(100):
            block_order = list(range(rows // block))
            rng.shuffle(block_order)
            placement = []
            for k in block_order:
                labels = list(range(k * block, (k + 1) * block))
                rng.shuffle(labels)
                placement.extend(labels)
            # Dropping region-crossing columns makes this count conservative.
            safe = sum(set(placement[j:j + p + 1]) <= safe_groups
                       for j in range(rows - p))
            self.assertGreaterEqual(safe, len(selected) * (block - p))
        result = block_placement_failure_lower(64, 32768, 16)
        self.assertEqual(result["a"], 26208)
        self.assertEqual(result["safe"], 22932)
        self.assertEqual(result["distance_upper_on_event"], 838657)
        self.assertGreater(result["failure_lower"], Fraction("0.29956"))
        self.assertEqual(block_placement_failure_lower(64, 32768, 8)["failure_lower"], 0)

    def test_grouped_degree_including_global_region_crossings(self):
        rng = random.Random(100938)
        b, rows, q = 4, 32, 1
        for block in (2, 4, 8):
            for p in (2, block, block + 1):
                placements = []
                for _ in range(b):
                    block_order = list(range(rows // block))
                    rng.shuffle(block_order)
                    for k in block_order:
                        labels = list(range(k * block, (k + 1) * block))
                        rng.shuffle(labels)
                        placements.extend(labels)
                route = list(range(b * rows))
                rng.shuffle(route)  # Arbitrary global route, not just regional.
                degree = ((p + block - 1) // block + 1) * (q + 1)
                for col in range(b * rows):
                    touched = set()
                    for j in range(col, min(b * rows, col + q + 1)):
                        middle = route[j]
                        for source in range(middle, min(b * rows, middle + p + 1)):
                            touched.add(placements[source] // block)
                    self.assertLessEqual(len(touched), degree)
        expected = {2: 832321, 4: 832385, 8: 832513, 16: 832513,
                    32: 833537, 64: 833537, 128: 835585}
        for block, distance in expected.items():
            result = block_degree_ceiling(64, 32768, block)
            self.assertEqual(result["degree"], 4)
            self.assertEqual(result["distance_upper"], distance)
            self.assertLess(distance, 838861)
        self.assertEqual(block_degree_ceiling(64, 32768, 4, 4)["degree"], 4)
        self.assertEqual(block_degree_ceiling(64, 32768, 4, 5)["degree"], 6)


if __name__ == "__main__":
    unittest.main()
