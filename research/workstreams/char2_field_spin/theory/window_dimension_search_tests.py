"""Exact objective, monotonicity, tiny exhaustive, and confinement checks."""

from itertools import product
import unittest

import numpy as np

from window_dimension_search import (generate_route, validate_route, count_hits, dimension,
                                     aligned_scores, coordinate_ascent, search_route,
                                     strict_window_geometry, exchange_lengths)


def kernel_basis_prime(matrix, columns, prime=5):
    rows = [list(row) for row in matrix]
    pivots = []
    for column in range(columns):
        current = len(pivots)
        pivot = next((i for i in range(current, len(rows)) if rows[i][column] % prime), None)
        if pivot is None:
            continue
        rows[current], rows[pivot] = rows[pivot], rows[current]
        inverse = pow(int(rows[current][column]), -1, prime)
        rows[current] = [(inverse*x) % prime for x in rows[current]]
        for i in range(len(rows)):
            if i != current:
                factor = rows[i][column]
                rows[i] = [(x-factor*y) % prime for x, y in zip(rows[i], rows[current])]
        pivots.append(column)
    basis = []
    for free in range(columns):
        if free in pivots:
            continue
        vector = [0]*columns
        vector[free] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = (-rows[row][free]) % prime
        basis.append(vector)
    return basis


class WindowDimensionTests(unittest.TestCase):
    def test_strict_ten_percent_geometry(self):
        width, threshold = strict_window_geometry(47680, 32, 16)
        self.assertEqual((width, threshold), (4752, 152575))
        self.assertEqual(width*32, 152064)
        self.assertGreater((width+16)*32, threshold)

    def test_count_objective_and_score_exactness(self):
        route = generate_route(12, 5, 501)
        starts = np.asarray((0, 2, 4, 6, 8))
        counts = count_hits(route, starts, 4)
        brute = [sum(i in region[first:first+4] for region, first in zip(route, starts)) for i in range(12)]
        self.assertEqual(counts.tolist(), brute)
        self.assertEqual(dimension(counts, 2), sum(max(c-2, 0) for c in brute))
        for region_index in range(5):
            without = counts.copy()
            without[route[region_index, starts[region_index]:starts[region_index]+4]] -= 1
            scores = aligned_scores(route[region_index], without, 2, 4, 2)
            for choice, first in enumerate(range(0, 9, 2)):
                candidate = without.copy()
                candidate[route[region_index, first:first+4]] += 1
                self.assertEqual(int(scores[choice]), dimension(candidate, 2)-dimension(without, 2))

    def test_tiny_exhaustive_coordinate_optima(self):
        route = generate_route(6, 3, 702)
        starts_list = list(product(range(5), repeat=3))
        values = [dimension(count_hits(route, starts, 2), 1) for starts in starts_list]
        maximum = max(values)
        for initial in starts_list:
            result = coordinate_ascent(route, 2, 1, initial, alignment=1, max_sweeps=20)
            self.assertTrue(result.converged)
            self.assertEqual(sorted(result.history), list(result.history))
            self.assertLessEqual(result.dimension, maximum)
            for coordinate in range(3):
                for first in range(5):
                    candidate = result.starts.copy()
                    candidate[coordinate] = first
                    self.assertLessEqual(dimension(count_hits(route, candidate, 2), 1), result.dimension)
        optimal = starts_list[values.index(maximum)]
        self.assertEqual(coordinate_ascent(route, 2, 1, optimal, alignment=1).dimension, maximum)

    def test_ties_preserve_current_starts(self):
        route = generate_route(12, 3, 47)
        initial = np.asarray((0, 2, 4))
        result = coordinate_ascent(route, 4, 4, initial, alignment=2)
        self.assertEqual(result.dimension, 0)
        self.assertEqual(result.coordinate_updates, 0)
        self.assertTrue(np.array_equal(result.starts, initial))

    def test_restarts_are_bounded_replayable_and_verified(self):
        route = generate_route(24, 6, 19)
        settings = dict(alignment=2, search_seed=101, restarts=3, perturbations=2,
                        perturb_regions=2, max_sweeps=8)
        first = search_route(route, 6, 2, **settings)
        second = search_route(route, 6, 2, **settings)
        self.assertEqual(len(first[2]), 9)
        self.assertEqual(first[1:], second[1:])
        self.assertTrue(np.array_equal(first[0].starts, second[0].starts))
        self.assertTrue(np.array_equal(first[0].counts, count_hits(route, first[0].starts, 6)))

    def test_unequal_lengths_and_bounded_exchanges(self):
        route = generate_route(24, 6, 93)
        widths = np.asarray((2, 4, 6, 8, 10, 12))
        initial = coordinate_ascent(route, widths, 2, (0, 2, 4, 6, 8, 10), alignment=2)
        exact = [sum(i in region[first:first+length]
                     for region, first, length in zip(route, initial.starts, widths)) for i in range(24)]
        self.assertEqual(initial.counts.tolist(), exact)
        settings = dict(alignment=2, seed=510, proposals=24, amounts=(2, 4, 8), max_sweeps=8)
        result, log = exchange_lengths(route, initial, 2, **settings)
        repeated, repeated_log = exchange_lengths(route, initial, 2, **settings)
        self.assertGreaterEqual(result.dimension, initial.dimension)
        self.assertEqual(int(result.widths.sum()), int(widths.sum()))
        self.assertTrue(np.all(result.widths >= 2))
        self.assertTrue(np.all(result.widths % 2 == 0))
        self.assertEqual(log, repeated_log)
        self.assertTrue(np.array_equal(result.starts, repeated.starts))
        self.assertTrue(np.array_equal(result.widths, repeated.widths))
        self.assertTrue(np.array_equal(result.counts, count_hits(route, result.starts, result.widths)))

    def test_shortened_dimensions_and_zero_exit_confinement(self):
        # Actual GRS[4,3] over GF5, six independently shortened rows,
        # arbitrary one-symbol linear state, and non-common coordinate supports.
        route = generate_route(6, 4, 877)
        starts = (0, 1, 2, 3)
        width, parity, prime = 3, 1, 5
        counts = count_hits(route, starts, width)
        row_bases = []
        source_dimension = 0
        for row in range(6):
            omitted = [j for j in range(4) if row not in route[j, starts[j]:starts[j]+width]]
            basis = kernel_basis_prime([[1, j, j*j % prime] for j in omitted], 3)
            self.assertEqual(len(basis), max(int(counts[row])-parity, 0))
            row_bases.append((source_dimension, basis))
            source_dimension += len(basis)
        self.assertEqual(source_dimension, dimension(counts, parity))
        self.assertGreater(source_dimension, 4)
        generator = np.random.Generator(np.random.PCG64(723))
        for _ in range(12):
            state = [0]*source_dimension
            exits, outputs, inside = [], [], []
            for region in range(4):
                for position, row in enumerate(route[region]):
                    x = [0]*source_dimension
                    offset, basis = row_bases[row]
                    for k, polynomial in enumerate(basis):
                        x[offset+k] = sum(c*pow(region, j, prime) for j, c in enumerate(polynomial)) % prime
                    a, c, d = (int(v) for v in generator.integers(0, prime, size=3))
                    outputs.append([(u+a*v) % prime for u, v in zip(x, state)])
                    state = [(d*v+c*u) % prime for u, v in zip(x, state)]
                    inside.append(starts[region] <= position < starts[region]+width)
                    if position == starts[region]+width-1:
                        exits.append(state[:])
            witness_basis = kernel_basis_prime(exits, source_dimension)
            self.assertGreaterEqual(len(witness_basis), source_dimension-4)
            witness = witness_basis[0]
            word = [sum(c*x for c, x in zip(form, witness)) % prime for form in outputs]
            self.assertTrue(any(word))
            self.assertTrue(all(is_inside or not value for is_inside, value in zip(inside, word)))

    def test_guards(self):
        route = generate_route(8, 3, 1)
        validate_route(route)
        with self.assertRaises(ValueError):
            validate_route(np.zeros_like(route))
        with self.assertRaises(ValueError):
            coordinate_ascent(route, 4, 1, (0, 1, 2), alignment=2)

    def test_frozen_large_witnesses_without_search(self):
        from window_dimension_witnesses import WITNESSES, verify
        reports = [verify(witness) for witness in WITNESSES]
        self.assertEqual([(x['route_seed'], x['dimension'], x['tail_excluded_dimension']) for x in reports],
                         [(1, 136, 136), (2, 133, 132)])
        self.assertTrue(all(x['surviving_dimension_lower'] > 0 for x in reports))


if __name__ == '__main__':
    unittest.main()
