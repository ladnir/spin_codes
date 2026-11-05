"""Tiny exhaustive indexed-word and exact recurrence checks; no K20 gate."""

from collections import Counter
from itertools import product
from math import factorial, prod
import unittest

import hd_field128_gamma_exact_profile_counts as exact


def profile(bits):
    pairs = len(bits)//2
    return (sum(bits), sum(bits[2*i]*bits[2*i+1] for i in range(pairs)),
            sum(bits[2*i+1]*bits[(2*i+2) % len(bits)] for i in range(pairs)))


def tree_oracle(pairs, weight, first, second, u):
    v = weight-2*first-u
    edges = {(0, 2): pairs-weight+first, (0, 3): v, (1, 2): u, (1, 3): first,
             (2, 0): pairs-weight+second, (2, 1): first+u-second,
             (3, 0): first+v-second, (3, 1): second}
    if min(edges.values()) < 0:
        return 0
    degrees = [sum(n for (a, b), n in edges.items() if a == vertex) for vertex in range(4)]
    active = [v for v in range(4) if degrees[v]]
    root = next(v for v in (0, 1) if degrees[v])
    other = [v for v in active if v != root]
    choices = [[(b, n) for (a, b), n in edges.items() if a == v and n] for v in other]
    tau = 0
    for selected in product(*choices):
        following = dict(zip(other, (edge[0] for edge in selected)))
        good = True
        for vertex in other:
            seen = set()
            while vertex != root:
                if vertex in seen:
                    good = False
                    break
                seen.add(vertex)
                vertex = following[vertex]
            if not good:
                break
        if good:
            tau += prod(edge[1] for edge in selected)
    value, remainder = divmod(pairs*tau*prod(factorial(degrees[v]-1) for v in active),
                               prod(factorial(n) for n in edges.values()))
    if remainder:
        raise AssertionError('Nonintegral independent tree oracle')
    return value


class ExactProfileCountTests(unittest.TestCase):
    def test_every_small_indexed_word_and_fixed_flow(self):
        for pairs in range(1, 8):
            actual, flows = Counter(), Counter()
            for word in range(1 << (2*pairs)):
                bits = tuple((word >> i) & 1 for i in range(2*pairs))
                key = profile(bits)
                u = sum(bits[2*i]*(1-bits[2*i+1]) for i in range(pairs))
                actual[key] += 1
                flows[key+(u,)] += 1
            for weight in range(2*pairs+1):
                for first in range(max(0, weight-pairs), weight//2+1):
                    for second in range(max(0, weight-pairs), weight//2+1):
                        key = pairs, weight, first, second
                        self.assertEqual(exact.joint_count(*key), actual[weight, first, second])
                        for u in range(weight-2*first+1):
                            value = exact.fixed_flow_count(*key, u)
                            self.assertEqual(value, flows[weight, first, second, u])
                            self.assertEqual(value, tree_oracle(*key, u))

    def test_recurrence_reflection_and_direct_initialization(self):
        recurrence_steps = 0
        for pairs in range(1, 12):
            for weight in range(2*pairs+1):
                for first in range(max(0, weight-pairs), weight//2+1):
                    for second in range(max(0, weight-pairs), weight//2+1):
                        value, diagnostics, terms = exact._evaluate(pairs, weight, first, second, record=True)
                        recurrence_steps += diagnostics[2]
                        self.assertEqual(value, exact.joint_count(pairs, weight, first, second))
                        singles = weight-2*first
                        self.assertEqual(value, sum(exact.fixed_flow_count(pairs, weight, first, second, u)
                                                    for u in range(singles+1)))
                        for u, multiplicity, term in terms:
                            self.assertEqual(term, exact.fixed_flow_count(pairs, weight, first, second, u))
                            self.assertEqual(term, exact.fixed_flow_count(pairs, weight, first, second, singles-u))
                            self.assertEqual(multiplicity, 1 if 2*u == singles else 2)
        self.assertGreater(recurrence_steps, 100)

    def test_indexed_permutation_swap_and_complement_symmetries(self):
        for pairs in range(1, 7):
            for word in range(1 << (2*pairs)):
                bits = tuple((word >> i) & 1 for i in range(2*pairs))
                weight, first, second = profile(bits)
                # Even translation preserves both matchings; odd translation exchanges them.
                self.assertEqual(profile(bits[2:]+bits[:2]), (weight, first, second))
                self.assertEqual(profile(bits[1:]+bits[:1]), (weight, second, first))
                reflected = tuple(bits[(1-i) % len(bits)] for i in range(len(bits)))
                self.assertEqual(profile(reflected), (weight, first, second))
                complement = (2*pairs-weight, pairs-weight+first, pairs-weight+second)
                self.assertEqual(profile(tuple(1-x for x in bits)), complement)
                count = exact.joint_count(pairs, weight, first, second)
                self.assertEqual(count, exact.joint_count(pairs, weight, second, first))
                self.assertEqual(count, exact.joint_count(pairs, *complement))

    def test_constants_disconnected_and_missing_phase_vertices(self):
        for pairs in range(1, 20):
            self.assertEqual(exact.joint_count(pairs, 0, 0, 0), 1)
            self.assertEqual(exact.joint_count(pairs, 2*pairs, pairs, pairs), 1)
            self.assertEqual(exact.joint_count(pairs, pairs, 0, 0), 2)
            for u in range(1, pairs):
                self.assertEqual(exact.fixed_flow_count(pairs, pairs, 0, 0, u), 0)
            for full in range(1, pairs):
                self.assertEqual(exact.joint_count(pairs, 2*full, full, full), 0)
        self.assertEqual(exact.fixed_flow_count(7, 3, 0, 0, 0), 35)
        self.assertEqual(exact.fixed_flow_count(7, 3, 0, 0, 3), 35)

    def test_guards_cache_and_source_scope(self):
        for args in ((0, 0, 0, 0), (4097, 0, 0, 0), (4, -1, 0, 0),
                     (4, 9, 0, 0), (True, 0, 0, 0), (4, 2, 1., 0)):
            with self.assertRaises(ValueError):
                exact.joint_count(*args)
        self.assertEqual(exact.joint_count(4, 2, -1, 0), 0)
        self.assertEqual(exact.joint_count(4, 2, 2, 0), 0)
        self.assertEqual(exact.fixed_flow_count(4, 2, 0, 0, -1), 0)
        with self.assertRaises(ValueError):
            exact.fixed_flow_count(4, 2, 0, 0, True)
        exact._cached.cache_clear()
        first = exact.count_diagnostics(11, 8, 2, 3)
        second = exact.count_diagnostics(11, 14, 5, 6)
        self.assertEqual(first, second)
        self.assertGreater(exact._cached.cache_info().hits, 0)
        self.assertEqual(exact._cached.cache_info().maxsize, 4096)
        pins = exact.source_pins()
        self.assertEqual(set(pins), {'hd_field128_gamma_exact_profile_counts.py',
                                    '../HD_FIELD128_EXACT_PROFILE_COUNTS.md'})


if __name__ == '__main__':
    unittest.main()
