"""Tiny exact polynomial checks; no regional provider or proof-grid run."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest


def add(left, right, scale=F(1)):
    out = dict(left)
    for powers, value in right.items():
        out[powers] = out.get(powers, F(0)) + scale*value
        if not out[powers]:
            del out[powers]
    return out


def mul(left, right):
    out = {}
    for x, a in left.items():
        for y, b in right.items():
            powers = tuple(i+j for i, j in zip(x, y))
            out[powers] = out.get(powers, F(0)) + a*b
    return {powers: value for powers, value in out.items() if value}


def one(n):
    return {(0,)*n: F(1)}


def zpoly(n, members, value, scale=3):
    out = one(n)
    for member in members:
        powers = tuple(int(i == member) for i in range(n))
        out = mul(out, add(one(n), {powers: scale*value}))
    return out


def fpoly(n, members, law, scale=3):
    out = {}
    for value, probability in law:
        out = add(out, zpoly(n, members, value, scale), probability)
    return out


def subsets(n):
    return [frozenset(i for i in range(n) if mask >> i & 1)
            for mask in range(1 << n)]


def poly_power(base, exponent, n):
    out = one(n)
    for _ in range(exponent):
        out = mul(out, base)
    return out


class GroupUncrossingTests(unittest.TestCase):
    def test_all_two_and_three_variable_pair_inequalities_and_identity(self):
        law = ((F(0), F(1, 4)), (F(1, 3), F(1, 2)), (F(2), F(1, 4)))
        for n in (2, 3):
            for left, right in product(subsets(n), repeat=2):
                common, u, v = left & right, left-right, right-left
                difference = add(mul(fpoly(n, left | right, law), fpoly(n, common, law)),
                                 mul(fpoly(n, left, law), fpoly(n, right, law)), -1)
                identity = {}
                for (x, px), (y, py) in product(law, repeat=2):
                    term = mul(zpoly(n, common, x), zpoly(n, common, y))
                    term = mul(term, add(zpoly(n, u, x), zpoly(n, u, y), -1))
                    term = mul(term, add(zpoly(n, v, x), zpoly(n, v, y), -1))
                    identity = add(identity, term, px*py/2)
                self.assertEqual(difference, identity)
                self.assertTrue(all(coefficient >= 0 for coefficient in difference.values()))

    def test_equal_degree_incidence_concentration(self):
        n, arity = 3, 2
        law = ((F(0), F(2, 5)), (F(2, 3), F(3, 5)))
        all_groups = frozenset(range(n))
        rhs = poly_power(fpoly(n, all_groups, law), arity, n)
        row_sets = list(combinations(range(n), arity))
        for assignment in product(row_sets, repeat=n):
            lhs = one(n)
            for region in range(n):
                incidence = frozenset(g for g in range(n) if region in assignment[g])
                lhs = mul(lhs, fpoly(n, incidence, law))
            difference = add(rhs, lhs, -1)
            self.assertTrue(all(value >= 0 for value in difference.values()))

    def test_two_atom_mixture_identity(self):
        arity, probability, atom, eta = 2, F(2, 5), F(3, 7), F(1, 9)
        law = ((F(0), 1-probability), (atom, probability))
        for n in (2, 3):
            expanded = {}
            for selected in subsets(n):
                weight = eta**len(selected)*(1-eta)**(n-len(selected))
                expanded = add(expanded, poly_power(fpoly(n, selected, law), arity, n), weight)
            compact = {}
            for k in range(arity+1):
                weight = comb(arity, k)*probability**k*(1-probability)**(arity-k)
                term = one(n)
                for group in range(n):
                    factor = add({}, one(n), 1-eta)
                    factor = add(factor, poly_power(zpoly(n, {group}, atom), k, n), eta)
                    term = mul(term, factor)
                compact = add(compact, term, weight)
            self.assertEqual(expanded, compact)

    def test_aggregate_floor_does_not_imply_each_group_floor(self):
        degrees = (128, 16)
        self.assertGreaterEqual(sum(degrees), 144)
        self.assertFalse(all(degree >= 72 for degree in degrees))


if __name__ == '__main__':
    unittest.main()
