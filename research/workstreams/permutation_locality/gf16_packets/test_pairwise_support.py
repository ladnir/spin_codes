import unittest
from collections import Counter
from fractions import Fraction as Q
from itertools import accumulate, combinations, permutations, product
from math import comb

import pairwise_support as pairwise
import pairwise_sparse
from bch_joint_support import rank_total
from joint_support import gf2_rank, span


def subset_law(n, a, b):
    counts = [0]*(n+1)
    for left in combinations(range(n), a):
        for right in combinations(range(n), b):
            counts[len(set(left) | set(right))] += 1
    return [Q(v, comb(n,a)*comb(n,b)) for v in counts]


class PairwiseSupportTests(unittest.TestCase):
    def test_union_transform_against_every_small_subset_pair(self):
        for n in range(7):
            for a in range(n+1):
                for b in range(n+1):
                    left = [int(j == a) for j in range(n+1)]
                    right = [int(j == b) for j in range(n+1)]
                    self.assertEqual(pairwise.union_shells(left, right), subset_law(n,a,b))

    def test_union_cdf_decreases_with_each_input_size(self):
        for n in range(1,7):
            laws = {(a,b): list(accumulate(subset_law(n,a,b)))
                    for a in range(n+1) for b in range(n+1)}
            for a in range(n):
                for b in range(n+1):
                    self.assertTrue(all(x >= y for x,y in zip(laws[a,b], laws[a+1,b])))

    def test_small_code_rescaling_and_expected_counts(self):
        words = span([0b011, 0b110])
        ranks = [[0]*4 for _ in range(4)]
        for rows in product(words, repeat=4):
            h = gf2_rank(rows)
            if h:
                ranks[h-1][(rows[0]|rows[1]|rows[2]|rows[3]).bit_count()] += 1
        ranks = [list(accumulate(row)) for row in ranks]
        pairs = pairwise.pair_cdf(ranks, k=2)
        pair_shells = [0]*4
        for a,b in product(words, repeat=2):
            pair_shells[(a|b).bit_count()] += 1
        self.assertEqual(pairs, list(accumulate(pair_shells)))
        expected = [Q(0)]*4
        for a,b,c,d in product(words, repeat=4):
            if a|b|c|d:
                law = subset_law(3,(a|b).bit_count(),(c|d).bit_count())
                expected = [v+p for v,p in zip(expected,law)]
        result = pairwise.combine_pair_cdfs(pairs, pairs)
        self.assertEqual(result, list(accumulate(expected)))
        # Move comparison mass earlier while preserving total and zero atom.
        upper = [pairs[0]]+[pairs[-1]]*3
        relaxed = pairwise.combine_pair_cdfs(upper, pairs)
        self.assertTrue(all(x >= y for x,y in zip(relaxed,result)))
        self.assertEqual(result[-1], len(words)**4-1)

    def test_conditional_support_uniformity_under_pair_shuffles(self):
        # Both pairs overlap internally, and the two permutations are
        # independent. Enumerate all routes, including collisions.
        for rows in ((1,3,2,6),(3,3,5,5),(0,0,2,2),(7,1,3,2)):
            by_size = {}
            for first,second in product(permutations(range(3)),repeat=2):
                mask = 0
                for i in range(3):
                    packet = sum(((word >> (first[i] if j < 2 else second[i])) & 1) << j
                                 for j,word in enumerate(rows))
                    if packet:
                        mask |= 1 << i
                by_size.setdefault(mask.bit_count(),Counter())[mask] += 1
            for u, counts in by_size.items():
                self.assertEqual(len(counts), comb(3,u))
                self.assertEqual(len(set(counts.values())), 1)

    def test_invalid_inputs(self):
        for left,right in (([],[]),([1],[1,0]),([-1],[1]),([Q(1)],[1])):
            with self.assertRaises(ValueError):
                pairwise.union_shells(left,right)
        for cdf in ([],[0,1],[1,0],[1,Q(2)]):
            with self.assertRaises(ValueError):
                pairwise.combine_pair_cdfs(cdf,cdf)
        with self.assertRaises(ValueError):
            pairwise.pair_cdf([[0,rank_total(2,4,h)] for h in (1,2)],k=2)

    def test_shell_caps_tighten_cdf_without_changing_total(self):
        actual=[1,1,7,16]
        original=[1,10,16,16]
        result=pairwise.tighten_cdf(original,[1,1,8,10])
        self.assertEqual(result,[1,2,10,16])
        self.assertTrue(all(a<=b<=c for a,b,c in zip(actual,result,original)))
        with self.assertRaises(ArithmeticError):
            pairwise.tighten_cdf(original,[1,0,0,0])

    def test_invalid_sparse_parameters(self):
        for change in (dict(occupancies=[]),dict(occupancies=[1,1]),dict(occupancies=[2049]),
                       dict(tilts=['0']),dict(threshold=-1),dict(updates=True),
                       dict(precision=64),dict(max_splits=-1),dict(target_bits=39),
                       dict(coupled_counts=True),dict(refined_counts='yes')):
            args = dict(occupancies=[1],tilts=['.01'])
            args.update(change)
            with self.assertRaises(ValueError):
                pairwise_sparse.run(**args)


if __name__ == '__main__':
    unittest.main()
