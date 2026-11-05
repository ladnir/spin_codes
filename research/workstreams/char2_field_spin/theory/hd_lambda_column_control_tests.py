"""Tiny score and category-adapter overlap checks, no production gates."""
from collections import Counter
from itertools import product
import unittest
import numpy as np

from hd_lambda_column_control import scored_weights
from hd_lambda_scheduled_relations import column_branches
from hd_scheduled_column_control import column_branches as old_branches


class LambdaColumnTests(unittest.TestCase):
    def test_score_is_original_mask_linear(self):
        np.testing.assert_array_equal(scored_weights(tuple(range(4)),(2,1)),[0,2,1,3])
        for score in ((1,1),(2,1),(1,2),(2,1,1),(1,1,1,1)):
            weights=scored_weights(tuple(range(1<<len(score))),score)
            for masks in product(range(1<<len(score)),repeat=3):
                expected=sum(c*sum(m>>i&1 for m in masks) for i,c in enumerate(score))
                self.assertEqual(int(sum(weights[m] for m in masks)),expected)

    def test_small_exact_branches_agree_after_merging(self):
        def merge(branches):
            result=Counter()
            for w,m,n in branches:result[w,m]+=n
            return result
        for source in ((0,),(13,),(0,7),(3,15)):
            for j in (0,2,7,15):
                for schedule in ('static','input-inverse','input-inverse-output-xor'):
                    self.assertEqual(merge(column_branches(source,schedule,j)),
                                     merge(old_branches(source,schedule,j,'lambda')))

    def test_scores_bounded_and_positive(self):
        for score in ((),(0,1),(3,2),(1.,1),(1,1,1,1,1)):
            with self.assertRaises(ValueError):scored_weights((0,),score)


if __name__=='__main__':unittest.main()
