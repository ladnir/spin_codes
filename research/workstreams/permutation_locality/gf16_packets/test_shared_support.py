import unittest
from unittest.mock import patch
from collections import Counter
from itertools import permutations, product
from fractions import Fraction as Q
from math import comb
import shared_support as shared
import shared_sparse
from scalar_cover import multiply
from joint_support import gf2_rank, span
from basis_lattice import bound
from bch_joint_support import rank_total


def exact_law(packets):
    """Enumerate shared coordinate permutations and independent labels."""
    counts = Counter()
    for order in permutations(range(len(packets))):
        routed = [packets[i] for i in order]
        active = [i for i, value in enumerate(routed) if value]
        for labels in product(range(1,16),repeat=len(active)):
            output = [0]*len(packets)
            for i, label in zip(active,labels):
                output[i] = multiply(label,routed[i])
            counts[tuple(output)] += 1
    total = sum(counts.values())
    return {key: Q(value,total) for key,value in counts.items()}


class SharedSupportTests(unittest.TestCase):
    def test_law_depends_only_on_union_support(self):
        # Repeated packets (rank one) and independent packet values induce
        # the same label/support law. Different support sizes do not.
        for left,right in (((1,1,0),(3,5,0)),((7,0,7),(0,9,2)),((1,0,0),(0,0,15))):
            law = exact_law(left)
            self.assertEqual(law,exact_law(right))
            u = sum(bool(x) for x in left)
            self.assertEqual(len(law),comb(3,u)*15**u)
            self.assertEqual(set(law.values()),{Q(1,comb(3,u)*15**u)})

    def test_rank_counts_and_cdf_on_small_code(self):
        words = span([0b011,0b110])
        spectrum = [sum(w.bit_count()==u for w in words) for u in range(4)]
        shells = [[0]*4 for _ in range(5)]
        for rows in product(words,repeat=4):
            union = rows[0]|rows[1]|rows[2]|rows[3]
            shells[gf2_rank(rows)][union.bit_count()] += 1
        for rank in (1,2):
            self.assertEqual(sum(shells[rank]),rank_total(2,4,rank))
            for u in range(4):
                self.assertLessEqual(sum(shells[rank][:u+1]),bound(spectrum,4,rank,u))
        self.assertEqual(shells[1],[15*c for c in [0]+spectrum[1:]])
        self.assertEqual(sum(map(sum,shells[1:])),4**4-1)

    def test_invalid_request_fails_before_expensive_work(self):
        for updates,thresholds,tilts,precision in (([1],[1],['.1'],128),([2],[-1],['.1'],128),
                                                   ([2],[1],['0'],128),([2],[1],['.1'],64)):
            with self.assertRaises(ValueError):
                shared.run(updates,thresholds,tilts,precision)

    def test_invalid_sparse_request(self):
        for change in (dict(threshold=-1),dict(threshold=1<<21),dict(tilts=[]),
                       dict(tilts=['0']),dict(probe_supports=[37]),dict(refined_counts='yes'),
                       dict(joint_return_through=5),dict(lazy_density_through=-1),
                       dict(lazy_density_through=True),
                       dict(coupled_counts=True),dict(coupled_counts='yes'),
                       dict(joint_counts=True),dict(joint_counts='yes'),
                       dict(wide_counts=True),dict(wide_counts='yes'),
                       dict(occupancies=[2,2]),dict(proposal_buffer_bits=-1),dict(proposal_buffer_bits=float('nan')),
                       dict(target_bits=0),dict(target_bits=-1),dict(target_bits=True),
                       dict(target_bits=20.0),dict(target_bits='20')):
            args=dict(updates=[4],occupancies=[2],threshold=209715,tilts=['.01'],
                      precision=256,max_splits=1,target_bits=48)
            args.update(change)
            with self.assertRaises(ValueError):shared_sparse.run(**args)

    def test_positive_integer_sparse_budgets_reach_count_reconstruction(self):
        # The margin is a caller-selected proof budget, not a construction
        # restriction. Preserve fresh counts and outward verification.
        for bits in (1,20,32,40,48):
            with patch.object(shared_sparse,'shared_counts',side_effect=RuntimeError('fresh count sentinel')) as counts:
                with self.assertRaisesRegex(RuntimeError,'fresh count sentinel'):
                    shared_sparse.run([2],[1],104857,['.01'],128,0,bits)
                counts.assert_called_once()


if __name__ == '__main__':
    unittest.main()
