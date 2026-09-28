"""Exact finite-support checks for interval density domination."""
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from itertools import permutations, product
from math import comb, log
from pathlib import Path
import sys
import unittest

from flint import ctx

from row_counts import (forced_reference, group_box_multiplicity,
                        packet_reference_factor, row_gamma_arb,
                        row_gamma_exact, row_gamma_function, row_gamma_log)


def rational(point):
    value = point.fmpq()
    return Q(int(value.p), int(value.q))


def exact_log(value):
    return log(value.numerator) - log(value.denominator)


class RowCountTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 192

    def test_exhaustive_support_domination(self):
        for n in range(1,6):
            caps = [Q((w + 1) if w % 2 == 0 else 0) for w in range(n + 1)]
            for lo in range(n + 1):
                for hi in range(lo, n + 1):
                    for exclude_zero in (False, True):
                        for p in (Q(1,5), Q(1,2), Q(4,5)):
                            gamma = row_gamma_exact(caps,lo,hi,p,exclude_zero=exclude_zero)
                            ratios = []
                            actual_moment = reference_moment = Q(0)
                            for mask in range(1 << n):
                                w = mask.bit_count()
                                actual = (caps[w] / comb(n,w)
                                          if lo <= w <= hi and (w or not exclude_zero) else Q(0))
                                reference = p**w * (1-p)**(n-w)
                                self.assertLessEqual(actual,gamma*reference)
                                ratios.append(actual/reference)
                                # A nonnegative function distinguishing masks of
                                # equal weight also obeys measure domination.
                                function = 1 + mask*mask
                                actual_moment += function*actual
                                reference_moment += function*reference
                            self.assertEqual(gamma,max(ratios))
                            self.assertLessEqual(actual_moment,gamma*reference_moment)

    def test_exact_outward_and_proposal(self):
        caps = [Q(1),Q(0),Q(7,3),Q(2),Q(0),Q(11,5)]
        for lo in range(6):
            for hi in range(lo,6):
                objective = row_gamma_function(caps,lo,hi,exclude_zero=True)
                for p in (Q(1,10**6),Q(1234567,10**9),Q(1,2),Q(999999,10**6)):
                    exact = row_gamma_exact(caps,lo,hi,p,exclude_zero=True)
                    outward = row_gamma_arb(caps,lo,hi,p,exclude_zero=True)
                    self.assertGreaterEqual(rational(outward),exact)
                    if exact:
                        self.assertAlmostEqual(objective(float(p)),exact_log(exact),places=9)
                    else:
                        self.assertEqual(objective(float(p)),float('-inf'))

    def test_zero_and_all_one_endpoints(self):
        caps = [1,0,3,0,1]
        self.assertEqual(forced_reference(caps,0,1),0)
        self.assertEqual(forced_reference(caps,3,4),1)
        self.assertIsNone(forced_reference(caps,0,1,exclude_zero=True))
        self.assertIsNone(forced_reference(caps,1,3))
        for lo,hi,p in ((0,0,Q(0)),(4,4,Q(1))):
            self.assertEqual(row_gamma_exact(caps,lo,hi,p),1)
            self.assertEqual(rational(row_gamma_arb(caps,lo,hi,p)),1)
            self.assertEqual(row_gamma_log(caps,lo,hi,float(p)),0)
        self.assertEqual(row_gamma_exact(caps,0,0,Q(1),exclude_zero=True),0)
        self.assertEqual(row_gamma_log(caps,0,0,1.,exclude_zero=True),float('-inf'))
        with self.assertRaises(ValueError):
            row_gamma_exact(caps,0,4,Q(0))
        with self.assertRaises(ValueError):
            row_gamma_exact(caps,0,4,Q(1))
        self.assertEqual(row_gamma_log(caps,0,4,0.),float('inf'))
        self.assertEqual(row_gamma_log(caps,0,4,1.),float('inf'))

    def test_interval_maximum_is_not_a_shell_sum(self):
        caps = [1,0,6,0,1]
        p = Q(1,2)
        gamma = row_gamma_exact(caps,0,4,p)
        point_factors = [row_gamma_exact(caps,w,w,p) for w in (0,2,4)]
        self.assertEqual(gamma,16)
        self.assertEqual(gamma,max(point_factors))
        self.assertLess(gamma,sum(point_factors))

    def test_group_label_multiplicity(self):
        for rows in range(1,6):
            for active in range(rows+1):
                for pattern in product(((1,2),(3,4),(5,5)),repeat=active):
                    labels = tuple(sorted(pattern)) + ((0,0),)*(rows-active)
                    expected = len(set(permutations(labels)))
                    self.assertEqual(group_box_multiplicity(pattern,rows=rows),expected)
                self.assertEqual(group_box_multiplicity([(1,2)]*active,rows=rows),comb(rows,active))
        with self.assertRaises(ValueError):
            group_box_multiplicity([(1,3),(2,4)])
        with self.assertRaises(ValueError):
            group_box_multiplicity([(0,0)])

    def test_group_product_for_all_label_assignments(self):
        # Two active rows and one inactive row. Sum over all row labels;
        # the construction includes a uniform row-label permutation here
        # solely to make every assignment share the same reference moment.
        caps,n,p = [1,0,3,0],3,Q(2,5)
        gamma = row_gamma_exact(caps,1,3,p,exclude_zero=True)
        multiplicity = group_box_multiplicity([(1,3)]*2,rows=3)
        actual = reference = Q(0)
        for words in product(range(1<<n),repeat=3):
            active = [w for w in words if w]
            if len(active)==2 and all(w.bit_count()==2 for w in active):
                actual += 1 + (words[0]^words[1]^words[2])**2
            # Reference: average over the inactive row's three positions.
            for zero_index in range(3):
                if words[zero_index]:
                    continue
                mass = Q(1,3)
                for i,word in enumerate(words):
                    if i!=zero_index:
                        weight=word.bit_count()
                        mass *= p**weight*(1-p)**(n-weight)
                reference += mass*(1+(words[0]^words[1]^words[2])**2)
        self.assertLessEqual(actual,multiplicity*gamma**2*reference)

    def test_common_packet_reference(self):
        for r in range(5):
            for reference_rows in range(r,5):
                for p in (Q(1,5),Q(1,2),Q(4,5)):
                    self.assertEqual(packet_reference_factor([p]*r,reference_rows,p),
                                     (1-p)**(r-reference_rows))
        for probabilities in ((Q(1,5),Q(4,5)),(Q(0),Q(1,3),Q(1)),(Q(1,2),)*4):
            p = Q(3,5)
            factor = packet_reference_factor(probabilities,4,p)
            source = [Q(0)]*5
            for bits in product((0,1),repeat=len(probabilities)):
                mass=Q(1)
                for bit,probability in zip(bits,probabilities):
                    mass *= probability if bit else 1-probability
                source[sum(bits)] += mass
            ratios=[]
            for mask in range(16):
                weight=mask.bit_count()
                original=source[weight]/comb(4,weight)
                reference=p**weight*(1-p)**(4-weight)
                self.assertLessEqual(original,factor*reference)
                ratios.append(original/reference)
            self.assertEqual(factor,max(ratios))
        self.assertEqual(packet_reference_factor([],0,Q(0)),1)
        self.assertEqual(packet_reference_factor([Q(1)]*4,4,Q(1)),1)
        with self.assertRaises(ValueError):
            packet_reference_factor([Q(1,2)]*4,3,Q(1,2))
        with self.assertRaises(ValueError):
            packet_reference_factor([Q(1,2)],4,Q(0))

    def test_authenticated_bch_caps(self):
        sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
        from bch_joint_support import authenticated_caps
        with redirect_stdout(StringIO()):
            caps=authenticated_caps()
        self.assertEqual(len(caps),257)
        self.assertEqual(caps[0],1)
        self.assertEqual(caps[256],1)
        self.assertTrue(all(caps[w]==0 for w in range(1,256,2)))
        self.assertEqual(row_gamma_exact(caps,1,37,Q(1,2)),0)
        self.assertEqual(row_gamma_exact(caps,256,256,Q(1)),1)
        self.assertEqual(row_gamma_exact(caps,0,0,Q(0)),1)
        for lo,hi in ((38,64),(66,128),(130,218),(220,255)):
            p=Q(3,5)
            expected=max((Q(caps[w])/(comb(256,w)*p**w*(1-p)**(256-w))
                          for w in range(lo,hi+1) if caps[w]),default=Q(0))
            self.assertEqual(row_gamma_exact(caps,lo,hi,p,exclude_zero=True),expected)
            self.assertGreaterEqual(rational(row_gamma_arb(caps,lo,hi,p,exclude_zero=True)),expected)


if __name__ == '__main__':
    unittest.main()
