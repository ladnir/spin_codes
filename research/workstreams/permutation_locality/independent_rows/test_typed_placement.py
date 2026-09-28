"""Exact labeled-slot enumeration and outward checks for two group types."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import comb
from pathlib import Path
import sys
import unittest

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

from typed_placement import typed_placement


def exact(operators, counts, epochs, slots, **kwargs):
    return typed_placement(operators,counts,epochs=epochs,slots=slots,
                           matrix=fmpq_mat,rounding=lambda value:value,**kwargs)


def identity(size):
    return fmpq_mat([[int(i==j) for j in range(size)] for i in range(size)])


def rational(value):
    value=value.fmpq() if isinstance(value,arb) else value
    return Q(int(value.p),int(value.q))


def enumerate_slots(operators, epochs, slots):
    size=next(iter(operators.values())).nrows()
    sums=defaultdict(lambda:fmpq_mat(size,size))
    numbers=defaultdict(int)
    for assignment in product((0,1,2),repeat=epochs*slots):
        value=identity(size)
        for e in range(epochs):
            local=assignment[e*slots:(e+1)*slots]
            value=value*operators[local.count(1),local.count(2)]
        counts=assignment.count(1),assignment.count(2)
        sums[counts]+=value
        numbers[counts]+=1
    result={}
    for (n1,n2),value in sums.items():
        denominator=comb(epochs*slots,n1)*comb(epochs*slots-n1,n2)
        assert numbers[n1,n2]==denominator
        result[n1,n2]=value/denominator
    return result


class TypedPlacementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec=192

    def test_exhaustive_noncommuting_slot_assignments(self):
        for slots in (1,2):
            operators={(a,b):fmpq_mat([[1+a,1+b],[a+b,2+a*b]])
                       for a in range(slots+1) for b in range(slots-a+1)}
            self.assertNotEqual(operators[1,0]*operators[0,1],operators[0,1]*operators[1,0])
            for epochs in (1,2,3):
                reference=enumerate_slots(operators,epochs,slots)
                for counts,value in reference.items():
                    self.assertEqual(exact(operators,counts,epochs,slots),value)

    def test_homogeneous_reduction(self):
        sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
        from occupancy_model import placement
        one_type=[fmpq_mat([[1,1],[0,1]]),fmpq_mat([[1,0],[1,1]]),
                  fmpq_mat([[2,1],[0,1]])]
        operators={(a,b):one_type[a+b] for a in range(3) for b in range(3-a)}
        for epochs in (1,2,3):
            reference=placement(one_type,epochs=epochs,windows=2,matrix=fmpq_mat,
                                rounding=lambda value:value,maximum_groups=epochs*2)
            for total in range(epochs*2+1):
                for first in range(total+1):
                    self.assertEqual(exact(operators,(first,total-first),epochs,2),reference[total])

    def test_grid_and_pruned_recurrence_agree(self):
        operators={(a,b):fmpq_mat([[1+a,fmpq(1,3+b)],[fmpq(1,2+a),2+b]])
                   for a in range(3) for b in range(3-a)}
        grid=exact(operators,(2,2),3,2,return_grid=True)
        reference=enumerate_slots(operators,3,2)
        self.assertEqual(set(grid),set(product(range(3),repeat=2)))
        for counts,value in grid.items():
            self.assertEqual(value,reference[counts])
            self.assertEqual(value,exact(operators,counts,3,2))

    def test_selected_zero_packets_are_not_unselected_slots(self):
        # Each selected type supplies a zero packet with its own probability.
        # Average those packet realizations inside the epoch operator first.
        # Type-count placement must retain the original selected groups.
        p1,p2=fmpq(1,3),fmpq(3,4)
        by_nonzero=[fmpq_mat([[1,1],[0,1]]),fmpq_mat([[1,0],[1,1]]),
                    fmpq_mat([[2,1],[0,1]])]
        operators={}
        for n1 in range(3):
            for n2 in range(3-n1):
                value=fmpq_mat(2,2)
                for k1 in range(n1+1):
                    for k2 in range(n2+1):
                        mass=(comb(n1,k1)*p1**k1*(1-p1)**(n1-k1)
                              *comb(n2,k2)*p2**k2*(1-p2)**(n2-k2))
                        value+=by_nonzero[k1+k2]*mass
                operators[n1,n2]=value
        reference=enumerate_slots(operators,2,2)
        for counts in ((1,1),(2,1),(1,2),(2,2)):
            self.assertEqual(exact(operators,counts,2,2),reference[counts])
            direct=fmpq_mat(2,2)
            placements=0
            for labels in product((0,1,2),repeat=4):
                if (labels.count(1),labels.count(2)) != counts:
                    continue
                placements+=1
                selected=[i for i,label in enumerate(labels) if label]
                for bits in product((0,1),repeat=len(selected)):
                    probability=fmpq(1)
                    nonzero=[0]*4
                    for position,bit in zip(selected,bits):
                        p=p1 if labels[position]==1 else p2
                        probability*=p if bit else 1-p
                        nonzero[position]=bit
                    direct+=(by_nonzero[sum(nonzero[:2])]
                             *by_nonzero[sum(nonzero[2:])])*probability
            self.assertEqual(exact(operators,counts,2,2),direct/placements)
        # If zero packets were silently treated as unselected positions,
        # these nontrivial typed operators would not be used correctly.
        self.assertNotEqual(operators[1,0],operators[0,0])
        self.assertNotEqual(operators[1,0],operators[0,1])

    def test_outward_enclosure(self):
        exact_operators={(a,b):fmpq_mat([[fmpq(1+a,7),fmpq(1+b,9)],
                                       [fmpq(a+b,11),fmpq(2+a*b,13)]])
                         for a in range(3) for b in range(3-a)}
        outward_operators={key:arb_mat([[arb(int(value[i,j].p))/int(value[i,j].q)
                                         for j in range(2)] for i in range(2)])
                            for key,value in exact_operators.items()}
        for counts in ((0,0),(1,0),(0,2),(2,2),(3,2)):
            reference=exact(exact_operators,counts,3,2)
            actual=typed_placement(outward_operators,counts,epochs=3,slots=2)
            for i in range(2):
                for j in range(2):
                    self.assertTrue(actual[i,j].is_exact())
                    self.assertGreaterEqual(rational(actual[i,j]),rational(reference[i,j]))

    def test_normalization_and_type_swap(self):
        ones={(a,b):identity(2) for a in range(5) for b in range(5-a)}
        self.assertEqual(exact(ones,(8,8),8,4),identity(2))
        operators={(a,b):fmpq_mat([[1+a,1+b],[a+b,2+a*b]])
                   for a in range(3) for b in range(3-a)}
        swapped={(b,a):value for (a,b),value in operators.items()}
        for counts in ((0,0),(1,2),(3,1),(2,4)):
            self.assertEqual(exact(operators,counts,3,2),
                             exact(swapped,counts[::-1],3,2))

    def test_forced_counts_sparse_operators_and_validation(self):
        operator=fmpq_mat([[1,1],[0,1]])
        self.assertEqual(exact({(2,0):operator},(6,0),3,2),operator**3)
        self.assertEqual(exact({(0,0):operator},(0,0),3,2),operator**3)
        self.assertEqual(exact({(0,0):operator},(0,0),0,2),identity(2))
        with self.assertRaises(ValueError):
            exact({(0,0):operator},(1,0),1,2)
        with self.assertRaises(ValueError):
            exact({(0,0):operator},(3,0),1,2)
        with self.assertRaises(ValueError):
            exact({(0,0):fmpq_mat([[-1]])},(0,0),1,2)
        with self.assertRaises(ValueError):
            exact({(0,0):operator},(-1,0),1,2)


if __name__=='__main__':
    unittest.main()
