from fractions import Fraction as Q
from math import factorial
from unittest.mock import patch
import unittest

from flint import arb,ctx

from lifted_mixture import affine_upper,outer_log,LiftedModel
from mixture_cover import cover,replay
from row_mixture import envelope,pair_components


class LiftedMixtureTests(unittest.TestCase):
    def test_three_coordinate_bounds_against_exact_compositions(self):
        ctx.prec=192;n=6;q_min=2
        features=[(Q(0),Q(0),Q(0)),(Q(1,2),Q(1,4),Q(1,2)),
                  (Q(1,3),Q(1,3),Q(1)),(Q(0),Q(1),Q(0))]
        active=[0,1,1,1];coefs=[Q(1),Q(3,2),Q(5,3),Q(2)];values=[Q(2),Q(1,3),Q(3,2),Q(1)]
        cell=(Q(1,8),Q(3,8),Q(1,8),Q(3,4),Q(1,4),Q(3,4));compositions=[];total=Q(0)
        for a in range(n+1):
            for b in range(n-a+1):
                for c in range(n-a-b+1):
                    counts=[a,b,c,n-a-b-c]
                    if n-a<q_min:continue
                    means=[sum(Q(k,n)*f[j] for k,f in zip(counts,features)) for j in range(3)]
                    if any(not lo<=m<=hi for m,lo,hi in zip(means,cell[::2],cell[1::2])):continue
                    compositions.append(counts);term=Q(factorial(n))
                    for k,w in zip(counts,coefs):term*=w**k/factorial(k)
                    total+=term
        self.assertTrue(compositions)
        for dual in ((0,0,0,0),(2,-3,4,1),(-5,7,-1,3)):
            outer=outer_log(coefs,features,active,n,q_min,cell,dual).exp()
            self.assertTrue(outer>=arb(total.numerator)/total.denominator)
            upper=affine_upper(values,features,active,Q(q_min,n),cell,dual)
            for counts in compositions:self.assertGreaterEqual(upper,sum(Q(k,n)*v for k,v in zip(counts,values)))

    def test_exact_separation_does_not_drop_feasible_compositions(self):
        components=pair_components(envelope([1,2,3,2,1],4,Q(1,4)))
        instance=LiftedModel(components,None,2048,[1,Q(1,4),Q(1,4)])
        impossible=(Q(14,100),Q(16,100),Q(19,100),Q(21,100),Q(95,100),Q(1))
        self.assertFalse(instance.base.empty(impossible[:4]))
        self.assertTrue(instance.empty(impossible))
        for i in range(1,len(components)):
            for j in (0,i,len(components)-1):
                mean=[(a+b)/2 for a,b in zip(instance.features[i],instance.features[j])]
                point=tuple(x for value in mean for x in (value,value))
                self.assertFalse(instance.empty(point),(i,j,point))

    def test_three_dimensional_partition_replays(self):
        class Toy:
            root=(Q(0),Q(1))*3
            def empty(self,cell):return cell[0]>=Q(1,2)
            def proposal(self,cell):return (-20 if cell[1]-cell[0]>Q(1,2) else -100),[1,0,0,0]
            def outward(self,cell,witness):return arb(2)**-100
        instance=Toy()
        with patch('builtins.print'):record=cover(instance,max_cells=10)
        self.assertFalse(record['unresolved']);self.assertEqual(replay(instance,record),arb(2)**-100)
        del record['leaves']['1']
        with self.assertRaises(ValueError):replay(instance,record)

    def test_integer_contraction_retains_counts_and_boundaries(self):
        components=pair_components(envelope([1,2,3,2,1],4,Q(1,4)))
        instance=LiftedModel(components,None,1,[1,Q(1,4),Q(1,4)])
        for i in range(len(components)):
            for j in (0,i,len(components)-1):
                for count in (1,127,2048,4095):
                    mean=[(count*a+(4096-count)*b)/4096 for a,b in zip(instance.features[i],instance.features[j])]
                    cell=tuple(y for x in mean for y in (max(Q(0),x-Q(1,20000)),min(Q(1),x+Q(1,20000))))
                    contracted=instance.integer_cell(cell)
                    self.assertIsNotNone(contracted)
                    for x,lo,hi in zip(mean,contracted[::2],contracted[1::2]):self.assertTrue(lo<=x<=hi)
        gap=(Q(0),Q(1),Q(0),Q(1),Q(1,32768),Q(1,16384))
        self.assertIsNone(instance.integer_cell(gap))
        self.assertTrue(instance.empty(gap))
        for j,minimum in enumerate(instance.minimum_positive):
            cell=list(instance.root);cell[2*j+1]=minimum/2
            contracted=instance.integer_cell(tuple(cell))
            self.assertEqual(contracted[2*j+1],0)
            cell[2*j]=minimum/4
            self.assertIsNone(instance.integer_cell(tuple(cell)))
            cell[2*j]=0;cell[2*j+1]=minimum
            self.assertEqual(instance.integer_cell(tuple(cell))[2*j+1],minimum)


if __name__=='__main__':unittest.main()
