import itertools
import math
import unittest
from fractions import Fraction as F

import numpy as np

import dense_occupation_ranges as dense
import typed_dense_boxes as typed
import activation_occupation as general
from test_activation_q1 import mul,identity
from composition_occupation import terminal_logs


def lattice(total,categories):
    return [c for c in itertools.product(range(total+1),repeat=categories) if sum(c)==total]


class TypedDenseTest(unittest.TestCase):
    def test_marked_moment_and_proposal_gradient(self):
        epoch=general.epoch_logs(8,4,{4:14,8:1},[1,0,0,0,14,0,0,0,1],.7,8)
        corners=typed.vertices([0,0,0],[16,16,16],16)
        logits=np.array([.2,-.4,.3]);ps=np.array([0.,.5,1.]);costs=np.array([0.,2.,0.])
        value,gradient=typed.proposal_objective(logits,corners,16,4,ps,costs,epoch,8)
        for j in range(3):
            plus=logits.copy();minus=logits.copy();plus[j]+=1e-5;minus[j]-=1e-5
            a=typed.proposal_objective(plus,corners,16,4,ps,costs,epoch,8)[0]
            b=typed.proposal_objective(minus,corners,16,4,ps,costs,epoch,8)[0]
            self.assertAlmostEqual(gradient[j],(a-b)/2e-5,places=7)
        self.assertAlmostEqual(float(gradient.sum()),0.,places=12)

    def test_vertices_bound_every_integer_count_and_splits_preserve_coverage(self):
        points=np.array(lattice(7,3))
        for lower,upper in (([0,0,0],[7,7,7]),([0,1,0],[5,6,4]),([2,2,2],[3,3,3])):
            lo=np.array(lower);hi=np.array(upper)
            feasible=points[np.all((points>=lo)&(points<=hi),axis=1)]
            corners=typed.vertices(lo,hi,7)
            self.assertTrue(all(tuple(c) in map(tuple,feasible) for c in corners))
            self.assertLessEqual(len(feasible),math.exp(typed.lattice_log_count(lo,hi))+1e-10)
            proposal=np.array([.2,.3,.5]);gamma=np.array([0.,1.,2.])
            actual=typed.point_logs(feasible,7,8,gamma,proposal,-3.,4,.2)
            bound=max(typed.point_logs(corners,7,8,gamma,proposal,-3.,4,.2))
            self.assertLessEqual(float(max(actual)),bound+1e-10)
            children=typed.split_box(lo,hi,corners)
            for point in feasible:
                self.assertEqual(sum(np.all(point>=a) and np.all(point<=b) for a,b in children),1)

    def test_multivariate_coefficient_bound_against_typed_arrangements(self):
        # Three categories: zero, Bernoulli(1/2), all-one. Two 2-bit epochs.
        p=[F(0),F(1,2),F(1)];counts=(2,1,1);length=4;t=2;block=4
        epoch=[[[F(1+i+j+d,24) for j in range(3)] for i in range(3)] for d in range(3)]
        assignments=set(itertools.permutations([0,0,1,2]))
        region=[[F(0)]*3 for _ in range(3)]
        for labels in assignments:
            product=identity(3)
            for start in (0,2):
                a,b=p[labels[start]],p[labels[start+1]]
                weights=[(1-a)*(1-b),a*(1-b)+(1-a)*b,a*b]
                matrix=[[sum(weights[d]*epoch[d][i][j] for d in range(3)) for j in range(3)] for i in range(3)]
                product=mul(product,matrix)
            region=[[x+y/len(assignments) for x,y in zip(row,other)] for row,other in zip(region,product)]
        value=identity(3)
        for _ in range(block):value=mul(value,region)
        gamma=[F(1),F(16,3),F(1)]
        exact=len(assignments)*math.prod(g**c for g,c in zip(gamma,counts))*sum(value[0])
        for proposal in ([.3,.3,.4],[.5,.25,.25]):
            theta=float(np.array(proposal)@np.array(p,dtype=float))
            matrix=dense.epoch_mixture_logs(np.log(np.array(epoch,dtype=float)),t,np.array([theta]))
            moment=float(terminal_logs(matrix,block*length//t)[0])
            result=typed.point_logs([counts],length,block,np.log(np.array(gamma,dtype=float)),proposal,moment,0,.2)[0]
            self.assertGreaterEqual(math.exp(result)+1e-12,float(exact))

    def test_adaptive_selected_boxes_cover_the_whole_dense_domain(self):
        engine=typed.TypedDense({2:2,4:1},4,8,4,{4:14,8:1},[1,0,0,0,14,0,0,0,1],16,[0.,.8])
        result=engine.search(minimum=3,maximum_nodes=31,target_bits=1000)
        points=[c for c in lattice(16,3) if c[0]<=13]
        boxes=result['selected_boxes']
        for point in points:
            self.assertEqual(sum(all(a<=x<=b for a,x,b in zip(row['lower'],point,row['upper'])) for row in boxes),1)
        union=float(np.logaddexp.reduce([row['own_log_bound'] for row in boxes]))
        self.assertAlmostEqual(union,result['log_union_upper'],places=10)


if __name__=='__main__':unittest.main()
