from fractions import Fraction as Q
from itertools import product
from math import factorial
import unittest

from flint import arb,ctx

from affine_mixture import AffineModel,vertices,tangent_weights,weight_at,dyadic_upper,coupled_outer_log,affine_maps,tangent_anchors,resolve_minimum_groups
from probe import aq
from polytope import ClippedBox
import model


class AffineMomentTests(unittest.TestCase):
    def test_occupancy_retarget_can_only_restrict_the_saved_domain(self):
        record={'minimum_groups':401}
        self.assertEqual(resolve_minimum_groups(None),512)
        self.assertEqual(resolve_minimum_groups(None,record),401)
        self.assertEqual(resolve_minimum_groups(409,record,retarget=True),409)
        for args,options in (((409,record),{}),((400,record),{'retarget':True}),
                             ((0,),{}),((4097,),{}),((True,),{})):
            with self.assertRaises(ValueError):resolve_minimum_groups(*args,**options)

    def test_tilt_polishing_retains_the_best_rational_witness(self):
        class Toy:
            clip_moments=False
            def plane_proposal(self,cell,tilt,lam,density,*,clipped):
                self_lam=lam
                if not isinstance(self_lam,Q):raise AssertionError('rational proposal required')
                return float((lam-Q(1,8))**2)-80,{'plane':{'tilt':str(lam)}}
        old=(-70,{'plane':{'tilt':'1/10','input_tilt':['1','1/4','1/2']}})
        result=AffineModel.polish_proposal(Toy(),(),old)
        self.assertLess(result[0],-79.99)
        self.assertLess(abs(Q(result[1]['plane']['tilt'])-Q(1,8)),Q(1,10000))
        retained=(-1000,old[1])
        self.assertEqual(AffineModel.polish_proposal(Toy(),(),retained),retained)

    def test_clipped_vertex_plane_bounds_the_positive_moment(self):
        ctx.prec=192
        domain=ClippedBox([(0,0,0),(1,0,0),(0,1,0),(0,0,1)])
        corners=domain.vertices((0,Q(3,4))*3)
        maps=[(Q(1),(Q(-1,4),Q(-1,4),Q(0))),
              (Q(1,8),(Q(1),Q(0),Q(1,4))),
              (Q(1,16),(Q(0),Q(1),Q(1,8)))]
        anchors=[Q(1),Q(1,2),Q(1,2)];slope=[Q(1),Q(-2),Q(3)]
        def polynomial(w):return w[0]**4+2*w[0]*w[1]**2*w[2]+w[2]**4
        intercept=max(model.up(polynomial(tangent_weights(maps,anchors,p)).log()-aq(sum(s*x for s,x in zip(slope,p)))) for p in corners)
        for p in product([Q(k,8) for k in range(7)],repeat=3):
            if sum(p)>1:continue
            original=[aq(weight_at(f,p)) for f in maps]
            self.assertTrue(polynomial(original)<=(intercept+aq(sum(s*x for s,x in zip(slope,p)))).exp())

    def test_anchors_handle_a_nonpositive_center_without_overflow(self):
        cell=(Q(0),Q(1))*3
        maps=[(Q(-1),(Q(2),Q(0),Q(0))),(Q(-10),(Q(0),Q(21),Q(0)))]
        anchors=tangent_anchors(maps,cell)
        self.assertEqual(anchors,[Q(1,2),Q(11,2)])
        for corner in vertices(cell):
            for affine,anchor in zip(maps,anchors):
                self.assertLessEqual((weight_at(affine,corner)-anchor)/anchor,1)

    def test_logconvex_vertex_plane_bounds_a_positive_polynomial(self):
        ctx.prec=192
        cell=(Q(0),Q(1,4),Q(0),Q(1,4),Q(0),Q(1,4))
        maps=[(Q(1),(Q(-1),Q(-1),Q(0))),
              (Q(1,8),(Q(1),Q(0),Q(1,4))),
              (Q(1,16),(Q(0),Q(1),Q(-1,8)))]
        anchors=[weight_at(f,(Q(1,8),)*3) for f in maps];slope=[Q(1),Q(-2),Q(3)]
        def polynomial(w):return 3*w[0]**3+w[0]*w[1]*w[2]+2*w[1]**2*w[2]+w[2]**3
        offset=max(model.up(polynomial(tangent_weights(maps,anchors,p)).log()-aq(sum(s*x for s,x in zip(slope,p))))
                   for p in vertices(cell))
        for point in product([Q(i,16) for i in range(5)],repeat=3):
            original=[aq(weight_at(f,point)) for f in maps]
            tangent=tangent_weights(maps,anchors,point)
            self.assertTrue(all(t>=o for t,o in zip(tangent,original)))
            self.assertTrue(polynomial(original)<=(offset+aq(sum(s*x for s,x in zip(slope,point)))).exp())
        for value in tangent_weights(maps,anchors,(Q(1,4),)*3):
            self.assertTrue(aq(dyadic_upper(value))>=value)

    def test_affine_maps_bound_every_small_feasible_composition(self):
        features=[(Q(0),Q(0),Q(0)),(Q(1,2),Q(1,4),Q(1,2)),(Q(1,4),Q(1,2),Q(1))]
        active=[0,1,1];values=[[Q(1),Q(2),Q(3)],[Q(3),Q(1),Q(2)]]
        maps=affine_maps(values,features,active,Q(1,3),[[2,-3,4,1],[-5,7,-1,3]])
        for a in range(7):
            for b in range(7-a):
                counts=[a,b,6-a-b]
                if 6-a<2:continue
                mean=[sum(Q(k,6)*f[j] for k,f in zip(counts,features)) for j in range(3)]
                for affine,column in zip(maps,values):
                    self.assertGreaterEqual(weight_at(affine,mean),sum(Q(k,6)*v for k,v in zip(counts,column)))

    def test_coupled_outer_counts_the_complete_cell(self):
        ctx.prec=192;n=6;q_min=2
        features=[(Q(0),Q(0),Q(0)),(Q(1,2),Q(1,4),Q(1,2)),(Q(1,4),Q(1,2),Q(1))]
        active=[0,1,1];coefs=[Q(1),Q(3,2),Q(5,3)]
        cell=(Q(1,8),Q(3,8),Q(1,8),Q(3,4),Q(1,4),Q(3,4));slope=[Q(7),Q(-13),Q(11)]
        actual=arb(0);terms=0
        for a in range(n+1):
            for b in range(n-a+1):
                counts=[a,b,n-a-b]
                if n-a<q_min:continue
                mean=[sum(Q(k,n)*f[j] for k,f in zip(counts,features)) for j in range(3)]
                if any(not lo<=x<=hi for x,lo,hi in zip(mean,cell[::2],cell[1::2])):continue
                coefficient=Q(factorial(n));terms+=1
                for k,c in zip(counts,coefs):coefficient*=c**k/factorial(k)
                actual+=aq(coefficient)*aq(sum(s*x for s,x in zip(slope,mean))).exp()
        self.assertGreater(terms,0)
        for dual in ((0,0,0,0),(2,-3,4,1),(-5,7,-1,3)):
            bound=coupled_outer_log(coefs,features,active,n,q_min,cell,slope,dual).exp()
            self.assertTrue(bound>=actual)


if __name__=='__main__':unittest.main()
