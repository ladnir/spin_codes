from fractions import Fraction as Q
import unittest
import numpy as np
from flint import arb,ctx

import kernel
from mixture import group_components,G
from atlas import Model,resolve_minimum_groups


class AtlasTests(unittest.TestCase):
    def test_domain_retarget_may_only_restrict_coverage(self):
        record={'minimum_groups':59}
        self.assertEqual(resolve_minimum_groups(None),59)
        self.assertEqual(resolve_minimum_groups(None,record),59)
        self.assertEqual(resolve_minimum_groups(65,record,retarget=True),65)
        for q in (58,0,2049,True):
            with self.assertRaises(ValueError):resolve_minimum_groups(q,record,retarget=True)
        with self.assertRaises(ValueError):resolve_minimum_groups(65,record)

    @classmethod
    def setUpClass(cls):
        rows=[(Q(1),Q(0)),(Q(1),Q(1)),(Q(2),Q(1,2)),(Q(1),Q(2,5)),(Q(1),Q(3,5))]
        cls.components=group_components(rows)
        cls.model=Model(cls.components,None,59,10485,[1,Q(1,32),Q(1,64),Q(1,128),Q(1,256)])

    def test_alternate_majorants_on_feasible_compositions(self):
        m=self.model;names=[r[0] for r in self.components]
        for pattern in [('0003',59,'1111',8),('2222',600,'0033',50),('0002',59,'0004',9)]:
            counts=[0]*len(names);counts[names.index(pattern[0])]=pattern[1];counts[names.index(pattern[2])]+=pattern[3]
            counts[0]=G-sum(counts)
            point=[sum(Q(n,G)*f[j] for n,f in zip(counts,m.features)) for j in range(1,5)]
            cell=tuple(v for x in point for v in (max(Q(0),x-Q(1,65536)),min(Q(1),x+Q(1,65536))))
            self.assertFalse(m.empty(cell))
            tilt=(1,Q(1,16),Q(1,4),Q(1,2),1)
            upper,_,_,duals=m.alternate(cell,tilt)
            replay,_,_,_=m.alternate(cell,tilt,duals)
            self.assertEqual(upper,replay)
            maps,_,_,_=m.affine_weights(cell,tilt,duals)
            affine=m.mapped_weights(maps,point)
            for j in range(5):
                exact=sum(Q(n,G)*r[2][j]/sum(p*t for p,t in zip(r[2],tilt)) for n,r in zip(counts,self.components))
                self.assertLessEqual(exact,upper[j])
                self.assertLessEqual(exact,affine[j])
                self.assertLessEqual(affine[j],upper[j])

    def test_alternate_plane_replays_and_validates_duals(self):
        ctx.prec=192
        data=kernel.prepare(list(range(8)),[1,2,4,3]*32,3)
        m=Model(self.components,data,59,10485,self.model.tilt)
        index=[r[0] for r in self.components].index('2222')
        point=[Q(600,G)*x for x in m.features[index][1:]]
        cell=tuple(v for x in point for v in (x-Q(1,65536),x+Q(1,65536)))
        tilt=[1,Q(1,16),Q(1,4),Q(1,2),1]
        legacy_score,legacy=m.plane_proposal(cell,Q(1,8))
        self.assertNotIn('input_tilt',legacy['plane'])
        self.assertAlmostEqual(legacy_score,float(m.outward(cell,legacy).log()/arb(2).log()),places=5)
        for density in (None,(2000,20,20,4,4)):
            score,witness=m.plane_proposal(cell,Q(1,8),density,input_tilt=tilt)
            upper=m.outward(cell,witness)
            self.assertTrue(upper>0)
            self.assertAlmostEqual(score,float(upper.log()/arb(2).log()),places=5)
        witness['plane']['weights_dual'][0][-1]='-1'
        with self.assertRaises(ValueError):m.outward(cell,witness)

    def test_verified_hull_contains_source_domain(self):
        m=self.model
        for f,a in zip(m.features,m.active):
            if not a:continue
            for factor in (Q(59,G),Q(1)):
                point=[factor*x for x in f[1:]]
                for eta,bound in m.facets:self.assertLessEqual(sum(e*x for e,x in zip(eta,point)),bound)

    def test_complete_composition_cell_replays_and_rejects_omission(self):
        ctx.prec=192
        data=kernel.prepare(list(range(8)),[1,2,4,3]*32,3)
        m=Model(self.components,data,65,10485,self.model.tilt)
        anchor=[r[0] for r in self.components].index('0003')
        first=Q(65,G)*m.features[anchor][1];epsilon=Q(1,2**32)
        cell=(first-epsilon,first+Q(1,64*G),Q(0),epsilon,Q(0),epsilon,Q(0),epsilon)
        complete=m.boundary(cell)
        self.assertEqual(len(complete),2)
        parameters=['1/8','1/8','1/4','1/2','1']
        rows=[dict(counts=list(c),parameters=parameters,posterior_shuffle=True) for c in complete]
        upper=m.outward(cell,dict(compositions=rows))
        self.assertTrue(upper>0)
        with self.assertRaises(ValueError):m.outward(cell,dict(compositions=rows[:-1]))
        with self.assertRaises(ValueError):m.outward(cell,dict(compositions=rows+[rows[0]]))


if __name__=='__main__':unittest.main()
