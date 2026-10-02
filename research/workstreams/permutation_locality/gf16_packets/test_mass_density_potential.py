from fractions import Fraction as Q
from itertools import product
from math import comb
import unittest
import numpy as np
from flint import ctx,arb
import birth_classes
import occupancy_birth_classes
import fiber_density
import mass_density_potential as potential
from test_scalar_cover import endpoint


class SharedMassTests(unittest.TestCase):
    def test_budget_matches_exhaustive_integer_allocations(self):
        ctx.prec=192
        for caps in ((1,2,3),(0,3,1),(2,2,2)):
            for mass in range(8):
                for values in ((1,4,2),(0,1,1),(3,2,4)):
                    truth=max(sum(x*v for x,v in zip(xs,values)) for xs in product(*(range(c+1) for c in caps)) if sum(xs)<=mass)
                    actual=potential.shared_budget(arb(mass),list(map(arb,caps)),list(map(arb,values)))
                    self.assertEqual(endpoint(actual),truth)
        with self.assertRaises(ValueError):potential.shared_budget(arb(-1),[],[])
        with self.assertRaises(ValueError):potential.shared_budget(arb(1),[arb(1)],[])

    def test_nonlinear_steps_bound_exact_weighted_functions(self):
        ctx.prec=192;z=Q(3,4)
        images=[0x93*(a&1)^0x65*((a>>1)&1)^0x3c*((a>>2)&1) for a in range(8)]
        for columns in ([1,2,4,3,6,5,7,1],[0]*8,[1]*8):
            for updates in (2,3):
                data=fiber_density.attach(birth_classes.prepare(images,columns,3,updates))
                source=occupancy_birth_classes.outward_at_z(data,potential.aq(z))
                classes=fiber_density.candidate(data,source,potential.aq(z),1,2,allocation='classes',include_uniform=True)
                bound=potential.Bound(data,source,potential.aq(z),[classes])
                levels=list(data['birth_class_levels']);alpha=Q(1,2**updates)
                exact=[[[Q(0)]*8 for _ in range(8)] for _ in range(3)]
                for x in range(256):
                    j=int(bool(x&15))+int(bool(x>>4));chance=Q(1,comb(2,j)*15**j);feedback=0
                    for bit,c in enumerate(columns):
                        if x>>bit&1:feedback^=c
                    for a in range(8):
                        weight=chance*z**(images[a]^x).bit_count()
                        if not a:exact[j][a][feedback]+=weight
                        else:
                            exact[j][a][a^feedback]+=alpha*weight
                            for b in range(1,8):exact[j][a][b^feedback]+=(1-alpha)*weight/7
                for sequence in product(range(3),repeat=3):
                    values=[Q(1+a*a,7) for a in range(8)]
                    initial=[values[0],max(values[1:]),sum(values[1:])/7,
                        *(max(values[a] for a in range(1,8) if images[a].bit_count()==level) for level in levels)]
                    upper=list(map(potential.aq,initial))
                    for j in sequence:
                        numeric=bound.floating(j,np.array(list(map(float,upper))))
                        larger=bound.floating(j,np.array(list(map(float,upper)))*2)
                        np.testing.assert_allclose(larger,numeric*2,rtol=1e-12,atol=1e-13)
                        upper=bound.outward(j,upper)
                        np.testing.assert_allclose(numeric,list(map(float,upper)),rtol=1e-12,atol=1e-13)
                        values=[sum(c*v for c,v in zip(row,values)) for row in exact[j]]
                        self.assertGreaterEqual(endpoint(upper[0]),values[0])
                        self.assertGreaterEqual(endpoint(upper[1]),max(values[1:]))
                        self.assertGreaterEqual(endpoint(upper[2]),sum(values[1:])/7)
                        for i,level in enumerate(levels,3):
                            self.assertGreaterEqual(endpoint(upper[i]),max(values[a] for a in range(1,8) if images[a].bit_count()==level))
                vectors=np.array([[1,2,3,*(4+i for i in range(len(levels)))],
                                  [2,4,5,*(2+i for i in range(len(levels)))]],dtype=float)
                for j in range(3):
                    np.testing.assert_allclose(bound.floating(j,vectors),np.array([bound.floating(j,v) for v in vectors]))
                    for linear in bound.float_linear[:,j]:self.assertTrue((bound.floating(j,vectors)<=vectors@linear.T+1e-12).all())


if __name__=='__main__':unittest.main()
